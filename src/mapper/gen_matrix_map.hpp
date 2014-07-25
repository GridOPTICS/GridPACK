/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/* *****************************************************************************
 * gen_matrix_map.hpp
 * gridpack
 * Bruce Palmer
 * July 1, 2014
 *
 * Documentation status
 * Revision status
 * File contents
 **************************************************************************** */

#ifndef GENMATRIXMAP_HPP_
#define GENMATRIXMAP_HPP_

#define NZ_PER_ROW

#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
#include "gridpack/parallel/parallel.hpp"
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/component/base_component.hpp>
#include <gridpack/network/base_network.hpp>
#include <gridpack/math/matrix.hpp>

//#define DBG_CHECK

namespace gridpack {
namespace mapper {

template <class _network>
class GenMatrixMap {
  public: 
/**
 * Initialize mapper for the given network and the current mode. Create global
 * arrays that contain offsets that will be used to create matrix from the
 * network component objects
 * @param network network that will generate matrix
 */
GenMatrixMap(boost::shared_ptr<_network> network)
  : p_network(network)
{
  p_row_Offsets = NULL;
  p_col_Offsets = NULL;
#ifdef NZ_PER_ROW
  p_nz_per_row = NULL;
#endif

  p_timer = NULL;
  //p_timer = gridpack::utility::CoarseTimer::instance();

  p_GAgrp = network->communicator().getGroup();
  p_me = GA_Pgroup_nodeid(p_GAgrp);
  p_nNodes = GA_Pgroup_nnodes(p_GAgrp);

  p_row_Offsets = new int[p_nNodes];
  p_col_Offsets = new int[p_nNodes];

  p_nBuses = p_network->numBuses();
  p_nBranches = p_network->numBranches();

  getDimensions();
  setOffsets();
  setIndices();
  numberNonZeros();
  GA_Pgroup_sync(p_GAgrp);
}

~GenMatrixMap()
{
  if (p_row_Offsets != NULL) delete [] p_row_Offsets;
  if (p_col_Offsets != NULL) delete [] p_col_Offsets;
#ifdef NZ_PER_ROW
  if (p_nz_per_row != NULL) delete [] p_nz_per_row;
#endif
  GA_Pgroup_sync(p_GAgrp);
}

/**
 * Generate matrix from current component state on network
 * @return return a pointer to new matrix
 */
boost::shared_ptr<gridpack::math::Matrix> mapToMatrix(void)
{
  gridpack::parallel::Communicator comm = p_network->communicator();
  int blockSize = p_maxRowIndex-p_minRowIndex+1;
  boost::shared_ptr<gridpack::math::Matrix>
    Ret(new gridpack::math::Matrix(comm, blockSize, blockSize, p_nz_per_row));
  loadBusData(*Ret,false);
  loadBranchData(*Ret,false);
  GA_Pgroup_sync(p_GAgrp);
  Ret->ready();
  return Ret;
}


/**
 * Reset existing matrix from current component state on network
 * @param matrix existing matrix (should be generated from same mapper)
 */
void mapToMatrix(gridpack::math::Matrix &matrix)
{
  int t_set, t_bus, t_branch;
  matrix.zero();
  loadBusData(matrix,false);
  loadBranchData(matrix,false);
  GA_Pgroup_sync(p_GAgrp);
  matrix.ready();
}

/**
 * Reset existing matrix from current component state on network
 * @param matrix existing matrix (should be generated from same mapper)
 */
void mapToMatrix(boost::shared_ptr<gridpack::math::Matrix> &matrix)
{
  mapToMatrix(*matrix);
}

/**
 * Overwrite elements of existing matrix. This can be used to overwrite selected
 * elements of a matrix
 * @param matrix existing matrix (should be generated from same mapper)
 */
void overwriteMatrix(gridpack::math::Matrix &matrix)
{
  loadBusData(matrix,false);
  loadBranchData(matrix,false);
  GA_Pgroup_sync(p_GAgrp);
  matrix.ready();
}

/**
 * Overwrite elements of existing matrix. This can be used to overwrite selected
 * elements of a matrix
 * @param matrix existing matrix (should be generated from same mapper)
 */
void overwriteMatrix(boost::shared_ptr<gridpack::math::Matrix> &matrix)
{
  overwriteMatrix(*matrix);
}

/**
 * Increment elements of existing matrix. This can be used to increment selected
 * elements of a matrix
 * @param matrix existing matrix (should be generated from same mapper)
 */
void incrementMatrix(gridpack::math::Matrix &matrix)
{
  loadBusData(matrix,true);
  loadBranchData(matrix,true);
  GA_Pgroup_sync(p_GAgrp);
  matrix.ready();
}

/**
 * Increment elements of existing matrix. This can be used to increment selected
 * elements of a matrix
 * @param matrix existing matrix (should be generated from same mapper)
 */
void incrementMatrix(boost::shared_ptr<gridpack::math::Matrix> &matrix)
{
  incrementMatrix(*matrix);
}

private:

/**
 * Check to see of both buses at either end of a branch belong to this processor
 * @return true if both buses belong to this processor, false if one does not
 */
bool isLocalBranch(int idx)
{
  int jdx1, jdx2;
  p_network->getBranchEndpoints(idx, &jdx1, &jdx2);
  bool check = true;
  check = check && p_network->getActiveBus(jdx1);
  check = check && p_network->getActiveBus(jdx2);
  return check;
}

/**
 * Evaluate dimensions of matrix and find out how many rows and columns are
 * contributed from each processor
 */
void getDimensions(void)
{
  int i, nval;
  // Find out how many rows and columns are contributed by this processor
  int nRows = 0;
  int nCols = 0;
  // Get number of rows and columns contributed by each bus
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      nval = p_network->getBus(i)->matrixNumRows();
      nRows += nval;
      nval = p_network->getBus(i)->matrixNumCols();
      nCols += nval;
    }
  }
  // Get number of rows and columns contributed by each branch
  for (i=0; i<p_nBranches; i++) {
    if (p_network->getActiveBranch(i)) {
      nval = p_network->getBranch(i)->matrixNumRows();
      nRows += nval;
      nval = p_network->getBranch(i)->matrixNumCols();
      nCols += nval;
    }
  }
  // Evaluate offsets for each processor
  int *sizebuf = new int[2*p_nNodes];
  for (i=0; i<2*p_nNodes; i++) {
    sizebuf[i] = 0;
  }
  sizebuf[2*p_me] = nRows;
  sizebuf[2*p_me+1] = nCols;
  GA_Pgroup_igop(p_GAgrp, sizebuf, 2*p_nNodes, "+");
  // Get total matrix dimensions and evaluate offsets for processor
  p_iDim = sizebuf[0];
  p_jDim = sizebuf[1];
  p_row_Offsets[0] = 0;
  p_col_Offsets[0] = 0;
  for (i=1; i<p_nNodes; i++) {
    p_iDim += sizebuf[2*i];
    p_jDim += sizebuf[2*i+1];
    p_row_Offsets[i] = p_row_Offsets[i-1] + sizebuf[2*(i-1)];
    p_col_Offsets[i] = p_col_Offsets[i-1] + sizebuf[2*(i-1)+1];
  }
  delete [] sizebuf;
}

/**
 * Evaluate offsets for each network component
 */
void setOffsets(void)
{
  // Interleave contributions from buses and branches to increase concentration
  // of matrix elements along the diagonal.
  int i,j,jdx,jdx1,jdx2;
  int *i_bus_offsets = new int[p_nBuses];
  int *i_branch_offsets = new int[p_nBranches];
  int *j_bus_offsets = new int[p_nBuses];
  int *j_branch_offsets = new int[p_nBranches];
  for (i=0; i<p_nBuses; i++) {
    i_bus_offsets[i] = 0;
    j_bus_offsets[i] = 0;
  }
  for (i=0; i<p_nBranches; i++) {
    i_branch_offsets[i] = 0;
    j_branch_offsets[i] = 0;
  }
  int icnt = 0;
  int jcnt = 0;
  int nsize;
  // Evaluate offsets for individual network components
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      i_bus_offsets[i] = icnt;
      icnt += p_network->getBus(i)->matrixNumRows();
      j_bus_offsets[i] = jcnt;
      jcnt += p_network->getBus(i)->matrixNumCols();
      std::vector<int> nghbrs = p_network->getConnectedBranches(i);
      nsize = nghbrs.size();
      for (j=0; j<nsize; j++) {
        // Need to avoid double counting of branches when evaluating offsets.
        // If branch is non-local and it is active, then include it in offsets.
        // Otherwise, if branch is local and bus i is equal to the "from" bus,
        // then include it in the offsets.
        jdx = nghbrs[j];
        if (isLocalBranch(jdx)) {
          p_network->getBranchEndpoints(jdx,&jdx1,&jdx2);
          if (jdx1 == i) {
            i_branch_offsets[jdx] = icnt;
            icnt += p_network->getBranch(jdx)->matrixNumRows();
            j_branch_offsets[jdx] = jcnt;
            jcnt += p_network->getBranch(jdx)->matrixNumCols();
          }
        } else {
          if (p_network->getActiveBranch(jdx)) {
            i_branch_offsets[jdx] = icnt;
            icnt += p_network->getBranch(i)->matrixNumRows();
            j_branch_offsets[jdx] = jcnt;
            jcnt += p_network->getBranch(i)->matrixNumCols();
          }
        }
      }
    }
  }
  // Total number of rows and columns from this processor have been evaluated,
  // now create buffers that can scatter individual offsets to global arrays
  int **i_bus_index = new int*[p_nBuses];
  int **j_bus_index = new int*[p_nBuses];
  int **i_branch_index = new int*[p_nBranches];
  int **j_branch_index = new int*[p_nBranches];
  int *i_bus_index_buf = new int[p_nBuses];
  int *j_bus_index_buf = new int[p_nBuses];
  int *i_branch_index_buf = new int[p_nBranches];
  int *j_branch_index_buf = new int[p_nBranches];
  int *i_bus_value_buf = new int[p_nBuses];
  int *j_bus_value_buf = new int[p_nBuses];
  int *i_branch_value_buf = new int[p_nBranches];
  int *j_branch_value_buf = new int[p_nBranches];
  int i_bus_cnt = 0;
  int j_bus_cnt = 0;
  int i_branch_cnt = 0;
  int j_branch_cnt = 0;
  int row_offset = p_row_Offsets[p_me];
  int col_offset = p_col_Offsets[p_me];
  int nbus = 0;
  int nbranch = 0;
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      nbus++;
      i_bus_value_buf[i_bus_cnt] = i_bus_offsets[i]+row_offset;
      i_bus_index_buf[i_bus_cnt] = p_network->getGlobalBusIndex(i);
      i_bus_index[i_bus_cnt] = &i_bus_index_buf[i_bus_cnt];
      i_bus_cnt++;

      j_bus_value_buf[j_bus_cnt] = j_bus_offsets[i]+col_offset;
      j_bus_index_buf[j_bus_cnt] = p_network->getGlobalBusIndex(i);
      j_bus_index[j_bus_cnt] = &j_bus_index_buf[j_bus_cnt];
      j_bus_cnt++;
    }
  }
  for (i=0; i<p_nBranches; i++) {
    if (p_network->getActiveBranch(i)) {
      nbranch++;
      i_branch_value_buf[i_branch_cnt] = i_branch_offsets[i]+row_offset;
      i_branch_index_buf[i_branch_cnt] = p_network->getGlobalBranchIndex(i);
      i_branch_index[i_branch_cnt] = &i_branch_index_buf[i_branch_cnt];
      i_branch_cnt++;

      j_branch_value_buf[j_branch_cnt] = j_branch_offsets[i]+col_offset;
      j_branch_index_buf[j_branch_cnt] = p_network->getGlobalBranchIndex(i);
      j_branch_index[j_branch_cnt] = &j_branch_index_buf[j_branch_cnt];
      j_branch_cnt++;
    }
  }
  delete [] i_bus_offsets;
  delete [] j_bus_offsets;
  delete [] i_branch_offsets;
  delete [] j_branch_offsets;
  // Create global arrays that hold column and row offsets for all buses and
  // branches in the network. First create map array for global arrays
  int *t_busMap = new int[p_nNodes];
  int *t_branchMap = new int[p_nNodes];
  for (i=0; i<p_nNodes; i++) {
    t_busMap[i] = 0;
    t_branchMap[i] = 0;
  }
  t_busMap[p_me] = nbus;
  t_branchMap[p_me] = nbranch;
  GA_Pgroup_igop(p_GAgrp, t_busMap, p_nNodes, "+");
  GA_Pgroup_igop(p_GAgrp, t_branchMap, p_nNodes, "+");
  int *busMap = new int[p_nNodes];
  int *branchMap = new int[p_nNodes];
  busMap[0] = 0;
  branchMap[0] = 0;
  int total_buses = t_busMap[0];
  int total_branches = t_branchMap[0];
  for (i=1; i<p_nNodes; i++) {
    busMap[i] = busMap[i-1] + t_busMap[i-1];
    total_buses += t_busMap[i];
    branchMap[i] = branchMap[i-1] + t_branchMap[i-1];
    total_branches += t_branchMap[i];
  }
  delete [] t_busMap;
  delete [] t_branchMap;

  int one = 1;
  g_bus_row_offsets = GA_Create_handle();
  GA_Set_data(g_bus_row_offsets, one, &total_buses, C_INT);
  GA_Set_irreg_distr(g_bus_row_offsets, busMap, &p_nNodes);
  GA_Set_pgroup(g_bus_row_offsets, p_GAgrp);
  if (!GA_Allocate(g_bus_row_offsets)) {
    // TODO: Some kind of error
  }
  GA_Zero(g_bus_row_offsets);

  g_bus_column_offsets = GA_Create_handle();
  GA_Set_data(g_bus_column_offsets, one, &total_buses, C_INT);
  GA_Set_irreg_distr(g_bus_column_offsets, busMap, &p_nNodes);
  GA_Set_pgroup(g_bus_column_offsets, p_GAgrp);
  if (!GA_Allocate(g_bus_column_offsets)) {
    // TODO: Some kind of error
  }
  GA_Zero(g_bus_column_offsets);

  g_branch_row_offsets = GA_Create_handle();
  GA_Set_data(g_branch_row_offsets, one, &total_branches, C_INT);
  GA_Set_irreg_distr(g_branch_row_offsets, branchMap, &p_nNodes);
  GA_Set_pgroup(g_branch_row_offsets, p_GAgrp);
  if (!GA_Allocate(g_branch_row_offsets)) {
    // TODO: Some kind of error
  }
  GA_Zero(g_branch_row_offsets);

  g_branch_column_offsets = GA_Create_handle();
  GA_Set_data(g_branch_column_offsets, one, &total_branches, C_INT);
  GA_Set_irreg_distr(g_branch_column_offsets, branchMap, &p_nNodes);
  GA_Set_pgroup(g_branch_column_offsets, p_GAgrp);
  if (!GA_Allocate(g_branch_column_offsets)) {
    // TODO: Some kind of error
  }
  GA_Zero(g_branch_column_offsets);

  delete [] busMap;
  delete [] branchMap;

  // Scatter offsets to global arrays
  NGA_Scatter(g_bus_row_offsets, i_bus_value_buf, i_bus_index, i_bus_cnt);
  NGA_Scatter(g_bus_column_offsets, j_bus_value_buf, j_bus_index, j_bus_cnt);
  NGA_Scatter(g_branch_row_offsets, i_branch_value_buf, i_branch_index, i_branch_cnt);
  NGA_Scatter(g_branch_column_offsets, j_branch_value_buf, j_branch_index, j_branch_cnt);

  delete [] i_bus_index;
  delete [] j_bus_index;
  delete [] i_branch_index;
  delete [] j_branch_index;

  delete [] i_bus_index_buf;
  delete [] j_bus_index_buf;
  delete [] i_branch_index_buf;
  delete [] j_branch_index_buf;
  delete [] i_bus_value_buf;
  delete [] j_bus_value_buf;
  delete [] i_branch_value_buf;
  delete [] j_branch_value_buf;
}

/**
 * Based on previously calculated offsets, set indices for a buses and branches.
 * It is up to the individual bus and branch implementations to store these
 * values.
 */
void setIndices(void)
{
  // Construct lists of indices that need to be collected
  int **bus_index = new int*[p_nBuses];
  int **branch_index = new int*[p_nBranches];
  int *bus_index_buf = new int[p_nBuses];
  int *branch_index_buf = new int[p_nBranches];
  int *i_bus_value_buf = new int[p_nBuses];
  int *j_bus_value_buf = new int[p_nBuses];
  int *i_branch_value_buf = new int[p_nBranches];
  int *j_branch_value_buf = new int[p_nBranches];
  int i, j;
  // Get offsets for all buses and branches;
  for (i=0; i<p_nBuses; i++) {
    bus_index_buf[i] = p_network->getGlobalBusIndex(i);
    bus_index[i] = &bus_index_buf[i];
  }
  for (i=0; i<p_nBranches; i++) {
    branch_index_buf[i] = p_network->getGlobalBranchIndex(i);
    branch_index[i] = &branch_index_buf[i];
  }
  NGA_Gather(g_bus_row_offsets, i_bus_value_buf, bus_index, p_nBuses);
  NGA_Gather(g_bus_column_offsets, j_bus_value_buf, bus_index, p_nBuses);
  NGA_Gather(g_branch_row_offsets, i_branch_value_buf, branch_index, p_nBranches);
  NGA_Gather(g_branch_column_offsets, j_branch_value_buf, branch_index, p_nBranches);

  // Offsets are now available. Set indices in all network components
  int offset, nrows, ncols, idx;
  for (i=0; i<p_nBuses; i++) {
    nrows = p_network->getBus(i)->matrixNumRows();
    if (nrows > 0) {
      offset = i_bus_value_buf[i];
      for (j=0; j<nrows; j++) {
        idx = offset+j;
        p_network->getBus(i)->matrixSetRowIndex(j,idx);
      }
    }
    ncols = p_network->getBus(i)->matrixNumCols();
    if (ncols > 0) {
      offset = j_bus_value_buf[i];
      for (j=0; j<ncols; j++) {
        idx = offset+j;
        p_network->getBus(i)->matrixSetColIndex(j,idx);
      }
    }
  }
  for (i=0; i<p_nBranches; i++) {
    nrows = p_network->getBranch(i)->matrixNumRows();
    if (nrows > 0) {
      offset = i_branch_value_buf[i];
      for (j=0; j<nrows; j++) {
        idx = offset+j;
        p_network->getBranch(i)->matrixSetRowIndex(j,idx);
      }
    }
    ncols = p_network->getBranch(i)->matrixNumCols();
    if (ncols > 0) {
      offset = j_branch_value_buf[i];
      for (j=0; j<ncols; j++) {
        idx = offset+j;
        p_network->getBranch(i)->matrixSetColIndex(j,idx);
      }
    }
  }

  delete [] bus_index;
  delete [] branch_index;

  delete [] bus_index_buf;
  delete [] branch_index_buf;
  delete [] i_bus_value_buf;
  delete [] j_bus_value_buf;
  delete [] i_branch_value_buf;
  delete [] j_branch_value_buf;

  // Global arrays are no longer needed so we can get rid of them
  GA_Destroy(g_bus_row_offsets);
  GA_Destroy(g_bus_column_offsets);
  GA_Destroy(g_branch_row_offsets);
  GA_Destroy(g_branch_column_offsets);
}

/**
 * Determine how many columns have non-zero values for each row in the matrix
 */
void numberNonZeros(void)
{
  // Evaluate max and min row indices for this processor
  p_minRowIndex = p_row_Offsets[p_me];
  if (p_me < p_nNodes-1) {
    p_maxRowIndex = p_row_Offsets[p_me+1] - 1;
  } else {
    p_maxRowIndex = p_iDim-1;
  }
  int dim = p_maxRowIndex - p_minRowIndex + 1;
  int *row_idx_buf = new int[dim];
  p_nz_per_row = new int[dim];
  int i, j;
  for (i=0; i<dim; i++) {
    row_idx_buf[i] = i+p_minRowIndex;
    p_nz_per_row[i] = 0;
  }
  delete [] row_idx_buf;
  int nvals;
  p_maxValues = 0;
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      nvals = p_network->getBus(i)->matrixNumValues();
      if (nvals > p_maxValues) p_maxValues = nvals;
      if (nvals > 0) {
        gridpack::ComplexType *values = new gridpack::ComplexType[nvals];
        int *rows = new int[nvals];
        int *cols = new int[nvals];
        p_network->getBus(i)->matrixGetValues(values, rows, cols);
        for (j=0; j<nvals; j++) {
          if (rows[j] >= p_minRowIndex && rows[j] <= p_maxRowIndex) {
            p_nz_per_row[rows[j]-p_minRowIndex]++;
          }
        }
        delete [] rows;
        delete [] cols;
        delete [] values;
      }
    }
  }
  for (i=0; i<p_nBranches; i++) {
    if (p_network->getActiveBranch(i)) {
      nvals = p_network->getBranch(i)->matrixNumValues();
      if (nvals > p_maxValues) p_maxValues = nvals;
      if (nvals > 0) {
        gridpack::ComplexType *values = new gridpack::ComplexType[nvals];
        int *rows = new int[nvals];
        int *cols = new int[nvals];
        p_network->getBranch(i)->matrixGetValues(values, rows, cols);
        for (j=0; j<nvals; j++) {
          if (rows[j] >= p_minRowIndex && rows[j] <= p_maxRowIndex) {
            p_nz_per_row[rows[j]-p_minRowIndex]++;
          }
        }
        delete [] rows;
        delete [] cols;
        delete [] values;
      }
    }
  }
}

/**
 * Add contributions from buses to matrix
 * @param matrix matrix to which contributions are added
 * @param flag flag to distinguish new matrix (true) from old (false)
 */
void loadBusData(gridpack::math::Matrix &matrix, bool flag)
{
  int i, j, nvals;
  ComplexType *values = new ComplexType[p_maxValues];
  int *rows = new int[p_maxValues];
  int *cols = new int[p_maxValues];
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      nvals = p_network->getBus(i)->matrixNumValues();
      p_network->getBus(i)->matrixGetValues(values,rows,cols);
      for (j=0; j<nvals; j++) {
        if (flag) {
          matrix.addElement(rows[j],cols[j],values[j]);
        } else {
          matrix.setElement(rows[j],cols[j],values[j]);
        }
      }
    }
  }
  delete [] values;
  delete [] rows;
  delete [] cols;
}

/**
 * Add contributions from branches to matrix
 * @param matrix matrix to which contributions are added
 * @param flag flag to distinguish new matrix (true) from old (false)
 */
void loadBranchData(gridpack::math::Matrix &matrix, bool flag)
{
  int i, j, nvals;
  ComplexType *values = new ComplexType[p_maxValues];
  int *rows = new int[p_maxValues];
  int *cols = new int[p_maxValues];
  for (i=0; i<p_nBranches; i++) {
    if (p_network->getActiveBranch(i)) {
      nvals = p_network->getBranch(i)->matrixNumValues();
      p_network->getBranch(i)->matrixGetValues(values,rows,cols);
      for (j=0; j<nvals; j++) {
        if (rows[j] >= p_minRowIndex && rows[j] <= p_maxRowIndex) {
          if (flag) {
            matrix.addElement(rows[j],cols[j],values[j]);
          } else {
            matrix.setElement(rows[j],cols[j],values[j]);
          }
        }
      }
    }
  }
  delete [] values;
  delete [] rows;
  delete [] cols;
}

    // Configuration information
int                         p_me;
int                         p_nNodes;

    // network information
boost::shared_ptr<_network> p_network;
int                         p_nBuses;
int                         p_nBranches;

    // matrix information
int                         p_iDim;
int                         p_jDim;
int                         p_minRowIndex;
int                         p_maxRowIndex;
int                         p_maxValues;
#ifdef NZ_PER_ROW
int*                        p_nz_per_row;
#endif

int*                        p_row_Offsets;
int*                        p_col_Offsets;

    // global matrix offset arrays
int                         g_bus_row_offsets;
int                         g_bus_column_offsets;
int                         g_branch_row_offsets;
int                         g_branch_column_offsets;
int                         p_GAgrp;

    // pointer to timer
gridpack::utility::CoarseTimer *p_timer;

};

} /* namespace mapper */
} /* namespace gridpack */

#endif //GENMATRIXMAP_HPP_
