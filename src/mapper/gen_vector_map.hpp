/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/* *****************************************************************************
 * gen_vector_map.hpp
 * gridpack
 * Bruce Palmer
 * July 23, 2014
 *
 * Documentation status
 * Revision status
 * File contents
 **************************************************************************** */

#ifndef GENVECTORMAP_HPP_
#define GENVECTORMAP_HPP_

#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
#include "gridpack/parallel/parallel.hpp"
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/component/base_component.hpp>
#include <gridpack/network/base_network.hpp>
#include <gridpack/math/vector.hpp>

//#define DBG_CHECK

namespace gridpack {
namespace mapper {

template <class _network>
class GenVectorMap {
  public: 
/**
 * Initialize mapper for the given network and the current mode. Create global
 * arrays that contain offsets that will be used to create vector from the
 * network component objects
 * @param network network that will generate vector
 */
GenVectorMap(boost::shared_ptr<_network> network)
  : p_network(network)
{
  p_Offsets = NULL;

  p_timer = NULL;
  //p_timer = gridpack::utility::CoarseTimer::instance();

  p_GAgrp = network->communicator().getGroup();
  p_me = GA_Pgroup_nodeid(p_GAgrp);
  p_nNodes = GA_Pgroup_nnodes(p_GAgrp);

  p_Offsets = new int[p_nNodes];

  p_nBuses = p_network->numBuses();
  p_nBranches = p_network->numBranches();

  getDimensions();
  setOffsets();
  setIndices();
  GA_Pgroup_sync(p_GAgrp);
}

~GenVectorMap()
{
  if (p_Offsets != NULL) delete [] p_Offsets;
  GA_Pgroup_sync(p_GAgrp);
}

/**
 * Generate vector from current component state on network
 * @return return a pointer to new vector
 */
boost::shared_ptr<gridpack::math::Vector> mapToVector(void)
{
  gridpack::parallel::Communicator comm = p_network->communicator();
  int blockSize = p_maxIndex-p_minIndex+1;
  boost::shared_ptr<gridpack::math::Vector>
    Ret(new gridpack::math::Vector(comm, blockSize));
  loadBusData(*Ret,false);
  loadBranchData(*Ret,false);
  GA_Pgroup_sync(p_GAgrp);
  Ret->ready();
  return Ret;
}


/**
 * Reset existing vector from current component state on network
 * @param vector existing vector (should be generated from same mapper)
 */
void mapToVector(gridpack::math::Vector &vector)
{
  int t_set, t_bus, t_branch;
  vector.zero();
  loadBusData(vector,false);
  loadBranchData(vector,false);
  GA_Pgroup_sync(p_GAgrp);
  vector.ready();
}

/**
 * Reset existing vector from current component state on network
 * @param vector existing vector (should be generated from same mapper)
 */
void mapToVector(boost::shared_ptr<gridpack::math::Vector> &vector)
{
  mapToVector(*vector);
}

/**
 * Push data from vector onto buses. Vector must be created with the
 * mapToVector method using the same GenVectorMap
 * @param vector vector containing data to be pushed to buses
 */
void mapToNetwork(const gridpack::math::Vector &vector)
{
  int i, j, nvals;
  ComplexType *values = new ComplexType[p_maxValues];
  int *idx = new int[p_maxValues];
  // get values from buses
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      nvals = p_network->getBus(i)->vectorNumElements();
      p_network->getBus(i)->vectorGetElementIndices(idx);
      for (j=0; j<nvals; j++) {
        vector.getElement(idx[j],values[j]);
      }
      p_network->getBus(i)->vectorSetElementValues(values);
    }
  }
  // get values from branches
  for (i=0; i<p_nBranches; i++) {
    if (p_network->getActiveBranch(i)) {
      nvals = p_network->getBranch(i)->vectorNumElements();
      p_network->getBranch(i)->vectorGetElementIndices(idx);
      for (j=0; j<nvals; j++) {
        vector.getElement(idx[j],values[j]);
      }
      p_network->getBranch(i)->vectorSetElementValues(values);
    }
  }
  delete [] values;
  delete [] idx;
  GA_Pgroup_sync(p_GAgrp);
}

/**
 * Push data from vector onto buses. Vector must be created with the
 * mapToVector method using the same GenVectorMap
 * @param vector vector containing data to be pushed to buses
 */
void mapToNetwork(boost::shared_ptr<gridpack::math::Vector> &vector)
{
  mapToNetwork(*vector);
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
 * Evaluate dimension of vector and find out how many elements are
 * contributed from each processor
 */
void getDimensions(void)
{
  int i, nval;
  // Find out how many rows and columns are contributed by this processor
  int nRows = 0;
  // Get number of rows and columns contributed by each bus
  p_maxValues = 0;
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      nval = p_network->getBus(i)->vectorNumElements();
      if (p_maxValues < nval) p_maxValues = nval;
      nRows += nval;
    }
  }
  // Get number of rows and columns contributed by each branch
  for (i=0; i<p_nBranches; i++) {
    if (p_network->getActiveBranch(i)) {
      nval = p_network->getBranch(i)->vectorNumElements();
      if (p_maxValues < nval) p_maxValues = nval;
      nRows += nval;
    }
  }
  // Evaluate offsets for each processor
  int *sizebuf = new int[p_nNodes];
  for (i=0; i<p_nNodes; i++) {
    sizebuf[i] = 0;
  }
  sizebuf[p_me] = nRows;
  GA_Pgroup_igop(p_GAgrp, sizebuf, p_nNodes, "+");
  // Get total vector dimension and evaluate offsets for processor
  p_Dim = sizebuf[0];
  p_Offsets[0] = 0;
  for (i=1; i<p_nNodes; i++) {
    p_Dim += sizebuf[i];
    p_Offsets[i] = p_Offsets[i-1] + sizebuf[i-1];
  }
  p_minIndex = p_Offsets[p_me];
  if (p_me < p_nNodes-1) {
    p_maxIndex = p_Offsets[p_me+1] - 1;
  } else {
    p_maxIndex = p_Dim - 1;
  }
  delete [] sizebuf;
}

/**
 * Evaluate offsets for each network component
 */
void setOffsets(void)
{
  // Interleave contributions from buses and branches to match matrices
  int i,j,jdx,jdx1,jdx2;
  int *i_bus_offsets = new int[p_nBuses];
  int *i_branch_offsets = new int[p_nBranches];
  for (i=0; i<p_nBuses; i++) {
    i_bus_offsets[i] = 0;
  }
  for (i=0; i<p_nBranches; i++) {
    i_branch_offsets[i] = 0;
  }
  int icnt = 0;
  int nsize;
  // Evaluate offsets for individual network components
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      i_bus_offsets[i] = icnt;
      icnt += p_network->getBus(i)->vectorNumElements();
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
            icnt += p_network->getBranch(jdx)->vectorNumElements();
          }
        } else {
          if (p_network->getActiveBranch(jdx)) {
            i_branch_offsets[jdx] = icnt;
            icnt += p_network->getBranch(i)->vectorNumElements();
          }
        }
      }
    }
  }
  // Total number of rows and columns from this processor have been evaluated,
  // now create buffers that can scatter individual offsets to global arrays
  int **i_bus_index = new int*[p_nBuses];
  int **i_branch_index = new int*[p_nBranches];
  int *i_bus_index_buf = new int[p_nBuses];
  int *i_branch_index_buf = new int[p_nBranches];
  int *i_bus_value_buf = new int[p_nBuses];
  int *i_branch_value_buf = new int[p_nBranches];
  int i_bus_cnt = 0;
  int i_branch_cnt = 0;
  int row_offset = p_Offsets[p_me];
  int nbus = 0;
  int nbranch = 0;
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      nbus++;
      i_bus_value_buf[i_bus_cnt] = i_bus_offsets[i]+row_offset;
      i_bus_index_buf[i_bus_cnt] = p_network->getGlobalBusIndex(i);
      i_bus_index[i_bus_cnt] = &i_bus_index_buf[i_bus_cnt];
      i_bus_cnt++;
    }
  }
  for (i=0; i<p_nBranches; i++) {
    if (p_network->getActiveBranch(i)) {
      nbranch++;
      i_branch_value_buf[i_branch_cnt] = i_branch_offsets[i]+row_offset;
      i_branch_index_buf[i_branch_cnt] = p_network->getGlobalBranchIndex(i);
      i_branch_index[i_branch_cnt] = &i_branch_index_buf[i_branch_cnt];
      i_branch_cnt++;
    }
  }
  delete [] i_bus_offsets;
  delete [] i_branch_offsets;
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
  g_bus_offsets = GA_Create_handle();
  GA_Set_data(g_bus_offsets, one, &total_buses, C_INT);
  GA_Set_irreg_distr(g_bus_offsets, busMap, &p_nNodes);
  GA_Set_pgroup(g_bus_offsets, p_GAgrp);
  if (!GA_Allocate(g_bus_offsets)) {
    // TODO: Some kind of error
  }
  GA_Zero(g_bus_offsets);

  g_branch_offsets = GA_Create_handle();
  GA_Set_data(g_branch_offsets, one, &total_branches, C_INT);
  GA_Set_irreg_distr(g_branch_offsets, branchMap, &p_nNodes);
  GA_Set_pgroup(g_branch_offsets, p_GAgrp);
  if (!GA_Allocate(g_branch_offsets)) {
    // TODO: Some kind of error
  }
  GA_Zero(g_branch_offsets);

  delete [] busMap;
  delete [] branchMap;

  // Scatter offsets to global arrays
  NGA_Scatter(g_bus_offsets, i_bus_value_buf, i_bus_index, i_bus_cnt);
  NGA_Scatter(g_branch_offsets, i_branch_value_buf, i_branch_index, i_branch_cnt);

  delete [] i_bus_index;
  delete [] i_branch_index;

  delete [] i_bus_index_buf;
  delete [] i_branch_index_buf;
  delete [] i_bus_value_buf;
  delete [] i_branch_value_buf;
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
  int *i_branch_value_buf = new int[p_nBranches];
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
  NGA_Gather(g_bus_offsets, i_bus_value_buf, bus_index, p_nBuses);
  NGA_Gather(g_branch_offsets, i_branch_value_buf, branch_index, p_nBranches);

  // Offsets are now available. Set indices in all network components
  int offset, nrows, ncols, idx;
  for (i=0; i<p_nBuses; i++) {
    nrows = p_network->getBus(i)->vectorNumElements();
    if (nrows > 0) {
      offset = i_bus_value_buf[i];
      for (j=0; j<nrows; j++) {
        idx = offset+j;
        p_network->getBus(i)->vectorSetElementIndex(j,idx);
      }
    }
  }
  for (i=0; i<p_nBranches; i++) {
    nrows = p_network->getBranch(i)->vectorNumElements();
    if (nrows > 0) {
      offset = i_branch_value_buf[i];
      for (j=0; j<nrows; j++) {
        idx = offset+j;
        p_network->getBranch(i)->vectorSetElementIndex(j,idx);
      }
    }
  }

  delete [] bus_index;
  delete [] branch_index;

  delete [] bus_index_buf;
  delete [] branch_index_buf;
  delete [] i_bus_value_buf;
  delete [] i_branch_value_buf;

  // Global arrays are no longer needed so we can get rid of them
  GA_Destroy(g_bus_offsets);
  GA_Destroy(g_branch_offsets);
}

/**
 * Add contributions from buses to vector
 * @param vector vector to which contributions are added
 * @param flag flag to distinguish new vector (true) from old (false)
 */
void loadBusData(gridpack::math::Vector &vector, bool flag)
{
  int i, j, nvals;
  ComplexType *values = new ComplexType[p_maxValues];
  int *idx = new int[p_maxValues];
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      nvals = p_network->getBus(i)->vectorNumElements();
      p_network->getBus(i)->vectorGetElementValues(values, idx);
      for (j=0; j<nvals; j++) {
//        if (flag) {
          vector.addElement(idx[j],values[j]);
//        } else {
//          vector.setElement(idx[j],values[j]);
//        }
      }
    }
  }
  delete [] values;
  delete [] idx;
}

/**
 * Add contributions from branches to vector
 * @param vector vector to which contributions are added
 * @param flag flag to distinguish new vector (true) from old (false)
 */
void loadBranchData(gridpack::math::Vector &vector, bool flag)
{
  int i, j, nvals;
  ComplexType *values = new ComplexType[p_maxValues];
  int *idx = new int[p_maxValues];
  for (i=0; i<p_nBranches; i++) {
    if (p_network->getActiveBranch(i)) {
      nvals = p_network->getBranch(i)->vectorNumElements();
      p_network->getBranch(i)->vectorGetElementValues(values,idx);
      for (j=0; j<nvals; j++) {
        if (idx[j] >= p_minIndex && idx[j] <= p_maxIndex) {
//          if (flag) {
            vector.addElement(idx[j],values[j]);
//          } else {
//            vector.setElement(idx[j],values[j]);
//          }
        }
      }
    }
  }
  delete [] values;
  delete [] idx;
}

    // Configuration information
int                         p_me;
int                         p_nNodes;

    // network information
boost::shared_ptr<_network> p_network;
int                         p_nBuses;
int                         p_nBranches;

    // vector information
int                         p_Dim;
int                         p_minIndex;
int                         p_maxIndex;
int                         p_maxValues;
#ifdef NZ_PER_ROW
int*                        p_nz_per_row;
#endif

int*                        p_Offsets;

    // global vector offset arrays
int                         g_bus_offsets;
int                         g_branch_offsets;
int                         p_GAgrp;

    // pointer to timer
gridpack::utility::CoarseTimer *p_timer;

};

} /* namespace mapper */
} /* namespace gridpack */

#endif //GENVECTORMAP_HPP_
