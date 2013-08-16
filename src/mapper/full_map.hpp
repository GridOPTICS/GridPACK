/* *****************************************************************************
 * full_map.hpp
 * gridpack
 * kglass, bjpalmer
 * Jul 22, 2013
 *
 * Documentation status
 * Revision status
 * File contents
 **************************************************************************** */

#ifndef FULLMATRIXMAP_HPP_
#define FULLMATRIXMAP_HPP_

#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
#include "gridpack/parallel/parallel.hpp"
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/component/base_component.hpp>
#include <gridpack/network/base_network.hpp>
#include <gridpack/math/matrix.hpp>

namespace gridpack {
namespace mapper {

template <class _network>
class FullMatrixMap {
  public: 
/**
 * Initialize mapper for the given network and the current mode. Create global
 * arrays that contain offsets that will be used to create matrix from the
 * network component objects
 * @param network: network that will generate matrix
 */
FullMatrixMap(boost::shared_ptr<_network> network)
  : p_me (GA_Nodeid()), p_nNodes(GA_Nnodes()), p_network(network)
{
  int                     iSize    = 0;
  int                     jSize    = 0;

  p_nBuses = p_network->numBuses();
  p_nBranches = p_network->numBranches();

  p_activeBuses         = getActiveBuses();

  setupGlobalArrays(p_activeBuses);  // allocate globalIndex arrays

  setupIndexingArrays();

  setupOffsetArrays();

  contributions();
}

~FullMatrixMap()
{
  GA_Destroy(gaOffsetI);
  GA_Destroy(gaOffsetJ);
}

/**
 * Return the number of active buses on this process
 * @return: number of active buses
 */
int getActiveBuses(void)
{
  int nActiveBuses = 0;
  for (int i = 0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      nActiveBuses++;
    }
  }
  return nActiveBuses;
}

/**
 * Allocate the gaMatBlksI and gaMatBlksJ global arrays
 * @param nActiveBuses: the number of active buses on this process
 */
void setupGlobalArrays(int nActiveBuses)
{
  int one = 1;

  p_totalBuses = nActiveBuses;

  GA_Igop(&p_totalBuses,one,"+");

  // the gaMatBlksI and gaMatBlksJ are the
  createIndexGA(&gaMatBlksI, p_totalBuses);
  createIndexGA(&gaMatBlksJ, p_totalBuses);
}

/**
 * Create a global array of integers
 * @param size: size of global array
 */
void createIndexGA(int * handle, int size)
{
  int one = 1;
  *handle = GA_Create_handle();
  GA_Set_data(*handle, one, &size, C_INT);
  if (!GA_Allocate(*handle)) {
    // TODO: some kind of error
  }
  GA_Zero(*handle);
}

/**
 * Set up global arrays that contain all matrix block sizes along the I and J
 * axes. These will then be used to create offset arrays.
 */
void setupIndexingArrays()
{
  int                    * iSizeArray     = NULL;
  int                    * jSizeArray     = NULL;
  int                   ** iIndexArray    = NULL;
  int                   ** jIndexArray    = NULL;
  int                      count          = 0;

  // set up bus indexing
  allocateIndexArray(p_nBuses, &iSizeArray, &jSizeArray, &iIndexArray, NULL, 1);

  loadBusArrays(iSizeArray, jSizeArray, iIndexArray, &count);
  scatterIndexingArrays(iSizeArray, jSizeArray, iIndexArray, NULL, count, 1);
  deleteIndexArrays(p_nBuses, iSizeArray, jSizeArray, iIndexArray, NULL, 1);
  GA_Sync();

  // set up branch indexing
  count               = 0;
  allocateIndexArray(p_nBranches, &iSizeArray, &jSizeArray, &iIndexArray,
      &jIndexArray, 2);
  loadForwardBranchArrays(iSizeArray, jSizeArray, iIndexArray, jIndexArray, &count);
  scatterIndexingArrays(iSizeArray, jSizeArray, iIndexArray, jIndexArray, count, 2);

  count               = 0;
  loadReverseBranchArrays(iSizeArray, jSizeArray, iIndexArray, jIndexArray, &count);
  scatterIndexingArrays(iSizeArray, jSizeArray, iIndexArray, jIndexArray, count, 2);

  deleteIndexArrays(p_nBranches, iSizeArray, jSizeArray, iIndexArray,
      jIndexArray, 2);
  GA_Sync();
}

/**
 * Allocate arrays that hold sizes and approximate indices of matrix elements
 * @param n: number of elements in array
 * @param iSizeArray: array containing size of matrix block along i axis
 * @param jSizeArray: array containing size of matrix block along j axis
 * @param iIndexArray: array containing i index of matrix block
 * @param jIndexArray: array containing j index of matrix block
 * @param nflag: number of indices being used (1 or 2)
 */
void allocateIndexArray(int n, int ** iSizeArray, int ** jSizeArray,
        int *** iIndexArray, int *** jIndexArray, int nflag)
{
  *iSizeArray         = new int[n];
  *jSizeArray         = new int[n];
  *iIndexArray        = new int*[n];

  for(int i = 0; i < n; i++) {
    (*iIndexArray)[i]  = new int;
  }

  if (nflag == 2) {
    *jIndexArray        = new int*[n];
    for(int i = 0; i < n; i++) {
      (*jIndexArray)[i]  = new int;
    }
  }
}

/**
 * Load arrays containing matrix block sizes and indices along diagonal of matrix.
 * These come from buses
 * @param iSizeArray: array containing size of matrix block along i axis
 * @param jSizeArray: array containing size of matrix block along j axis
 * @param iIndexArray: array containing i index of matrix block
 * @param count: total number of non-zero blocks
 */
void loadBusArrays(int * iSizeArray, int * jSizeArray,
        int ** iIndexArray, int *count)
{
  int                      index          = 0;
  int                      iSize          = 0;
  int                      jSize          = 0;
  bool                     status         = true;

  *count = 0;
  for (int i = 0; i < p_nBuses; i++) {
    status = p_network->getBus(i)->matrixDiagSize(&iSize, &jSize);
    if (status) {
      p_network->getBus(i)->getMatVecIndex(&index);
      iSizeArray[*count]     = iSize;
      jSizeArray[*count]     = jSize;
      *(iIndexArray[*count])  = index;
      (*count)++;
    }
  }
}

/**
 * Load arrays containing matrix block sizes and indices for off-diagonal
 * matrix blocks. These come from branches.
 * @param iSizeArray: array containing size of matrix block along i axis
 * @param jSizeArray: array containing size of matrix block along j axis
 * @param iIndexArray: array containing i index of matrix block
 * @param jIndexArray: array containing j index of matrix block
 * @param count: total number of non-zero blocks
 */
void loadForwardBranchArrays(int * iSizeArray, int * jSizeArray,
        int ** iIndexArray, int ** jIndexArray, int * count)
{
  int                      iIndex         = 0;
  int                      jIndex         = 0;
  int                      iSize          = 0;
  int                      jSize          = 0;
  bool                     status         = true;

  *count = 0;
  for (int i = 0; i < p_nBranches; i++) {
    status = p_network->getBranch(i)->matrixForwardSize(&iSize, &jSize);
    if (status) {
      p_network->getBranch(i)->getMatVecIndices(&iIndex, &jIndex);
      iSizeArray[*count]      = iSize;
      jSizeArray[*count]      = jSize;
      *(iIndexArray[*count])  = iIndex;
      *(jIndexArray[*count])  = jIndex;
      (*count)++;
    }
  }
}

void loadReverseBranchArrays(int * iSizeArray, int * jSizeArray,
        int ** iIndexArray, int ** jIndexArray, int * count)
{
  int                      iIndex         = 0;
  int                      jIndex         = 0;
  int                      iSize          = 0;
  int                      jSize          = 0;
  bool                     status         = true;

  *count = 0;
  for (int i = 0; i < p_nBranches; i++) {
    status = p_network->getBranch(i)->matrixReverseSize(&iSize, &jSize);
    if (status) {
      p_network->getBranch(i)->getMatVecIndices(&iIndex, &jIndex);
      iSizeArray[*count]      = iSize;
      jSizeArray[*count]      = jSize;
      *(iIndexArray[*count])  = jIndex;
      *(jIndexArray[*count])  = iIndex;
      (*count)++;
    }
  }
}

/**
 *  Clean up index arrays
 *  @param n: array size
 *  @param iSizeArray: array containing size of matrix block along i axis
 *  @param jSizeArray: array containing size of matrix block along j axis
 *  @param iIndexArray: array containing i index of matrix block
 *  @param jIndexArray: array containing j index of matrix block
 * @param nflag: number of indices being used (1 or 2)
 */
void deleteIndexArrays(int n, int * iSizeArray, int * jSizeArray,
        int ** iIndexArray, int ** jIndexArray, int nflag)
{
  for(int i = 0; i < n; i++) {
    delete iIndexArray[i];
  }
  delete [] iIndexArray;

  if (nflag == 2) {
    for(int i = 0; i < n; i++) {
      delete jIndexArray[i];
    }
    delete [] jIndexArray;
  }

  delete [] iSizeArray;
  delete [] jSizeArray;
}

/**
 * Scatter elements into global arrays
 * @param iSizeArray: array containing size of matrix block along i axis
 * @param jSizeArray: array containing size of matrix block along j axis
 * @param iIndexArray: array containing i index of matrix block
 * @param jIndexArray: array containing j index of matrix block
 * @param count: number of elements to be scattered
 * @param nflag: number of indices being used (1 or 2)
 */
void scatterIndexingArrays(int * iSizeArray, int * jSizeArray,
                                  int ** iIndexArray, int ** jIndexArray,
                                  int count, int nflag)
{
  if (count > 0) NGA_Scatter(gaMatBlksI, iSizeArray, iIndexArray, count);
  if (count > 0 && nflag == 2)
     NGA_Scatter(gaMatBlksJ, jSizeArray, jIndexArray, count);
}

/**
 * Set up the offset arrays that will be used to find the exact location of
 * each matrix block in the matrix produced by the mapper
 */
void setupOffsetArrays()
{
  int *itmp = new int[p_nNodes];
  int *jtmp = new int[p_nNodes];

  int *rptr;
  int i, idx;
  int one = 1;

  // We need to decompose this matrix by rows so that each processor has a
  // contiguous set of rows. The MatVec indices are set up so that the active
  // buses on each processor are indexed consecutively. This index can be used
  // to set up this partition.
  p_minRowIndex = p_totalBuses;
  p_maxRowIndex = 0;
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      p_network->getBus(i)->getMatVecIndex(&idx);
      if (idx > p_maxRowIndex) p_maxRowIndex = idx;
      if (idx < p_minRowIndex) p_minRowIndex = idx;
    }
  }
//  printf("p[%d] (FullMatrixMap) minRow: %d maxRow: %d\n",p_me,p_minRowIndex,p_maxRowIndex);
  // Create array to hold information about desired rows
  int nRows = p_maxRowIndex-p_minRowIndex+1;
  int *iSizes = new int[nRows]; 
  int *jSizes = new int[nRows]; 
  GA_Sync();
  NGA_Get(gaMatBlksI,&p_minRowIndex,&p_maxRowIndex,iSizes,&one);
  NGA_Get(gaMatBlksJ,&p_minRowIndex,&p_maxRowIndex,jSizes,&one);

  // Calculate total number of elements associated with row block and column
  // block associated with this processor
  int iSize = 0;
  int jSize = 0;
  p_maxIBlock = 0;
  p_maxJBlock = 0;
  for (i=0; i<nRows; i++) {
    if (p_maxIBlock < iSizes[i]) p_maxIBlock = iSizes[i];
    if (p_maxJBlock < jSizes[i]) p_maxJBlock = jSizes[i];
    if (iSizes[i] > 0) iSize += iSizes[i];
    if (jSizes[i] > 0) jSize += jSizes[i];
    if (iSizes[i] == 0 || jSizes[i] == 0) {
//      printf("p[%d] Sizes[%d] I: %d J: %d\n",p_me,i,iSizes[i],jSizes[i]);
    }
  }
  p_rowBlockSize = iSize;
  GA_Igop(&p_maxIBlock,one,"max");
  GA_Igop(&p_maxJBlock,one,"max");

  for (i = 0; i<p_nNodes; i++) {
    itmp[i] = 0;
    jtmp[i] = 0;
  }
  itmp[p_me] = iSize;
  jtmp[p_me] = jSize;
//  printf("p[%d] (FullMatrixMap) iSize: %d jSize: %d\n",p_me,iSize,jSize);

  GA_Igop(itmp, p_nNodes, "+");
  GA_Igop(jtmp, p_nNodes, "+");

  int offsetArrayISize = 0;
  int offsetArrayJSize = 0;
  for (i = 0; i < p_me; i++) {
    offsetArrayISize += itmp[i];
    offsetArrayJSize += jtmp[i];
  }

  // Calculate matrix dimension
  p_iDim = 0;
  p_jDim = 0;
  for (i=0; i<p_nNodes; i++) {
    p_iDim += itmp[i];
    p_jDim += jtmp[i];
  }
//  printf("p[%d] (FullMatrixMap) iDim: %d jDim: %d\n",p_me,p_iDim,p_jDim);

  // Create map array so that offset arrays can be created with a specified
  // distribution
  int *offset = new int[p_nNodes];
  for (i=0; i<p_nNodes; i++) {
    offset[i] = 0;
  }
  offset[p_me] = p_activeBuses;
//  printf("p[%d] (FullMatrixMap) activeBuses: %d\n",p_me,p_activeBuses);
  GA_Igop(offset,p_nNodes,"+");

  int *mapc = new int[p_nNodes];
  mapc[0]=0;
  for (i=1; i<p_nNodes; i++) {
    mapc[i] = mapc[i-1] + offset[i-1];
  }

  delete [] offset;
  delete [] itmp;
  delete [] jtmp;

  // Create arrays that will hold final offsets
  gaOffsetI = GA_Create_handle();
  GA_Set_data(gaOffsetI, one, &p_totalBuses, C_INT);
  GA_Set_irreg_distr(gaOffsetI, mapc, &p_nNodes);
  if (!GA_Allocate(gaOffsetI)) {
    // TODO: some kind of error
  }
  GA_Zero(gaOffsetI);

  gaOffsetJ = GA_Create_handle();
  GA_Set_data(gaOffsetJ, one, &p_totalBuses, C_INT);
  GA_Set_irreg_distr(gaOffsetJ, mapc, &p_nNodes);
  if (!GA_Allocate(gaOffsetJ)) {
    // TODO: some kind of error
  }
  GA_Zero(gaOffsetJ);

  // Evaluate offsets for this processor
  int *iOffsets = new int[nRows]; 
  int *jOffsets = new int[nRows]; 
  iOffsets[0] = offsetArrayISize;
  jOffsets[0] = offsetArrayJSize;
  for (i=1; i<nRows; i++) {
    iOffsets[i] = iOffsets[i-1] + iSizes[i-1];
    jOffsets[i] = jOffsets[i-1] + jSizes[i-1];
  }

  // Put offsets into global arrays
  if (nRows > 0) {
    NGA_Put(gaOffsetI,&p_minRowIndex,&p_maxRowIndex,iOffsets,&one);
    NGA_Put(gaOffsetJ,&p_minRowIndex,&p_maxRowIndex,jOffsets,&one);
  }
  GA_Sync();

  // Clean up arrays that are no longer needed
  GA_Destroy(gaMatBlksI);
  GA_Destroy(gaMatBlksJ);

  delete [] mapc;
  delete [] iSizes;
  delete [] jSizes;
  delete [] iOffsets;
  delete [] jOffsets;
}

// map to matrix
//void mapToMatrix(gridpack::math::matrix    & matrix)
boost::shared_ptr<gridpack::math::Matrix> mapToMatrix(void)
{
  gridpack::parallel::Communicator comm = p_network->communicator();
  boost::shared_ptr<gridpack::math::Matrix>
             Ret(new gridpack::math::Matrix(comm,
             p_rowBlockSize, p_jDim, gridpack::math::Matrix::Sparse));
  loadBusData(Ret);
  loadBranchData(Ret);
  GA_Sync();
  Ret->ready();
  return Ret;
}

/**
 * Add diagonal block contributions from buses to matrix
 * @param matrix: matrix to which contributions are added
 */
void loadBusData(boost::shared_ptr<gridpack::math::Matrix> matrix)
{
  int i,idx,jdx,isize,jsize,icnt;
  int **indices = new int*[p_busContribution];
  icnt = 0;
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      if (p_network->getBus(i)->matrixDiagSize(&isize,&jsize)) {
        indices[icnt] = new int;
        p_network->getBus(i)->getMatVecIndex(&idx);
        *(indices[icnt]) = idx;
        icnt++;
      }
    }
  }
  if (icnt != p_busContribution) {
    // TODO: some kind of error
  }

  // Gather matrix offsets
  int *offsets = new int[p_busContribution];
  if (p_busContribution > 0) {
    NGA_Gather(gaOffsetI,offsets,indices,p_busContribution);
  }

  // Add matrix elements
  ComplexType *values = new ComplexType[p_maxIBlock*p_maxJBlock];
  int j,k;
  int jcnt = 0;
  int acnt = 0;
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      if (p_network->getBus(i)->matrixDiagSize(&isize,&jsize)) {
        p_network->getBus(i)->matrixDiagValues(values);
        icnt = 0;
        for (k=0; k<jsize; k++) {
          jdx = offsets[jcnt] + k;
          for (j=0; j<isize; j++) {
            idx = offsets[jcnt] + j;
            matrix->add_element(idx, jdx, values[icnt]);
            icnt++;
          }
        }
        delete indices[jcnt];
        jcnt++;
      }
      acnt++;
    }
  }

  // Clean up arrays
  delete [] indices;
  delete [] offsets;
  delete [] values;
}

/**
 * Add off-diagonal block contributions from branches to matrix
 * @param matrix: matrix to which contributions are added
 */
void loadBranchData(boost::shared_ptr<gridpack::math::Matrix> matrix)
{
  int i,idx,jdx,isize,jsize,icnt;
  int **i_indices = new int*[p_branchContribution];
  int **j_indices = new int*[p_branchContribution];
  icnt = 0;
  for (i=0; i<p_nBranches; i++) {
    if (p_network->getBranch(i)->matrixForwardSize(&isize,&jsize)) {
      p_network->getBranch(i)->getMatVecIndices(&idx, &jdx);
      if (idx >= p_minRowIndex && idx <= p_maxRowIndex) {
        i_indices[icnt] = new int;
        j_indices[icnt] = new int;
        *(i_indices[icnt]) = idx;
        *(j_indices[icnt]) = jdx;
        icnt++;
      } else {
        // TODO: some kind of error
      }
    }
    if (p_network->getBranch(i)->matrixReverseSize(&isize,&jsize)) {
      p_network->getBranch(i)->getMatVecIndices(&idx, &jdx);
      if (jdx >= p_minRowIndex && jdx <= p_maxRowIndex) {
        i_indices[icnt] = new int;
        j_indices[icnt] = new int;
        *(i_indices[icnt]) = jdx;
        *(j_indices[icnt]) = idx;
        icnt++;
      }
    }
  }
  if (icnt != p_branchContribution) {
    printf("p[%d] Mismatch in loadBranchData icnt: %d branchContribution: %d\n",
        p_me,icnt,p_branchContribution);
    // TODO: some kind of error
  }

  // Gather matrix offsets
  int *i_offsets = new int[p_branchContribution];
  int *j_offsets = new int[p_branchContribution];
  if (p_branchContribution > 0) {
    NGA_Gather(gaOffsetI,i_offsets,i_indices,p_branchContribution);
    NGA_Gather(gaOffsetJ,j_offsets,j_indices,p_branchContribution);
  }

  // Add matrix elements
  ComplexType *values = new ComplexType[p_maxIBlock*p_maxJBlock];
  int j,k;
  int jcnt = 0;
  for (i=0; i<p_nBranches; i++) {
    if (p_network->getBranch(i)->matrixForwardSize(&isize,&jsize)) {
      p_network->getBranch(i)->getMatVecIndices(&idx, &jdx);
      if (idx >= p_minRowIndex && idx <= p_maxRowIndex) {
        p_network->getBranch(i)->matrixForwardValues(values);
        icnt = 0;
        for (k=0; k<jsize; k++) {
          jdx = j_offsets[jcnt] + k;
          for (j=0; j<isize; j++) {
            idx = i_offsets[jcnt] + j;
//            printf("p[%d] (1) idx: %d jdx: %d\n",p_me,idx,jdx);
//            if (idx >= p_iDim || idx < 0 || jdx >= p_jDim || jdx < 0 ) {
//              printf("p[%d] (1) idx: %d jdx: %d\n",p_me,idx,jdx);
//              matrix->add_element(idx, jdx, values[icnt]);
//              printf("p[%d] (1) finished idx: %d jdx: %d\n",p_me,idx,jdx);
//            } else {
              matrix->add_element(idx, jdx, values[icnt]);
//            }
            icnt++;
          }
        }
        delete i_indices[jcnt];
        delete j_indices[jcnt];
        jcnt++;
      }
    }
    if (p_network->getBranch(i)->matrixReverseSize(&isize,&jsize)) {
      p_network->getBranch(i)->getMatVecIndices(&idx, &jdx);
      if (jdx >= p_minRowIndex && jdx <= p_maxRowIndex) {
        p_network->getBranch(i)->matrixReverseValues(values);
        icnt = 0;
        for (k=0; k<jsize; k++) {
          jdx = j_offsets[jcnt] + k;
          for (j=0; j<isize; j++) {
            idx = i_offsets[jcnt] + j;
//            printf("p[%d] (2) idx: %d jdx: %d\n",p_me,idx,jdx);
//            if (idx >= p_iDim || idx < 0 || jdx >= p_jDim || jdx < 0 ) {
//              printf("p[%d] (2) idx: %d jdx: %d\n",p_me,idx,jdx);
//              matrix->add_element(idx, jdx, values[icnt]);
//              printf("p[%d] (2) finished idx: %d jdx: %d\n",p_me,idx,jdx);
//            } else {
              matrix->add_element(idx, jdx, values[icnt]);
//            }
            icnt++;
          }
        }
        delete i_indices[jcnt];
        delete j_indices[jcnt];
        jcnt++;
      }
    }
  }

  // Clean up arrays
  delete [] i_indices;
  delete [] j_indices;
  delete [] i_offsets;
  delete [] j_offsets;
  delete [] values;
}

/**
 * Calculate how many buses and branches contribute to matrix
 */
void contributions(void)
{
  int i;
  // Get number of contributions from buses
  int isize, jsize;
  p_busContribution = 0;
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      if (p_network->getBus(i)->matrixDiagSize(&isize, &jsize)) p_busContribution++;
    }
  }

  // Get number of contributions from branches
  int idx, jdx;
  p_branchContribution = 0;
  for (i=0; i<p_nBranches; i++) {
    if (p_network->getBranch(i)->matrixForwardSize(&isize, &jsize)) {
      p_network->getBranch(i)->getMatVecIndices(&idx, &jdx);
      if (idx >= p_minRowIndex && idx <= p_maxRowIndex) {
        p_branchContribution++;
      } else {
        // TODO: some kind of error
      }
    }
    if (p_network->getBranch(i)->matrixReverseSize(&isize, &jsize)) {
      p_network->getBranch(i)->getMatVecIndices(&idx, &jdx);
      if (jdx >= p_minRowIndex && jdx <= p_maxRowIndex) {
        p_branchContribution++;
      }
    }
  }
}

private:
    // GA information
int                         p_me;
int                         p_nNodes;

    // network information
boost::shared_ptr<_network> p_network;
int                         p_nBuses;
int                         p_nBranches;
int                         p_totalBuses;
int                         p_activeBuses;

    // matrix information
int                         p_iDim;
int                         p_jDim;
int                         p_minRowIndex;
int                         p_maxRowIndex;
int                         p_rowBlockSize;
int                         p_minColIndex;
int                         p_maxColIndex;
int                         p_busContribution;
int                         p_branchContribution;
int                         p_maxIBlock;
int                         p_maxJBlock;

    // global matrix block size array
int                         gaMatBlksI; // g_idx
int                         gaMatBlksJ; // g_jdx
int                         gaOffsetI; // g_ioff
int                         gaOffsetJ; // g_joff
};

} /* namespace mapper */
} /* namespace gridpack */

#endif //FULLMATRIXMAP_HPP_
