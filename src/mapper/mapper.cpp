/* *****************************************************************************
 * mapper.cpp
 * gridpack
 * kglass
 * Jul 22, 2013
 *
 * Documentation status
 * Revision status
 * File contents
 **************************************************************************** */

#include "mapper.hpp"

namespace gridpack {
namespace mapper {

/**
 * Initialize mapper for the given network and the current mode. Create global
 * arrays that contain offsets that will be used to create matrix from the
 * network component objects
 * @param network: network that will generate matrix
 */
Mapper::Mapper(boost::shared_ptr<MapNetwork> network)
  : p_me (GA_Nodeid()), p_nNodes(GA_Nnodes()), p_network(network)
{
  int                     iSize    = 0;
  int                     jSize    = 0;

  int p_nBuses = p_network->numBuses();
  int p_nBranches = p_network->numBranches();

  int nActiveBuses         = getActiveBuses();

  setupGlobalArrays(nActiveBuses);  // allocate globalIndex arrays

  setupIndexingArrays();

  setupOffsetArrays();
}

Mapper::~Mapper()
{
  GA_Destroy(gaOffsetI);
  GA_Destroy(gaOffsetJ);
}

/**
 * Return the number of active buses on this process
 * @return: number of active buses
 */
int Mapper::getActiveBuses(void)
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
void Mapper::setupGlobalArrays(int nActiveBuses)
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
void Mapper::createIndexGA(int * handle, int size)
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
void Mapper::setupIndexingArrays()
{
  int                    * iSizeArray     = NULL;
  int                    * jSizeArray     = NULL;
  int                   ** iIndexArray    = NULL;
  int                   ** jIndexArray    = NULL;
  int                      count          = 0;

  // set up bus indexing
  allocateIndexArray(p_nBuses, &iSizeArray, &jSizeArray, &iIndexArray, NULL);

  loadBusArrays(iSizeArray, jSizeArray, iIndexArray, &count);
  scatterIndexingArrays(iSizeArray, jSizeArray, iIndexArray, NULL, count);
  deleteIndexArrays(p_nBuses, iSizeArray, jSizeArray, iIndexArray, NULL);
  GA_Sync();

  // set up branch indexing
  count               = 0;
  allocateIndexArray(p_nBranches, &iSizeArray, &jSizeArray, &iIndexArray, &jIndexArray);
  loadForwardBranchArrays(iSizeArray, jSizeArray, iIndexArray, jIndexArray, &count);
  scatterIndexingArrays(iSizeArray, jSizeArray, iIndexArray, jIndexArray, count);

  count               = 0;
  loadReverseBranchArrays(iSizeArray, jSizeArray, iIndexArray, jIndexArray, &count);
  scatterIndexingArrays(iSizeArray, jSizeArray, iIndexArray, jIndexArray, count);

  deleteIndexArrays(p_nBranches, iSizeArray, jSizeArray, iIndexArray, jIndexArray);
  GA_Sync();
}

/**
 * Allocate arrays that hold sizes and approximate indices of matrix elements
 * @param n: number of elements in array
 * @param iSizeArray: array containing size of matrix block along i axis
 * @param jSizeArray: array containing size of matrix block along j axis
 * @param iIndexArray: array containing i index of matrix block
 * @param jIndexArray: array containing j index of matrix block
 */
void Mapper::allocateIndexArray(int n, int ** iSizeArray, int ** jSizeArray,
        int *** iIndexArray, int *** jIndexArray)
{
  *iSizeArray         = new int[n];
  *jSizeArray         = new int[n];
  *iIndexArray        = new int*[n];

  for(int i = 0; i < n; i++) {
    (*iIndexArray)[i]  = new int;
  }

  if (jSizeArray != NULL) {
    *jSizeArray         = new int[n];
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
void Mapper::loadBusArrays(int * iSizeArray, int * jSizeArray,
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
void Mapper::loadForwardBranchArrays(int * iSizeArray, int * jSizeArray,
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

void Mapper::loadReverseBranchArrays(int * iSizeArray, int * jSizeArray,
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
 */
void Mapper::deleteIndexArrays(int n, int * iSizeArray, int * jSizeArray,
        int ** iIndexArray, int ** jIndexArray)
{
  for(int i = 0; i < n; i++) {
    delete iIndexArray[i];
  }
  delete [] iIndexArray;

  if (jIndexArray != NULL) {
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
 */
void Mapper::scatterIndexingArrays(int * iSizeArray, int * jSizeArray,
                                  int ** iIndexArray, int ** jIndexArray, int count)
{
  if (count > 0) NGA_Scatter(gaMatBlksI, iSizeArray, iIndexArray, count);
  if (count > 0) NGA_Scatter(gaMatBlksJ, jSizeArray, jIndexArray, count);
}

/**
 * Set up the offset arrays that will be used to find the exact location of
 * each matrix block in the matrix produced by the mapper
 */
void Mapper::setupOffsetArrays()
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
    p_network->getBus(i)->getMatVecIndex(&idx);
    if (idx > p_maxRowIndex) p_maxRowIndex = idx;
    if (idx < p_minRowIndex) p_minRowIndex = idx;
  }
  // Create array to hold information about desired rows
  int nRows = p_maxRowIndex-p_minRowIndex+1;
  int *iSizes = new int[nRows]; 
  int *jSizes = new int[nRows]; 
  NGA_Get(gaMatBlksI,&p_minRowIndex,&p_maxRowIndex,iSizes,&one);
  NGA_Get(gaMatBlksJ,&p_minRowIndex,&p_maxRowIndex,jSizes,&one);

  // Calculate total number of elements associated with row block and column
  // block associated with this processor
  int iSize = 0;
  int jSize = 0;
  for (i=0; i<nRows; i++) {
    if (iSizes[i] > 0) iSize += iSizes[i];
    if (jSizes[i] > 0) jSize += jSizes[i];
  }

  for (i = 0; i<p_nNodes; i++) {
    itmp[i] = 0;
    jtmp[i] = 0;
  }
  itmp[p_me] = iSize;
  jtmp[p_me] = jSize;

  GA_Igop(itmp, p_nNodes, "+");
  GA_Igop(jtmp, p_nNodes, "+");

  int offsetArrayISize = 0;
  int offsetArrayJSize = 0;
  for (i = 0; i < p_me; i++) {
    offsetArrayISize += itmp[i];
    offsetArrayJSize += jtmp[i];
  }

  // Calculate matrix dimension
  int p_iDim = 0;
  int p_jDim = 0;
  for (i=0; i<p_nNodes; i++) {
    p_iDim += itmp[i];
    p_jDim += jtmp[i];
  }

  // Create map array so that offset arrays can be created with a specified
  // distribution
  int *mapc = new int[p_nNodes];
  mapc[0]=0;
  for (i=1; i<p_nNodes; i++) {
    mapc[i] = mapc[i-1] + itmp[i-1];
  }

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
void mapToMatrix(void)
{
//    loadBusData(matrix);
//    loadBranchData(matrix);

}

void Mapper::loadBusData(gridpack::math::Matrix    & matrix)
{
#if 0
    // loop on p_network to get submatrix
    for (int i = 0; i < p_totalBuses; i++) {
        // using the bus block set, get matrix offset for the block
        DataCollection  * data = p_network.getBusData();
//        int   rowOffset = ?; // get bus row offset from GA
//        int   colOffset = ?; // get bus col offset from GA
        // add values to matrix
        matrix.add(rowOffset, colOffset, data);
    }
#endif
}

void loadBranchData(gridpack::math::Matrix    & matrix)
{
#if 0
    for (int i = 0; i < p_totalBuses; i++) {
        // using the bus block set, get matrix offset for the block
        DataCollection  * data = p_network.getBranchData();
//        int   rowOffset = ?; // get branch row offset from GA
//        int   colOffset = ?; // get branch col offset from GA
        // add values to matrix
        matrix.add(rowOffset, colOffset, data);
    }
#endif
}

} /* namespace mapper */
} /* namespace gridpack */




