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

/*
 * Get the active buses and branches from the network
 *  The active buses and branches contribute values to the YBus matrix
 *
 */
int Mapper::setupMapper()
{
    int                     iSize    = 0;
    int                     jSize    = 0;

    int nActiveBuses         = getActiveBuses(nBuses);
    int nActiveBranches      = getActiveBranches(nBranches);

    setupGlobalArrays(nActiveBuses, nActiveBranches);  // allocate globalIndex arrays

    setupIndexingArrays();

    findArrayDimensions(&iSize, &jSize);

    setupOffsetArrays(iSize, jSize);

    return 0;
}

int Mapper::getActiveBuses(int nBuses)
{
    for (int i = 0; i<nBuses; i++) {
      if (network->getActiveBus(i)) {
          nActiveBuses++;
      }
    }

    return 0;
}

int Mapper::getActiveBranches(int nBranches)
{

    for (int i = 0; i<nBranches; i++) {
      if (network->getActiveBranch(i)) {
        nActiveBranches++;
      }
    }

    return 0;
}

/*
 * Allocate the gobalIndexI and gaMatBlksJ global arrays
 */
int Mapper::setupGlobalArrays(int nActiveBuses, int nActiveBranchess)
{
    int                      one            = 1;

    totalBuses           = nActiveBuses;
    totalBranches        = nActiveBranches;

    GA_Igop(&nActiveBuses,one,"+");
    GA_Igop(&nActiveBranches,one,"+");


    // the gaMatBlksI and gaMatBlksJ are the
    createGAHandle(&gaMatBlksI, totalBuses);
    createGAHandle(&gaMatBlksJ, totalBuses);

    return 0;
}

int Mapper::createGAHandle(int * handle, int size)
{
    *handle = GA_Create_handle();
    GA_Set_data(*handle, one, &size, C_INT);
    if (!GA_Allocate(*handle)) {
      // TODO: some kind of error
    }
    GA_Zero(*handle);

    return 0;
}

int Mapper::setupIndexingArrays()
{
    int                    * iSizeArray     = NULL;
    int                    * jSizeArray     = NULL;
    int                   ** iIndexArray    = NULL;
    int                   ** jIndexArray    = NULL;
    int                      count          = 0;

    // set up bus indexing
    allocateIndexArray(nBuses, &iSizeArray, &jSizeArray, &iIndexArray, NULL);

    loadBusArrays(iSizeArray, jSizeArray, iIndexArray, &count);
    scatterIndexingArrays(iSizeArray, jSizeArray, iIndexArray, NULL, count);
    deleteIndexArrays(nBuses, iSizeArray, jSizeArray, iIndexArray, NULL);
    GA_Sync();

    // set up branch indexing
    count               = 0;
    allocateIndexArray(nBranches, &iSizeArray, &jSizeArray, &iIndexArray, &jIndexArray);
    loadForwardBranchArrays(iSizeArray, jSizeArray, iIndexArray, jIndexArray, &count);
    scatterIndexingArrays(iSizeArray, jSizeArray, iIndexArray, jIndexArray, count);

    count               = 0;
    loadReverseBranchArrays(iSizeArray, jSizeArray, iIndexArray, &count);
    scatterIndexingArrays(iSizeArray, jSizeArray, iIndexArray, count);

    deleteIndexArrays(nBranches, iSizeArray, jSizeArray, iIndexArray, jIndexArray);
    GA_Sync();

    return 0;
}

int Mapper::allocateIndexArray(int n, int ** iSizeArray, int ** jSizeArray,
        int *** iIndexArray, int *** jIndexArray)
{
    *iSizeArray         = new int[n];
    *jSizeArray         = new int[n];
    *iIndexArray        = new int*[n];

    for(int i = 0; i < n; i++) {
        iIndexArray[i]  = new int;
    }

    if (jSizeArray != NULL) {
        *jSizeArray         = new int*[n];

        for(int i = 0; i < n; i++) {
            jIndexArray[i]  = new int;
        }
    }
    return 0;
}

/*
 *
 */
int Mapper::loadBusArrays(int * iSizeArray, int * jSizeArray,
        int *** iIndexArray, int *count)
{
    int                      index          = 0;
    int                      iSize          = 0;
    int                      jSize          = 0;
    bool                     status         = true;

    for (int i = 0; i < nBuses; i++) {
      status = network->getBus(i)->getMatrixSize(&iSize, &jSize);
      if (status) {
          network->getBus(i)->getMatVecIndex(&index);
          iSizeArray[*count]     = iSize;
          jSizeArray[*count]     = jSize;
          *(iIndexArray[*count])  = index;
          (*count)++;
      }
    }

    return 0;
}

int Mapper::loadForwardBranchArrays(int * iSizeArray, int * jSizeArray,
        int *** iIndexArray, int *** jIndexArray, int * count)
{
    int                      iIndex         = 0;
    int                      jIndex         = 0;
    int                      iSize          = 0;
    int                      jSize          = 0;
    bool                     status         = true;

    for (int i = 0; i < nBranches; i++) {
      status = network->getBranch(i)->getMatrixSize(&iSize, &jSize);
      if (status) {
          network->getBranch(i)->getMatrixForwardSize(&iIndex, &jIndex);
          iSizeArray[*count]      = iSize;
          jSizeArray[*count]      = jSize;
          *(iIndexArray[*count])  = iIndex;
          *(jIndexArray[*count])  = jIndex;
          (*count)++;
      }
    }

    return 0;
}

int Mapper::loadReverseBranchArrays(int * iSizeArray, int * jSizeArray,
        int *** iIndexArray, int *** jIndexArray, int * count)
{
    int                      iIndex         = 0;
    int                      jIndex         = 0;
    int                      iSize          = 0;
    int                      jSize          = 0;
    bool                     status         = true;

    for (int i = 0; i < nBranches; i++) {
      status = network->getBranch(i)->getMatrixReverseSize(&iSize, &jSize);
      if (status) {
          network->getBranch(i)->getMatVecIndex(&iIndex, &jIndex);
          iSizeArray[*count]      = iSize;
          jSizeArray[*count]      = jSize;
          *(iIndexArray[*count])  = iIndex;
          *(jIndexArray[*count])  = jIndex;
          (*count)++;
      }
    }

    return 0;
}

int Mapper::deleteIndexArrays(int n, int * iSizeArray, int * jSizeArray,
        int *** iIndexArray, int *** jIndexArray)
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

    return 0;
}

int Mapper::scatterIndexingArrays(int * iSizeArray, int * jSizeArray, int ** iIndexArray, int ** jIndexArray, int count)
{
    if (count > 0) NGA_Scatter(gaMatBlksI, iSizeArray, iIndexArray, count);
    if (count > 0) NGA_Scatter(gaMatBlksJ, jSizeArray, jIndexArray, count);

     return 0;
}

int Mapper::findArrayDimensions(int * iSize, int * jSize)
{
    int                      size           = 0;
    int                      one            = 1;

    *iSize    = computeArraySize(gaMatBlksI);
    GA_Igop(size, one, "+");
    *jSize    = computeArraySize(gaMatBlksJ);
    GA_Igop(size, one, "+");

    return 0;
}

int Mapper::setupOffsetArrays(int iSize, int jSize)
{
    int                    * itmp             = new int[nNodes];
    int                    * jtmp             = new int[nNodes];

    offsetArrayISize = 0;
    offsetArrayJSize = 0;

    for (i = 0; i<nNodes; i++) {
      itmp[i] = 0;
      jtmp[i] = 0;
    }
    itmp[me] = iSize;
    jtmp[me] = jSize;

    GA_Igop(itmp, nNodes, "+");
    GA_Igop(jtmp, nNodes, "+");

    offsetArrayISize = 0;
    offsetArrayJSize = 0;
    for (i = 1; i <= me; i++) {
        offsetArrayISize += itmp[i-1];
        offsetArrayJSize += jtmp[i-1];
    }
    delete [] itmp;
    delete [] jtmp;

    // Evaluate individual offsets for each processor along i axis
    NGA_Distribution(g_idx, me, &lo, &hi);
    NGA_Access(g_idx, &lo, &hi, ptr, &ld);
    int *offsetArrayI_idx = new int[hi-lo+1];
    offsetArrayI_idx[0] = offsetArrayISize;
    for (i=1; i<hi-lo+1; i++) {
      offsetArrayI_idx[i] = offsetArrayI_idx[i-1] + ((int*)ptr)[i-1];
    }
    NGA_Release(g_idx, &lo, &hi);

    // Create global array that holds all offsets
    int g_offsetArrayI = GA_Create_handle();
    GA_set_data(g_offsetArrayI, one, &totalBuses, C_INT);
    if (!GA_Allocate(g_offsetArrayI)) {
      // TODO: some kind of error
    }
    GA_Zero(g_offsetArrayI);
    NGA_Put(g_offsetArrayI, &lo, &hi, offsetArrayI_idx, &one);
    delete [] offsetArrayI_idx;

    // Evaluate individual offsets for each processor along j axis
    NGA_Distribution(g_jdx, me, &lo, &hi);
    NGA_Access(g_jdx, &lo, &hi, ptr, &ld);
    int *offsetArrayI_jdx = new int[hi-lo+1];
    offsetArrayJ_idx[0] = offsetArrayJSize;
    for (i=1; i<hi-lo+1; i++) {
      offsetArrayJ_idx[i] = offsetArrayJ_idx[i-1] + ((int*)ptr)[i-1];
    }
    NGA_Release(g_jdx, &lo, &hi);

    // Create global array that holds all offsets
    int g_offsetArrayJ = GA_Create_handle();
    GA_set_data(g_offsetArrayJ, one, &totalBuses, C_INT);
    if (!GA_Allocate(g_offsetArrayJ)) {
      // TODO: some kind of error
    }
    GA_Zero(g_offsetArrayJ);
    NGA_Put(g_offsetArrayJ, &lo, &hi, offsetArrayJ_idx, &one);
    delete [] offsetArrayJ_idx;

    return 0;
}

int Mapper::computeArraySize(int globalIndex)
{
    int                      sum            = 0;
    int                      lo             = 0;
    int                      hi             = 0;
    int                      ld             = 0;
    void                   * ptr            = NULL;

    NGA_Distribution(globalIndex, me, &lo, &hi);
    NGA_Access(globalIndex, &lo, &hi, ptr, &ld);

    sum = 0;
    for (i=0; i<hi-lo+1; i++) {
        sum += ((int*)ptr)[i];
    }
    NGA_Release(globalIndex, &lo, &hi);

    return  sum;
}

// map to matrix
void mapToMatrix(gridpack::math::matrix    & matrix)
{
    loadBusData(matrix);
    loadBranchData(matrix);

}

void Mapper::loadBusData(gridpack::math::matrix    & matrix)
{
    // loop on network to get submatrix
    for (int i = 0; i < totalBuses; i++) {
        // using the bus block set, get matrix offset for the block
        DataCollection  * data = network.getBusData();
        int   rowOffset = ?; // get bus row offset from GA
        int   colOffset = ?; // get bus col offset from GA
        // add values to matrix
        matrix.add(rowOffset, colOffset, data);
    }
}

svoid loadBranchData(gridpack::math::matrix    & matrix)
{
    for (int i = 0; i < totalBuses; i++) {
        // using the bus block set, get matrix offset for the block
        DataCollection  * data = network.getBranchData();
        int   rowOffset = ?; // get branch row offset from GA
        int   colOffset = ?; // get branch col offset from GA
        // add values to matrix
        matrix.add(rowOffset, colOffset, data);
    }
}



Mapper::Mapper(int gaMe, int gaNodes,
        boost::shared_ptr<gridpack::network::BaseNetwork> & network), me(gaMe),
                nNodes(gaNodes), nBuses(network->numBuses()),
                nBranches(network->numBranches()),
                network(network) {};

Mapper::~Mapper(){}

} /* namespace mapper */
} /* namespace gridpack */




