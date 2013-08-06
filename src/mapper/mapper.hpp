/* *****************************************************************************
 * mapper.hpp
 * gridpack
 * kglass
 * Jul 22, 2013
 *
 * Documentation status
 * Revision status
 * File contents
 **************************************************************************** */

#ifndef MAPPER_HPP_
#define MAPPER_HPP_

#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
#include <gridpack/component/base_component.hpp>
#include <gridpack/network/base_network.hpp>
#include <gridpack/math/matrix.hpp>

namespace gridpack {
namespace mapper {

typedef gridpack::network::BaseNetwork<
        gridpack::component::BaseComponent,
        gridpack::component::BaseComponent> MapNetwork;

class Mapper
{
public:
//    Mapper(boost::shared_ptr<MapNetwork> network);
    Mapper(boost::shared_ptr<MapNetwork> network):
            me (GA_Nodeid()), nNodes(GA_Nnodes()), network(network){};
    ~Mapper(){};
    int resetMapper();
    void mapToMatrix();
protected:
    int setupMapper();
    int getActiveBuses(int nBuses);
    int getActiveBranches(int nBranches);
    int setupGlobalArrays(int nActiveBuses, int nActiveBranchess);
    int createGAHandle(int * handle, int size);
    int setupIndexingArrays();
    int allocateIndexArray(int n, int ** iSizeArray, int ** jSizeArray,
            int *** iIndexArray, int *** jIndexArray);
    int loadBusArrays(int * iSizeArray, int * jSizeArray,
            int *** iIndexArray, int *count);
    int loadForwardBranchArrays(int * iSizeArray, int * jSizeArray,
            int *** iIndexArray, int *** jIndexArray, int * count);
    int loadReverseBranchArrays(int * iSizeArray, int * jSizeArray,
            int *** iIndexArray, int *** jIndexArray, int * count);
    int deleteIndexArrays(int n, int * iSizeArray, int * jSizeArray,
            int *** iIndexArray, int *** jIndexArray);
    int scatterIndexingArrays(int * iSizeArray, int * jSizeArray, int ** iIndexArray, int ** jIndexArray, int count);
    int findArrayDimensions(int * iSize, int * jSize);
    int setupOffsetArrays(int iSize, int jSize);
    int computeArraySize(int globalIndex);
    void mapToMatrix(gridpack::math::Matrix    & matrix);
    void loadBusData(gridpack::math::Matrix    & matrix);
    void loadBranchData(gridpack::math::Matrix    & matrix);


private:
    // GA information
    int                          me;
    int                          nNodes;

    // network information
    boost::shared_ptr<MapNetwork> network;
    int                          nBuses;
    int                          nBranches;

    int                          totalBuses;
    int                          totalBranches;

    // global matrix block size array
    int                          gaMatBlksI; // g_idx
    int                          gaMatBlksJ; // g_jdx
    int                          gaOffsetI;  // g_ioff
    int                          gaOffsetJ;  // g_joff

};

} /* namespace gridpack */
} /* namespace mapper */
#endif /* MAPPER_HPP_ */
