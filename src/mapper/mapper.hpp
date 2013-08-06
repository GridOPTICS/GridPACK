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
    Mapper(boost::shared_ptr<MapNetwork> network);
    ~Mapper();
    int resetMapper();
    void mapToMatrix();
protected:
    int getActiveBuses();
    void setupGlobalArrays(int nActiveBuses);
    void createIndexGA(int * handle, int size);
    void setupIndexingArrays();
    void allocateIndexArray(int n, int ** iSizeArray, int ** jSizeArray,
            int *** iIndexArray, int *** jIndexArray);
    void loadBusArrays(int * iSizeArray, int * jSizeArray,
            int ** iIndexArray, int *count);
    void loadForwardBranchArrays(int * iSizeArray, int * jSizeArray,
            int ** iIndexArray, int ** jIndexArray, int * count);
    void loadReverseBranchArrays(int * iSizeArray, int * jSizeArray,
            int ** iIndexArray, int ** jIndexArray, int * count);
    void deleteIndexArrays(int n, int * iSizeArray, int * jSizeArray,
            int ** iIndexArray, int ** jIndexArray);
    void scatterIndexingArrays(int * iSizeArray, int * jSizeArray,
            int ** iIndexArray, int ** jIndexArray, int count);
    void setupOffsetArrays();
    void mapToMatrix(gridpack::math::Matrix    & matrix);
    void loadBusData(gridpack::math::Matrix    & matrix);
    void loadBranchData(gridpack::math::Matrix    & matrix);


private:
    // GA information
    int                           p_me;
    int                           p_nNodes;

    // network information
    boost::shared_ptr<MapNetwork> p_network;
    int                           p_nBuses;
    int                           p_nBranches;

    int                           p_totalBuses;

    // matrix information
    int                           p_iDim;
    int                           p_jDim;
    int                           p_minRowIndex;
    int                           p_maxRowIndex;

    // global matrix block size array
    int                           gaMatBlksI; // g_idx
    int                           gaMatBlksJ; // g_jdx
    int                           gaOffsetI;  // g_ioff
    int                           gaOffsetJ;  // g_joff

};

} /* namespace gridpack */
} /* namespace mapper */
#endif /* MAPPER_HPP_ */
