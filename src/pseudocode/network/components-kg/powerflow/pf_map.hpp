/* *****************************************************************************
 * class static_map
 *
 * Responsibilities:
 * 1) create and populate matrix
 * 2)
 *
 **************************************************************************** */

#ifndef GENERICMAP_HPP_
#define GENERICMAP_HPP_

#include "gridpack/math/matrix.hpp"


namespace gridpack {
namespace map {

class PFMap : public GenericMap
{
public:
	PFMap(network::base_network * network) :
        network_(network){};

    virtual ~PFMap(){};

    /*
     * The problem object is problem dependent. The mapper
     * subclass will query the network using the appropriate
     * Map visitors (MapSize and MapData) to get the necessary
     * size and value information from the network components
     */
    virtual void mapNetworkToProblem(PROBLEM & problem)
    {
        int             * size[2]
        MapSize           mapSize;
        network_->networkMap(mapSize);
        mapSize->getSize(size);
        // allocate space for Y, I and V

        PFData           mapData;
        network_->networkMap(mapData);
        problem->setY(mapData->getY());
        problem->setI(mapData->getI());
        problem->setV(mapData->getV());
    }

    // the controlling method invokes the map method to start the mapping process
    void mapSolutionToNetwork(PROBLEM & problem)
    {
        // this will map matrix data back to the network
        this->mapToNetwork_();
    }

protected:
    /*
     * Get size of matrix from MatrixInterfaces in the network
     * Allocate the matrix
     * Populate the matrix with data from the network
     * Pass populated motrix to solver
     */
    virtual math::matrix * mapToMatrix_(void)   = 0;

    /*
     * Get component's matrix location and size
     * Retrieve matrix data at location and size
     *
     * pass matrix data to components MatrixInterface
     */
    virtual void mapToNetwork_(void)  = 0;

    /*
     * Get the grid components from the network, and iterate through them
     * to determine the size of the matrix
     */
    virtual int getSize_(void)   = 0;

    /*
     * Iterate through the network's components and get:
     *      size of matrix contribution
     *      local location of contribution
     *      shape of data
     * send collected information to Matrix
     */
    virtual void populateMatrix_(math::Matrix & matrix) = 0;

private:
    network::base_network  * network_
    // location within the source matrix
    int               source_x;
    int               source_y;
};

} // namespace math
} // namespace gridpack

#endif /* GENERICMAP_HPP_ */
