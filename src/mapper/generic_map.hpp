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

class GenericMap
{
public:
    GenericMap(network::base_network * network) :
        network_(network){};

    virtual ~GenericMap(){};

    // the controlling method invokes the map method to start the mapping process
    math::matrix * mapToMatrix(void)
    {
        // this will allocate and populate the matrices associated with the
        // application implementation
        math::matrix * matrix = this->mapToMatrix_();
	return matrix;
    }

    // the controlling method invokes the map method to start the mapping process
    void mapToNetwork(void)
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
