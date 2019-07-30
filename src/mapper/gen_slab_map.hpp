/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/* *****************************************************************************
 * gen_slab_map.hpp
 * gridpack
 * Bruce Palmer
 * July 23, 2014
 *
 * Documentation status
 * Revision status
 * File contents
 **************************************************************************** */

#ifndef GENSLABMAP_HPP_
#define GENSLABMAP_HPP_

#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
#include "gridpack/parallel/parallel.hpp"
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/component/base_component.hpp>
#include <gridpack/component/matvec_ifc.hpp>
#include <gridpack/mapper/base_gen_slab_map.hpp>
#include <gridpack/network/base_network.hpp>
#include <gridpack/utilities/exception.hpp>
#include <gridpack/math/vector.hpp>

//#define DBG_CHECK

namespace gridpack {
namespace mapper {

template <class _network>
class GenSlabMap : BaseGenSlabMap<_network> {
  public: 
/**
 * Initialize mapper for the given network and the current mode. Create global
 * arrays that contain offsets that will be used to create vector from the
 * network component objects
 * @param network network that will generate vector
 */
GenSlabMap(boost::shared_ptr<_network> network)
  : BaseGenSlabMap<_network>::BaseGenSlabMap(network)
{
  this->setMode(gridpack::component::BaseGenMatVecInterface::STANDARD);
}

~GenSlabMap()
{
}

/**
 * Generate a dense matrix from current component state on network
 * @return return a pointer to new matrix
 */
boost::shared_ptr<gridpack::math::Matrix> mapToMatrix(void)
{
  return this->baseMapToMatrix();
}

/**
 * Reset existing matrix from current component state on network
 * @param matrix existing matrix (should be generated from same mapper)
 */
void mapToMatrix(gridpack::math::Matrix &matrix)
{
  this->baseMapToMatrix(matrix);
}

/**
 * Reset existing matrix from current component state on network
 * @param matrix existing vector (should be generated from same mapper)
 */
void mapToMatrix(boost::shared_ptr<gridpack::math::Matrix> &matrix)
{
  this->baseMapToMatrix(matrix);
}

/**
 * Push data from matrix onto buses and branches. Matrix must be
 * created with the mapToMatrix method using the same GenSlabMap
 * @param matrix matrix containing data to be pushed to network
 */
void mapToNetwork(const gridpack::math::Matrix &matrix)
{
  this->baseMapToNetwork(matrix);
}

/**
 * Push data from matrix onto buses and branches. Matrix must be
 * created with the mapToMatrix method using the same GenSlabMap
 * @param matrix matrix containing data to be pushed to network
 */
void mapToNetwork(boost::shared_ptr<gridpack::math::Matrix> &matrix)
{
  this->baseMapToNetwork(matrix);
}

private:

};

} /* namespace mapper */
} /* namespace gridpack */

#endif //GENSLABMAP_HPP_
