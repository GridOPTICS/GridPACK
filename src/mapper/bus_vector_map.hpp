/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/* *****************************************************************************
 * bus_vector_map.hpp
 * gridpack
 * kglass, bjpalmer
 * Jul 22, 2013
 *
 * Documentation status
 * Revision status
 * File contents
 **************************************************************************** */

#ifndef BUSVECTORMAP_HPP_
#define BUSVECTORMAP_HPP_

#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
#include "gridpack/parallel/parallel.hpp"
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/component/base_component.hpp>
#include <gridpack/component/matvec_ifc.hpp>
#include <gridpack/mapper/base_vector_map.hpp>
#include <gridpack/network/base_network.hpp>
#include <gridpack/math/vector.hpp>

//#define DBG_CHECK

namespace gridpack {
namespace mapper {

template <class _network>
class BusVectorMap : BaseVectorMap<_network> {
  public: 
/**
 * Initialize mapper for the given network and the current mode. Create global
 * arrays that contain offsets that will be used to create vector from the
 * network component objects
 * @param network network that will generate vector
 */
BusVectorMap(boost::shared_ptr<_network> network)
  : BaseVectorMap<_network>::BaseVectorMap(network)
{
  this->setMode(gridpack::component::BaseMatrixInterface::STANDARD);
}

~BusVectorMap()
{
}

/**
 * Create a vector from the current bus state
 * @return return pointer to new vector
 */
boost::shared_ptr<gridpack::math::Vector> mapToVector(void)
{
  return this->baseMapToVector();
}

/**
 * Create a vector from the current bus state
 * @return return pointer to new vector
 */
boost::shared_ptr<gridpack::math::RealVector> mapToRealVector(void)
{
  return this->baseMapToRealVector();
}

/**
 * Reset a vector from the current bus state (vector should be created with same
 * mapper)
 * @param vector existing vector that should be reset
 */
void mapToVector(gridpack::math::Vector &vector)
{
  this->baseMapToVector(vector);
}

/**
 * Reset a vector from the current bus state (vector should be created with same
 * mapper)
 * @param vector existing vector that should be reset
 */
void mapToRealVector(gridpack::math::RealVector &vector)
{
  this->baseMapToRealVector(vector);
}

/**
 * Reset a vector from the current bus state (vector should be created with same
 * mapper). 
 * @param vector existing vector that should be reset
 */
void mapToVector(boost::shared_ptr<gridpack::math::Vector> &vector)
{
  this->baseMapToVector(*vector);
}

/**
 * Reset a vector from the current bus state (vector should be created with same
 * mapper). 
 * @param vector existing vector that should be reset
 */
void mapToRealVector(boost::shared_ptr<gridpack::math::RealVector> &vector)
{
  this->baseMapToRealVector(*vector);
}

/**
 * Push data from vector onto buses. Vector must be created with the
 * mapToVector method using the same BusVectorMap
 * @param vector vector containing data to be pushed to buses
 */
void mapToBus(const gridpack::math::Vector &vector)
{
  this->baseMapToBus(vector);
}

/**
 * Push data from vector onto buses. Vector must be created with the
 * mapToVector method using the same BusVectorMap
 * @param vector vector containing data to be pushed to buses
 */
void mapToBus(const gridpack::math::RealVector &vector)
{
  this->baseMapToBus(vector);
}

/**
 * Push data from vector onto buses. Vector must be created with the
 * mapToVector method using the same BusVectorMap
 * @param vector vector containing data to be pushed to buses
 */
void mapToBus(boost::shared_ptr<gridpack::math::Vector> &vector)
{
  this->baseMapToBus(*vector);
}

void mapToBus(boost::shared_ptr<gridpack::math::RealVector> &vector)
{
  this->baseMapToBus(*vector);
}

private:

};

} /* namespace mapper */
} /* namespace gridpack */

#endif //BUSVECTORMAP_HPP_
