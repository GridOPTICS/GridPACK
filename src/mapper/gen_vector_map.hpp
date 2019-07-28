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
#include <gridpack/component/matvec_ifc.hpp>
#include <gridpack/mapper/base_gen_vector_map.hpp>
#include <gridpack/network/base_network.hpp>
#include <gridpack/math/vector.hpp>
#include <gridpack/utilities/exception.hpp>

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
{
  p_base_mapper.reset(new BaseGenVectorMap<_network>(network));
  p_base_mapper->setMode(gridpack::component::BaseGenMatVecInterface::STANDARD);
}

~GenVectorMap()
{
}

/**
 * Generate vector from current component state on network
 * @return return a pointer to new vector
 */
boost::shared_ptr<gridpack::math::Vector> mapToVector(void)
{
  return p_base_mapper->baseMapToVector();
}

/**
 * Reset existing vector from current component state on network
 * @param vector existing vector (should be generated from same mapper)
 */
void mapToVector(gridpack::math::Vector &vector)
{
  p_base_mapper->baseMapToVector(vector);
}

/**
 * Reset existing vector from current component state on network
 * @param vector existing vector (should be generated from same mapper)
 */
void mapToVector(boost::shared_ptr<gridpack::math::Vector> &vector)
{
  p_base_mapper->baseMapToVector(vector);
}

/**
 * Push data from vector onto buses and branches. Vector must
 * be created with the mapToVector method using the same
 * GenVectorMap
 * @param vector vector containing data to be pushed to network
 */
void mapToNetwork(const gridpack::math::Vector &vector)
{
  p_base_mapper->baseMapToNetwork(vector);
}

/**
 * Push data from vector onto buses and branches. Vector must
 * be created with the mapToVector method using the same
 * GenVectorMap
 * @param vector vector containing data to be pushed to network
 */
void mapToNetwork(boost::shared_ptr<gridpack::math::Vector> &vector)
{
  p_base_mapper->baseMapToNetwork(vector);
}

private:

boost::shared_ptr<gridpack::mapper::BaseGenVectorMap<_network> > p_base_mapper;
 
};

} /* namespace mapper */
} /* namespace gridpack */

#endif //GENVECTORMAP_HPP_
