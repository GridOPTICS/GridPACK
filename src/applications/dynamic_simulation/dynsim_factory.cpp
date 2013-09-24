// -------------------------------------------------------------
/**
 * @file   dynsim_factory.cpp
 * @author Shuangshuang Jin 
 * @date   September 19, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/applications/dynamic_simulation/dynsim_components.hpp"
#include "gridpack/applications/dynamic_simulation/dynsim_factory.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"

namespace gridpack {
namespace dynsim {

// Powerflow factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
DynSimFactory::DynSimFactory(DynSimFactory::NetworkPtr network)
  : gridpack::factory::BaseFactory<DynSimNetwork>(network)
{
  p_network = network;
}

/**
 * Basic destructor
 */
gridpack::dynsim::DynSimFactory::~DynSimFactory()
{
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::dynsim::DynSimFactory::setYBus(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // Invoke setYBus method on all bus objects
  for (i=0; i<numBus; i++) {
    (dynamic_cast<gridpack::dynsim::DynSimBus*>(p_network->getBus(i).get()))->setYBus();
  }

  // Invoke setYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
    (dynamic_cast<gridpack::dynsim::DynSimBranch*>(p_network->getBranch(i).get()))->setYBus();
  }
}

} // namespace dynsim
} // namespace gridpack
