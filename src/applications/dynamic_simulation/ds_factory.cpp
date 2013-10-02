// -------------------------------------------------------------
/**
 * @file   ds_factory.cpp
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
#include "gridpack/applications/dynamic_simulation/ds_components.hpp"
#include "gridpack/applications/dynamic_simulation/ds_factory.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"

namespace gridpack {
namespace dynamic_simulation {

// Powerflow factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
DSFactory::DSFactory(DSFactory::NetworkPtr network)
  : gridpack::factory::BaseFactory<DSNetwork>(network)
{
  p_network = network;
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::DSFactory::~DSFactory()
{
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::dynamic_simulation::DSFactory::setYBus(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // Invoke setYBus method on all bus objects
  for (i=0; i<numBus; i++) {
    (dynamic_cast<DSBus*>(p_network->getBus(i).get()))->setYBus();
  }

  // Invoke setYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
    (dynamic_cast<DSBranch*>(p_network->getBranch(i).get()))->setYBus();
  }
}

} // namespace dynamic_simulation
} // namespace gridpack
