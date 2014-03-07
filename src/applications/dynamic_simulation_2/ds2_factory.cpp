/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds2_factory.cpp
 * @author Shuangshuang Jin 
 * @date   March 07, 2014
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
#include "gridpack/applications/dynamic_simulation_2/ds2_components.hpp"
#include "gridpack/applications/dynamic_simulation_2/ds2_factory.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"

namespace gridpack {
namespace dynamic_simulation {

// Dynamic simulation factory class implementations

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

/**
 * Get the updating factor for posfy11 stage ybus
 */
gridpack::ComplexType
gridpack::dynamic_simulation::DSFactory::setFactor(int sw2_2, int sw3_2)
{
  gridpack::ComplexType dummy(-999.0, -999.0);

  int numBranch = p_network->numBranches();
  int i;

  // Invoke getPosfy11YbusUpdateFactor method on all branch objects
  for (i=0; i<numBranch; i++) {
    gridpack::ComplexType ret = (dynamic_cast<DSBranch*>(p_network->getBranch(i).get()))
      ->getPosfy11YbusUpdateFactor(sw2_2, sw3_2);
    if (ret != dummy) {
      return ret;
    }
  }
}

/**
 * Apply an event to all branches in the system
 * @param event a struct describing a fault
 */
void gridpack::dynamic_simulation::DSFactory::setEvent(const
    gridpack::dynamic_simulation::DSBranch::Event &event)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;
  for (i=0; i<numBus; i++) {
    dynamic_cast<DSBus*>(p_network->getBus(i).get())->clearEvent();
  }
  for (i=0; i<numBranch; i++) {
    dynamic_cast<DSBranch*>(p_network->getBranch(i).get())->setEvent(event);
  }
}

/**
 * Check network to see if there is a process with no generators
 * @return true if all processors have at least on generator
 */
bool gridpack::dynamic_simulation::DSFactory::checkGen(void)
{
  int numBus = p_network->numBuses();
  int i, count;
  count = 0;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      count += dynamic_cast<DSBus*>(p_network->getBus(i).get())->getNumGen();
    }
  }
  printf("p[%d] number of generators: %d\n",p_network->communicator().rank(),count);
  int iok = 0;
  if (count > 0) iok = 1;
  int ok;
  boost::mpi::all_reduce(p_network->communicator(),iok,ok,std::plus<int>());
  bool ret = false;
  if (ok == p_network->communicator().size()) ret = true;
  return ret;
}

} // namespace dynamic_simulation
} // namespace gridpack
