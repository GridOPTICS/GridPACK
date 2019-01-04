/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_factory.cpp
 * @author Shuangshuang Jin 
 * @date   2014-12-09 14:20:51 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "ds_factory.hpp"

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
  return dummy;
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


/**
 * Initialize dynamic simulation data structures
 */
void gridpack::dynamic_simulation::DSFactory::setDSParams()
{
  int numBus = p_network->numBuses();
  int i;
  for (i=0; i<numBus; i++) {
    dynamic_cast<DSBus*>(p_network->getBus(i).get())->setDSParams();
  }
}

/**
 * Evaluate first part of a dynamic simulation step
 * @param flag false if step is not initial step
 */
void gridpack::dynamic_simulation::DSFactory::initDSStep(bool flag)
{
  int numBus = p_network->numBuses();
  int i;
  for (i=0; i<numBus; i++) {
    dynamic_cast<DSBus*>(p_network->getBus(i).get())->initDSStep(flag);
  }
}

/**
 * Evaluate predictor part of dynamic simulation step
 * @param t_inc time increment
 */
void gridpack::dynamic_simulation::DSFactory::predDSStep(double t_inc)
{
  int numBus = p_network->numBuses();
  int i;
  for (i=0; i<numBus; i++) {
    dynamic_cast<DSBus*>(p_network->getBus(i).get())->predDSStep(t_inc);
  }
}

/**
 * Evaluate corrector part of dynamic simulation step
 * @param t_inc time increment
 */
void gridpack::dynamic_simulation::DSFactory::corrDSStep(double t_inc)
{
  int numBus = p_network->numBuses();
  int i;
  for (i=0; i<numBus; i++) {
    dynamic_cast<DSBus*>(p_network->getBus(i).get())->corrDSStep(t_inc);
  }
}

} // namespace dynamic_simulation
} // namespace gridpack
