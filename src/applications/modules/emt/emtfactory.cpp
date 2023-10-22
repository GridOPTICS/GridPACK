/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   emtfactory.cpp
 * 
 * @brief  
 * Application factory implementation
 * 
 */
// -------------------------------------------------------------

#include <boost/smart_ptr/shared_ptr.hpp>
#include <gridpack/include/gridpack.hpp>
#include <emtnetwork.hpp>
#include <emtfactory.hpp>
#include <gridpack/math/dae_solver.hpp>

/**
 * Set the current time provided by TS onto bus and branch components
 *
 * This is expensive. Should rather set the value of shift in the factory and then
 * have the bus components point to it. Thus, we only need to set the shift value
 * value once in the factory and it will be seen by the bus components.
 * This is not implemented currently though
 */
void EmtFactory::setTime(double time) 
{
  int numBuses = p_network->numBuses();
  int numBranches = p_network->numBranches();
  int i;
  
  for(i=0; i < numBuses; i++) {
    dynamic_cast<EmtBus*>(p_network->getBus(i).get())->setTime(time); 
  }

  for(i=0; i < numBranches; i++) {
    dynamic_cast<EmtBranch*>(p_network->getBranch(i).get())->setTime(time); 
  }
}

/**
 * Set shift value provided by TS onto bus and branch components
 *
 * This is expensive. Should rather set the value of shift in the factory and then
 * have the bus components point to it. Thus, we only need to set the shift value
 * value once in the factory and it will be seen by the bus components.
 * This is not implemented currently though
 */
void EmtFactory::setTSshift(double shift) 
{
  int numBuses = p_network->numBuses();
  int numBranches = p_network->numBranches();
  int i;
  
  for(i=0; i < numBuses; i++) {
    dynamic_cast<EmtBus*>(p_network->getBus(i).get())->setTSshift(shift); 
  }

  for(i=0; i < numBranches; i++) {
    dynamic_cast<EmtBranch*>(p_network->getBranch(i).get())->setTSshift(shift); 
  }
}

/** 
 * Add events from buses and branches to the event manager 
*/
void EmtFactory::setEvents(gridpack::math::DAESolver::EventManagerPtr eman,gridpack::mapper::GenVectorMap<EmtNetwork> *vecmap)
{
  int numBuses = p_network->numBuses();
  int numBranches = p_network->numBranches();
  int i;
  int offset,size;
  EmtBus *bus;
  EmtBranch *branch;
  
  for(i=0; i < numBuses; i++) {
    bus = dynamic_cast<EmtBus*>(p_network->getBus(i).get());
    bus->setEvent(eman);

    // This should be moved to a separate function
    vecmap->getLocalBusOffset(i,&offset,&size);
    bus->setLocalOffset(offset);
  }

  for(i=0; i < numBranches; i++) {
    branch = dynamic_cast<EmtBranch*>(p_network->getBranch(i).get());

    // This should be moved to a separate function
    vecmap->getLocalBranchOffset(i,&offset,&size);
    branch->setLocalOffset(offset);
  }
}

/** 
 * Reset flags after event is handled
*/
void EmtFactory::resetEventFlags()
{
  int numBuses = p_network->numBuses();
  int i;
  int offset,size;
  EmtBus *bus;
  
  for(i=0; i < numBuses; i++) {
    bus = dynamic_cast<EmtBus*>(p_network->getBus(i).get());
    bus->resetEventFlags();
  }
}

void EmtFactory::setup(void) 
{
  int numBuses = p_network->numBuses();
  int numBranches = p_network->numBranches();
  int i;
  bool isactive;
  int rank = p_network->communicator().rank();

  for(i=0; i < numBranches; i++) {
    int idxf,idxt;
    isactive = p_network->getActiveBranch(i);
    EmtBranch* emtbranch = dynamic_cast<EmtBranch*>(p_network->getBranch(i).get());
    emtbranch->setGhostStatus(!isactive);
    emtbranch->setRank(rank);
    emtbranch->setup();
  }

  
  for(i=0; i < numBuses; i++) {
    isactive = p_network->getActiveBus(i);
    EmtBus* emtbus = dynamic_cast<EmtBus*>(p_network->getBus(i).get());
    emtbus->setGhostStatus(!isactive);
    emtbus->setRank(rank);
    emtbus->setup();
  }

}

void EmtFactory::setGlobalIndices()
{
  int numBuses = p_network->numBuses();
  int numBranches = p_network->numBranches();
  int i;

  for(i=0; i < numBranches; i++) {
    EmtBranch* emtbranch = dynamic_cast<EmtBranch*>(p_network->getBranch(i).get());
    emtbranch->setGlobalLocation();
  }

  
  for(i=0; i < numBuses; i++) {
    EmtBus* emtbus = dynamic_cast<EmtBus*>(p_network->getBus(i).get());
    emtbus->setGlobalLocation();
  }
}


