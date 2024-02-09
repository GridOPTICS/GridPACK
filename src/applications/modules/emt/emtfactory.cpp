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
void EmtFactory::setEvents(gridpack::math::RealDAESolver::EventManagerPtr eman,gridpack::mapper::GenVectorMap<EmtNetwork,gridpack::RealType,gridpack::math::RealVector> *vecmap)
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

void EmtFactory::setGlobalLocations()
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


/**
   Read events from configuration file
*/
void EmtFactory::readEvents(gridpack::utility::Configuration::CursorPtr cursor)
{
  gridpack::utility::Configuration::CursorPtr list;
  list = cursor->getCursor("Events");
  if(list) {
    gridpack::utility::Configuration::ChildElements events;
    list->children(events);
    int size = events.size();
    int idx;
    for(idx=0; idx < size; idx++) {
      if(strcmp(events[idx].name.c_str(),"BusFault") == 0) {
	int busnum;
	std::vector<int> bus_local_idx;
	EmtBus *bus;
	busnum = events[idx].cursor->get("bus",0);
	bus_local_idx = p_network->getLocalBusIndices(busnum);
	if(bus_local_idx.size()) {
	  bus = dynamic_cast<EmtBus*>(p_network->getBus(bus_local_idx[0]).get());
	
	  // Read fault parameters
	  double ton  = events[idx].cursor->get("begin",0.0);
	  double toff = events[idx].cursor->get("end",0.0);
	  std::string faulttype = events[idx].cursor->get("type","SLG");
	  std::string faultphases;
	  if(faulttype == "SLG") {
	    faultphases = events[idx].cursor->get("phases","A");
	  } else if(faulttype == "ThreePhase") {
	    faultphases = "ABC";
	  }
	  double Ron = events[idx].cursor->get("Ron",1e-2);
	  double Rgnd = events[idx].cursor->get("Rgnd",1e-3);

	  bus->setFault(ton,toff,faulttype,faultphases,Ron,Rgnd);
	}
      }
    }
  }
}
