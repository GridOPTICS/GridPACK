/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

#include "dsf_app_module.hpp"
#include <iostream>
#include <string>
#include <vector>
#include "gridpack/utilities/string_utils.hpp"

/**
   Set an event for the dynamic simulation
*/
void gridpack::dynamic_simulation::DSFullApp::
setEvent(gridpack::dynamic_simulation::Event event)
{
  p_events.clear();
  p_events.push_back(event);
}

/**
 * Read events starting at the config cursor  and form a list of events
 * @param  cursor to Events in input.xml file
 * @return a list of events
 * Note : The configuration needs to be already set
 : setNetwork method
 *
 */
std::vector<gridpack::dynamic_simulation::Event>
gridpack::dynamic_simulation::DSFullApp::getEvents(gridpack::utility::Configuration::CursorPtr cursor)
{
  gridpack::utility::Configuration::CursorPtr list;
  list = cursor->getCursor("Events");
  gridpack::utility::Configuration::ChildElements events;
  std::vector<gridpack::dynamic_simulation::Event> ret;
  if (list) {
    list->children(events);
    int size = events.size();
    int idx;
    // Parse events
    for (idx=0; idx<size; idx++) {
      gridpack::dynamic_simulation::Event event;
      event.isBus = false;
      event.isBusFault = false;
      event.isLine = false;
      event.isLineStatus = false;
      event.isGenerator = false;
      event.isGenStatus = false;
      if(strcmp(events[idx].name.c_str(),"faultEvent") == 0) {
        // This is kept for backward compatibility only and should be removed
        // once all the dynamic simulation input.xml files are updated
        // with the new format
        event.start = events[idx].cursor->get("beginFault",0.0);
        event.end = events[idx].cursor->get("endFault",0.0);
        std::string line_ends = events[idx].cursor->get("faultBranch","0 0");
        if(line_ends.c_str() != NULL) {
          sscanf(line_ends.c_str(),"%d %d",&event.from_idx, &event.to_idx);
        }
        event.bus_idx = event.from_idx;
        event.tag = events[idx].cursor->get("id","1");

        event.Gfault = 0.0;
        event.Bfault = 99999.0;

        event.isLine = true;
        event.isBusFault = true;
        ret.push_back(event);
      } else if(strcmp(events[idx].name.c_str(),"BusFault") == 0) {
        // Bus fault event
        event.start = events[idx].cursor->get("begin",0.0);
        event.end = events[idx].cursor->get("end",0.0);
        event.bus_idx = events[idx].cursor->get("bus",0);

        std::string Yfault = events[idx].cursor->get("yfault","0 0");
        if(Yfault.c_str() != NULL) {
          sscanf(Yfault.c_str(),"%lf %lf",&event.Gfault, &event.Bfault);
        }

        event.isBusFault = true;
        ret.push_back(event);
      } else if(strcmp(events[idx].name.c_str(),"LineStatus") == 0) {
        // Line status change event
        event.time = events[idx].cursor->get("time",0.0);
        std::string line = events[idx].cursor->get("line","0 0");
        if(line.c_str() != NULL) {
          sscanf(line.c_str(),"%d %d",&event.from_idx, &event.to_idx);
        }
        event.tag = events[idx].cursor->get("id","1");
        event.status = events[idx].cursor->get("status",1);

        event.isLineStatus = true;
        ret.push_back(event);
      } else if(strcmp(events[idx].name.c_str(),"GenStatus") == 0) {
        // Gen status change event
        event.time = events[idx].cursor->get("time",0.0);
        event.bus_idx = events[idx].cursor->get("bus",0);
        event.tag = events[idx].cursor->get("id","1");
        event.status = events[idx].cursor->get("status",1);


        event.isGenStatus = true;
        ret.push_back(event);
      }
    }
  }
  return ret;
}

/**
 * Read events set in the config file and form a list of events
 * @return a list of events
 * Note : The configuration needs to be already set
 : setNetwork method
 *
 */
std::vector<gridpack::dynamic_simulation::Event>
gridpack::dynamic_simulation::DSFullApp::getEvents(gridpack::utility::Configuration::CursorPtr cursor)
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");

  return getEvents(cursor);
}


/**
 * handleEvents - handle any events
**/
void gridpack::dynamic_simulation::DSFullApp::handleEvents()
{
  int nevents = p_events.size();
  int i;

  for(i = 0; i < nevents; i++) {
    gridpack::dynamic_simulation::Event event = p_events[i];

    if(event.isBusFault) {
      if(fabs(event.start - p_current_time) < 1e-6) {
	/* Fault start */
	std::vector<int> bus_internal_idx;
	gridpack::dynamic_simulation::DSFullBus *bus;
	bus_internal_idx = p_network->getLocalBusIndices(event.bus_idx);
	if(bus_internal_idx.size()) {
	  bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(p_network->getBus(bus_internal_idx[0]).get());
	  bus->setFault(event.Gfault,event.Bfault);
	}

	// Update Ybus
	p_factory->setMode(BUSFAULTON);
	ybusMap_sptr->incrementMatrix(ybus);
	
      } else if(fabs(event.end - p_current_time) < 1e-6) {
	/* Fault end */
	std::vector<int> bus_internal_idx;
	gridpack::dynamic_simulation::DSFullBus *bus;
	bus_internal_idx = p_network->getLocalBusIndices(event.bus_idx);
	if(bus_internal_idx.size()) {
	  bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(p_network->getBus(bus_internal_idx[0]).get());
	}
	
	// Update Ybus
	p_factory->setMode(BUSFAULTOFF);
	ybusMap_sptr->incrementMatrix(ybus);
      }
    } else if(event.isLineStatus) {
      if(fabs(event.time - p_current_time) < 1e-6) {
	setLineStatus(event.from_idx,event.to_idx,event.tag,event.status);
	p_factory->setMode(LINESTATUSCHANGE);
	ybusMap_sptr->incrementMatrix(ybus);
      }
    } else if(event.isGenStatus) {
      if(fabs(event.time - p_current_time) < 1e-6) {
	/* Generator status change */
	setGenStatus(event.bus_idx,event.tag,event.status);
	p_factory->setMode(GENSTATUSCHANGE);
	ybusMap_sptr->incrementMatrix(ybus);
      }
    }
  }
}

/**
 * Read events set in the config file and form a list of events
 * @return a list of events
 * Note : The configuration needs to be already set
 : setNetwork method
 *
 */
std::vector<gridpack::dynamic_simulation::Event>
gridpack::dynamic_simulation::DSFullApp::getEvents()
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");

  return getEvents(cursor);
}


/**
 * Read events set in the config file and form a list of events
 * @return a list of events
 * Note : The configuration needs to be already set
 : setNetwork method
 *
 */
std::vector<gridpack::dynamic_simulation::Event>
gridpack::dynamic_simulation::DSFullApp::getEvents()
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Dynamic_simulation");

  return getEvents(cursor);
}


/**
 * handleEvents - handle any events
**/
void gridpack::dynamic_simulation::DSFullApp::handleEvents()
{
  int nevents = p_events.size();
  int i;

  for(i = 0; i < nevents; i++) {
    gridpack::dynamic_simulation::Event event = p_events[i];

    if(event.isBusFault) {
      if(fabs(event.start - p_current_time) < 1e-6) {
	/* Fault start */
	std::vector<int> bus_internal_idx;
	gridpack::dynamic_simulation::DSFullBus *bus;
	bus_internal_idx = p_network->getLocalBusIndices(event.bus_idx);
	if(bus_internal_idx.size()) {
	  bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(p_network->getBus(bus_internal_idx[0]).get());
	  bus->setFault(event.Gfault,event.Bfault);
	}

	// Update Ybus
	p_factory->setMode(BUSFAULTON);
	ybusMap_sptr->incrementMatrix(ybus);
	
      } else if(fabs(event.end - p_current_time) < 1e-6) {
	/* Fault end */
	std::vector<int> bus_internal_idx;
	gridpack::dynamic_simulation::DSFullBus *bus;
	bus_internal_idx = p_network->getLocalBusIndices(event.bus_idx);
	if(bus_internal_idx.size()) {
	  bus = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(p_network->getBus(bus_internal_idx[0]).get());
	}
	
	// Update Ybus
	p_factory->setMode(BUSFAULTOFF);
	ybusMap_sptr->incrementMatrix(ybus);
      }
    } else if(event.isLineStatus) {
      if(fabs(event.time - p_current_time) < 1e-6) {
	setLineStatus(event.from_idx,event.to_idx,event.tag,event.status);
	p_factory->setMode(LINESTATUSCHANGE);
	ybusMap_sptr->incrementMatrix(ybus);
      }
    } else if(event.isGenStatus) {
      if(fabs(event.time - p_current_time) < 1e-6) {
	/* Generator status change */
	setGenStatus(event.bus_idx,event.tag,event.status);
	p_factory->setMode(GENSTATUSCHANGE);
	ybusMap_sptr->incrementMatrix(ybus);
      }
    }
  }
}
