/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_factory.cpp
 * @author Shrirang Abhyankar
 * @date   2015-01-26 08:25:26 d3g096
 * 
 * @brief  
 * Application factory implementation
 * 
 */
// -------------------------------------------------------------

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "ds_factory.hpp"
#include "ds_components.hpp"

namespace gridpack {
namespace dsimplicit {

  /**
   * Set shift value provided by TS onto bus components
   *
   * This is expensive. Should rather set the value of shift in the factory and then
   * have the bus components point to it. Thus, we only need to set the shift value
   * value once in the factory and it will be seen by the bus components.
   * This is not implemented currently though
   */
  void gridpack::dsimplicit::DSFactory::setTSshift(double shift) 
  {
    int numBuses = p_network->numBuses();
    int i;
    
    for(i=0; i < numBuses; i++) {
      dynamic_cast<DSBus*>(p_network->getBus(i).get())->setTSshift(shift);
    }
  }

  /**
   * Locates the faulted bus and modifies its shunt to insert the bus fault
   */
  void gridpack::dsimplicit::DSFactory::setfault(int faultbus,double Gfault,double Bfault) 
  {
    int numBuses = p_network->numBuses();
    int i,busnum;
    
    for(i=0; i < numBuses; i++) {
      gridpack::dsimplicit::DSBus *bus = dynamic_cast<gridpack::dsimplicit::DSBus*>(p_network->getBus(i).get());
      busnum = bus->getOriginalIndex();
      if(faultbus == busnum) {
	bus->addBusShunt(Gfault,Bfault);
	return;
      }
    }
  }

} // dsimplicit
} // gridpack

