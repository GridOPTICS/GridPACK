/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   emtfactory.hpp
 * 
 * @brief  
 * Application factory definitions
 * 
 */
// -------------------------------------------------------------

#ifndef _emtfactory_h_
#define _emtfactory_h_


#include <boost/smart_ptr/shared_ptr.hpp>
#include <gridpack/include/gridpack.hpp>
#include <gridpack/applications/modules/emt/emtnetwork.hpp>
#include <gridpack/math/dae_solver.hpp>

// This example only needs the functionality in the base factory class

class EmtFactory
  : public gridpack::factory::BaseFactory<EmtNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    EmtFactory(boost::shared_ptr<EmtNetwork> network)
      : gridpack::factory::BaseFactory<EmtNetwork>(network)
    {
      p_network = network;
    }

    /**
     * Basic destructor
     */
    ~EmtFactory() {}

  /**
   * setTime - set the current time onto bus and branch components
   */
  void setTime(double);
  
  /**
   * Set the shift value provided by TS onto bus and branch components 
   */
  void setTSshift(double);

  /**
     Read events from configuration file
  */
  void readEvents(gridpack::utility::Configuration::CursorPtr);

  /** 
    Set up components
  */
  void setup(void);

  /**
     Set the global location for the first variable for a bus
     and branch in the solution vector. This location is used when
     setting values in the matrix
  */
  void setGlobalLocations(void);
  
  /**
   * Set events 
   */
  void setEvents(gridpack::math::RealDAESolver::EventManagerPtr,
		 gridpack::mapper::GenVectorMap<EmtNetwork,gridpack::RealType,gridpack::math::RealVector>*);

  /** 
   * Reset flags after event is handled
   */
  void resetEventFlags();
  private:
  // NetworkPtr is a typedef for boost::shared_ptr<_network> defined in base_factory.hpp
    NetworkPtr p_network;
};

#endif
