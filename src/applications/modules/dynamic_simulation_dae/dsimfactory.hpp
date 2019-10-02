/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsimfactory.hpp
 * @author Shrirang Abhyankar
 * @date   02/10/19
 * 
 * @brief  
 * Application factory definitions
 * 
 */
// -------------------------------------------------------------

#ifndef _dsimfactory_h_
#define _dsimfactory_h_


#include <boost/smart_ptr/shared_ptr.hpp>
#include <gridpack/include/gridpack.hpp>
#include <dsimnetwork.hpp>

// This example only needs the functionality in the base factory class

class DSimFactory
  : public gridpack::factory::BaseFactory<DSimNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    DSimFactory(boost::shared_ptr<DSimNetwork> network)
      : gridpack::factory::BaseFactory<DSimNetwork>(network)
    {
      p_network = network;
    }

    /**
     * Basic destructor
     */
    ~DSimFactory() {}

  /**
   * Set the shift value provided by TS onto bus components 
   */
  void setTSshift(double);

  /**
   * Insert fault impedance 
   */
  void setfault(int,double,double);

  /** 
    Initialize components
  */
  void initialize(void);

  private:
  // NetworkPtr is a typedef for boost::shared_ptr<_network> defined in base_factory.hpp
    NetworkPtr p_network;
};

#endif
