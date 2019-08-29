/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_factory.hpp
 * @author Shrirang Abhyankar
 * @date   2015-01-23 08:25:26 d3g096
 * 
 * @brief  
 * Application factory definitions
 * 
 */
// -------------------------------------------------------------

#ifndef _ds_factory_h_
#define _ds_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "ds_components.hpp"

namespace gridpack {
namespace dsimplicit {

// This example only needs the functionality in the base factory class

class DSFactory
  : public gridpack::factory::BaseFactory<DSNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    DSFactory(boost::shared_ptr<DSNetwork> network)
      : gridpack::factory::BaseFactory<DSNetwork>(network)
    {
      p_network = network;
    }

    /**
     * Basic destructor
     */
    ~DSFactory() {}

  /**
   * Set the shift value provided by TS onto bus components 
   */
  void setTSshift(double);

  /**
   * Insert fault impedance 
   */
  void setfault(int,double,double);

  private:
  // NetworkPtr is a typedef for boost::shared_ptr<_network> defined in base_factory.hpp
    NetworkPtr p_network;
};

} // dsimplicit
} // gridpack
#endif
