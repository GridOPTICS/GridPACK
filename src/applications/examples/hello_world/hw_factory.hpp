/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hw_factory.hpp
 * @author Bruce Palmer
 * @date   2014-01-28 10:32:02 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _hw_factory_h_
#define _hw_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "hw_components.hpp"

namespace gridpack {
namespace hello_world {

// This example only needs the functionality in the base factory class

class HWFactory
  : public gridpack::factory::BaseFactory<HWNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    HWFactory(boost::shared_ptr<HWNetwork> network)
      : gridpack::factory::BaseFactory<HWNetwork>(network)
    {
      p_network = network;
    }

    /**
     * Basic destructor
     */
    ~HWFactory() {}

  private:

    NetworkPtr p_network;
};

} // hello_world
} // gridpack
#endif
