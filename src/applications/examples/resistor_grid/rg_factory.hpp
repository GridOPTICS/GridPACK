/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   rg_factory.hpp
 * @author Bruce Palmer
 * @date   2014-02-05 08:25:26 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _rg_factory_h_
#define _rg_factory_h_

#include "rg_components.hpp"

namespace gridpack {
namespace resistor_grid {

// This example only needs the functionality in the base factory class

class RGFactory
  : public gridpack::factory::BaseFactory<RGNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    RGFactory(boost::shared_ptr<RGNetwork> network)
      : gridpack::factory::BaseFactory<RGNetwork>(network)
    {
    }

    /**
     * Basic destructor
     */
    ~RGFactory() {}
};

} // resistor_grid
} // gridpack
#endif
