/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   uc_factory.hpp
 * @author 
 * @date   
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _uc_factory_h_
#define _uc_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "uc_components.hpp"

namespace gridpack {
namespace unit_commitment {

// This example only needs the functionality in the base factory class

class UCFactory
  : public gridpack::factory::BaseFactory<UCNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    UCFactory(boost::shared_ptr<UCNetwork> network)
      : gridpack::factory::BaseFactory<UCNetwork>(network)
    {
      p_network = network;
    }

    /**
     * Basic destructor
     */
    ~UCFactory() {}

  private:

    NetworkPtr p_network;
};

} // uc_commitment
} // gridpack
#endif
