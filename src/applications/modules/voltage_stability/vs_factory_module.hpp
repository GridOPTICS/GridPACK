/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_factory_module.hpp
 * @author Yousu Chen, Bruce Palmer 
 * @date   1/23/2015
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _vs_factory_module_h_
#define _vs_factory_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/applications/components/vs_matrix/vs_components.hpp"
/// The type of network used in the powerflow application

  
  
namespace gridpack {
namespace voltage_stability {

typedef gridpack::network::BaseNetwork<VSBus, VSBranch > VSNetwork;

class VSFactoryModule
  : public gridpack::factory::BaseFactory<VSNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    VSFactoryModule(NetworkPtr network);

    /**
     * Basic destructor
     */
    ~VSFactoryModule();


  private:

    NetworkPtr p_network;
};

} // voltage_stability
} // gridpack
#endif
