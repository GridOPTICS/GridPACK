/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_factory_module.cpp
 * @author Yousu Chen, Bruce Palmer
 * @date   1/23/2015
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/hash_distr.hpp"
#include "vs_factory_module.hpp"

namespace gridpack {
namespace voltage_stability {

// State estimation factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
VSFactoryModule::VSFactoryModule(VSFactoryModule::NetworkPtr network)
  : gridpack::factory::BaseFactory<VSNetwork>(network)
{
  p_network = network;
}

/**
 * Basic destructor
 */
gridpack::voltage_stability::VSFactoryModule::~VSFactoryModule()
{
}


} // namespace voltage_stability
} // namespace gridpack
