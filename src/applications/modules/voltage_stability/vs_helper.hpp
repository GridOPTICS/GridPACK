// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   vs_helper.cpp
 * @author William A. Perkins
 * @date   2014-11-25 07:12:27 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <ga++.h>

#include "vs_factory_module.hpp"


namespace gridpack {
namespace voltage_stability {

/**
 * A helper functor for the powerflow solver.
 * This is a utility functor that provides functions that build the
 * Jacobian and RHS from the network.  
 * 
 */
struct VSSolverHelper 
  : private utility::Uncopyable
{
  
  // pointer to factory
  boost::shared_ptr<gridpack::powerflow::PFFactoryModule> p_factory;

  // pointer to factory
  boost::shared_ptr<VSFactoryModule> v_factory;
  // The powerflow factory 
  //boost::shared_ptr<VSFactoryModule> p_factory;
  
  // pointer to network
  boost::shared_ptr<gridpack::powerflow::PFNetwork> p_network;
  // The powerflow network controlled by ::p_factory
  //boost::shared_ptr<VSNetwork> p_network;

  // A place to build/store the Jacobian
  boost::shared_ptr<math::RealMatrix> J;


  };
}}