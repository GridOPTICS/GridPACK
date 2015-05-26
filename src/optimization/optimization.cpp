/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_factory.cpp
 * @author Bruce Palmer
 * @date   2013-08-08 10:06:51 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/component/base_component.hpp"
#include "gridpack/optimization/optimization.hpp"


namespace gridpack {

// just to force instantiation of a base optimizer<>
static 
void a_function(void) {

  boost::mpi::communicator world;
  typedef network::BaseNetwork<component::BaseBusComponent,
                               component::BaseBranchComponent> ANetwork;
  typedef optimization::Optimizer<ANetwork> AOptimizer;
  AOptimizer::NetworkPtr a_network(new ANetwork(world));
  AOptimizer a_optimizer(a_network);
}

} // namespace gridpack
