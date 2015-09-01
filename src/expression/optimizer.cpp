// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   optimizer.cpp
 * @author William A. Perkins
 * @date   2015-09-01 14:22:48 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "optimizer.hpp"
#include "cplex_optimizer_implementation.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class Optimizer
// -------------------------------------------------------------

// -------------------------------------------------------------
// Optimizer:: constructors / destructor
// -------------------------------------------------------------
Optimizer::Optimizer(const parallel::Communicator& comm)
  : OptimizerInterface(),
    parallel::WrappedDistributed(),
    utility::WrappedConfigurable(),
    p_impl()
{
  p_setImpl(new CPlexOptimizerImplementation(comm));
}

Optimizer::~Optimizer(void)
{
}

} // namespace optimization
} // namespace gridpack




