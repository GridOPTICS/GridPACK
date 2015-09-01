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
 * @date   2015-08-28 15:41:56 d3g096
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
Optimizer::Optimizer()
  : OptimizerInterface()
{
  p_impl.reset(new CPlexOptimizerImplementation);
}

Optimizer::~Optimizer(void)
{
}

} // namespace optimization
} // namespace gridpack




