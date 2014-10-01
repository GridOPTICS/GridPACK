// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dae_solver.cpp
 * @author William A. Perkins
 * @date   2013-11-13 09:49:20 d3g096
 * 
 * @brief  
 * 
 * 
 */

#include "dae_solver.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class DAESolver
// -------------------------------------------------------------

// -------------------------------------------------------------
// DAESolver:: constructors / destructor
// -------------------------------------------------------------
DAESolver::~DAESolver(void)
{
}

// -------------------------------------------------------------
// DAESolver::p_setImpl
// -------------------------------------------------------------
void
DAESolver::p_setImpl(DAESolverImplementation *impl)
{
  p_impl.reset(impl);
  p_setDistributed(p_impl.get());
  p_setConfigurable(p_impl.get());
}


} // namespace math
} // namespace gridpack
