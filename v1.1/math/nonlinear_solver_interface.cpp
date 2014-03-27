// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver_interface.cpp
 * @author William A. Perkins
 * @date   2013-11-08 11:51:17 d3g096
 * 
 * @brief  Implementation of NonlinearSolverInterface
 * 
 * 
 */
// -------------------------------------------------------------

#include "nonlinear_solver_interface.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NonlinearSolverInterface
// -------------------------------------------------------------

// -------------------------------------------------------------
// NonlinearSolverInterface:: constructors / destructor
// -------------------------------------------------------------
NonlinearSolverInterface::NonlinearSolverInterface()
  : parallel::WrappedDistributed(), 
    utility::WrappedConfigurable(),
    utility::Uncopyable(),
    p_impl()
{
  
}

NonlinearSolverInterface::~NonlinearSolverInterface(void)
{
}

// -------------------------------------------------------------
// NonlinearSolverInterface::p_setImpl
// -------------------------------------------------------------
void
NonlinearSolverInterface::p_setImpl(NonlinearSolverImplementation *impl)
{
  p_impl.reset(impl);
  p_setDistributed(p_impl.get());
  p_setConfigurable(p_impl.get());
}



} // namespace math
} // namespace gridpack

