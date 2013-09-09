/**
 * @file   nonlinear_solver_interface.cpp
 * @author William A. Perkins
 * @date   2013-09-09 10:21:02 d3g096
 * 
 * @brief  Implementation of NonlinearSolverInterface
 * 
 * 
 */

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
    utility::Uncopyable(),
    p_impl()
{
  
}

NonlinearSolverInterface::~NonlinearSolverInterface(void)
{
}

// -------------------------------------------------------------
// NonlinearSolverInterface::p_set_impl
// -------------------------------------------------------------
void
NonlinearSolverInterface::p_set_impl(NonlinearSolverImplementation *impl)
{
  p_impl.reset(impl);
  p_set_distributed(p_impl.get());
}



} // namespace math
} // namespace gridpack

