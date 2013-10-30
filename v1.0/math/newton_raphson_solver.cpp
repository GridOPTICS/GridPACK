// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   newton_raphson_solver.cpp
 * @author William A. Perkins
 * @date   2013-10-09 12:24:32 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/assert.hpp>
#include "newton_raphson_solver.hpp"
#include "newton_raphson_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NewtonRaphsonSolver
// -------------------------------------------------------------

// -------------------------------------------------------------
// NewtonRaphsonSolver:: constructors / destructor
// -------------------------------------------------------------
NewtonRaphsonSolver::NewtonRaphsonSolver(const parallel::Communicator& comm,
                                         const int& local_size,
                                         JacobianBuilder form_jacobian,
                                         FunctionBuilder form_function)
  : NonlinearSolverInterface()
{
  p_set_impl(new NewtonRaphsonSolverImplementation(comm, local_size,
                                                   form_jacobian,
                                                   form_function));
  p_set_distributed(p_impl.get());
}

NewtonRaphsonSolver::~NewtonRaphsonSolver(void)
{
}

// -------------------------------------------------------------
// NewtonRaphsonSolver::tolerance
// -------------------------------------------------------------
double
NewtonRaphsonSolver::tolerance(void) const
{
  NewtonRaphsonSolverImplementation *s = 
    dynamic_cast<NewtonRaphsonSolverImplementation *>(p_impl.get());
  BOOST_ASSERT(s != NULL);
  return s->tolerance();
}

void
NewtonRaphsonSolver::tolerance(const double& tol)
{
  NewtonRaphsonSolverImplementation *s = 
    dynamic_cast<NewtonRaphsonSolverImplementation *>(p_impl.get());
  BOOST_ASSERT(s != NULL);
  s->tolerance(tol);
}

// -------------------------------------------------------------
// NewtonRaphsonSolver::maximum_iterations
// -------------------------------------------------------------
int
NewtonRaphsonSolver::maximum_iterations(void) const
{
  NewtonRaphsonSolverImplementation *s = 
    dynamic_cast<NewtonRaphsonSolverImplementation *>(p_impl.get());
  BOOST_ASSERT(s != NULL);
  return s->maximum_iterations();
}

void
NewtonRaphsonSolver::maximum_iterations(const int& n)
{
  NewtonRaphsonSolverImplementation *s = 
    dynamic_cast<NewtonRaphsonSolverImplementation *>(p_impl.get());
  BOOST_ASSERT(s != NULL);
  s->maximum_iterations(n);
}

  


} // namespace math
} // namespace gridpack
