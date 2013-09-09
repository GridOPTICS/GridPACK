/**
 * @file   newton_raphson_solver.cpp
 * @author William A. Perkins
 * @date   2013-09-09 11:37:38 d3g096
 * 
 * @brief  
 * 
 * 
 */

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
// NewtonRaphsonSolver
// -------------------------------------------------------------


} // namespace math
} // namespace gridpack
