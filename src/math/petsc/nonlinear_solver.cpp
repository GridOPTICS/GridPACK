/**
 * @file   nonlinear_solver.cpp
 * @author William A. Perkins
 * @date   2013-09-06 11:17:33 d3g096
 * 
 * @brief  PETSc-specific implementation of NonlinearSolver
 * 
 * 
 */

#include "nonlinear_solver.hpp"
#include "petsc/petsc_nonlinear_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NonlinearSolver
// -------------------------------------------------------------

// -------------------------------------------------------------
// NonlinearSolver:: constructors / destructor
// -------------------------------------------------------------
NonlinearSolver::NonlinearSolver(const parallel::Communicator& comm,
                                 const int& local_size,
                                 JacobianBuilder form_jacobian,
                                 FunctionBuilder form_function)
  : parallel::WrappedDistributed(),
    utility::Uncopyable(), 
    p_impl(new PetscNonlinearSolverImplementation(comm, local_size, form_jacobian, form_function))
{
  p_set_distributed(p_impl.get());
}

} // namespace math
} // namespace gridpack
