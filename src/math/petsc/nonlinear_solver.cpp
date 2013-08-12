/**
 * @file   nonlinear_solver.cpp
 * @author William A. Perkins
 * @date   2013-08-12 12:59:05 d3g096
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
  : utility::Uncopyable(), 
    p_impl(new PetscNonlinearSolverImplementation(comm, local_size, form_jacobian, form_function))
{
  
}

} // namespace math
} // namespace gridpack
