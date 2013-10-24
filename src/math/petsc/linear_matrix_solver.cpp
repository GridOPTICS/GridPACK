// -------------------------------------------------------------
/**
 * @file   linear_matrix_solver.cpp
 * @author William A. Perkins
 * @date   2013-10-23 15:09:54 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "linear_matrix_solver.hpp"
#include "basic_linear_matrix_solver_implementation.hpp"
#include "petsc_linear_matrix_solver_impl.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearMatrixSolver
// -------------------------------------------------------------

// -------------------------------------------------------------
// LinearMatrixSolver:: constructors / destructor
// -------------------------------------------------------------
LinearMatrixSolver::LinearMatrixSolver(const Matrix& A)
  : parallel::WrappedDistributed(),
    utility::Uncopyable()
{
#ifdef PETSC_HAVE_SUPERLU
  p_impl.reset(new PetscLinearMatrixSolverImplementation(A));
#else
  p_impl.reset(new BasicLinearMatrixSolverImplementation(A));
#endif
  p_set_distributed(p_impl.get());
}


} // namespace math
} // namespace gridpack
