// -------------------------------------------------------------
/**
 * @file   linear_matrix_solver.cpp
 * @author William A. Perkins
 * @date   2013-10-21 10:06:38 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "linear_matrix_solver.hpp"
#include "basic_linear_matrix_solver_implementation.hpp"

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
  p_impl.reset(new BasicLinearMatrixSolverImplementation(A));
  p_set_distributed(p_impl.get());
}


} // namespace math
} // namespace gridpack
