// -------------------------------------------------------------
/**
 * @file   linear_matrix_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2013-10-21 10:28:33 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "linear_matrix_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearMatrixSolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// LinearMatrixSolverImplementation:: constructors / destructor
// -------------------------------------------------------------
LinearMatrixSolverImplementation::LinearMatrixSolverImplementation(const Matrix& A)
  : parallel::Distributed(A.communicator()),
    utility::Configurable(),
    utility::Uncopyable(),
    p_A(A.clone())
{
  configurationKey("LinearMatrixSolver");
}

LinearMatrixSolverImplementation::~LinearMatrixSolverImplementation(void)
{
}

} // namespace math
} // namespace gridpack
