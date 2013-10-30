// -------------------------------------------------------------
/**
 * @file   basic_linear_matrix_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2013-10-21 09:43:04 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "basic_linear_matrix_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class BasicLinearMatrixSolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// BasicLinearMatrixSolverImplementation:: constructors / destructor
// -------------------------------------------------------------
BasicLinearMatrixSolverImplementation::BasicLinearMatrixSolverImplementation(const Matrix& A)
  : LinearMatrixSolverImplementation(A),
    p_solver(A)
{

}

BasicLinearMatrixSolverImplementation::~BasicLinearMatrixSolverImplementation(void)
{

}

// -------------------------------------------------------------
// BasicLinearMatrixSolverImplementation::p_solve
// -------------------------------------------------------------
Matrix *
BasicLinearMatrixSolverImplementation::p_solve(const Matrix& B) const
{
  return p_solver.solve(B);
}

// -------------------------------------------------------------
// BasicLinearMatrixSolverImplementation::p_configure
// -------------------------------------------------------------
void
BasicLinearMatrixSolverImplementation::p_configure(utility::Configuration::Cursor *props)
{
  p_solver.configure(props);
}


} // namespace math
} // namespace gridpack
