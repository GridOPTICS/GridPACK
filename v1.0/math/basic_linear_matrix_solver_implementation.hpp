// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   basic_linear_matrix_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2013-10-21 09:41:36 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _basic_linear_matrix_solver_implementation_hpp_
#define _basic_linear_matrix_solver_implementation_hpp_

#include "linear_solver.hpp"
#include "linear_matrix_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class BasicLinearMatrixSolverImplementation
// -------------------------------------------------------------
/// 
class BasicLinearMatrixSolverImplementation 
  : public LinearMatrixSolverImplementation
{
public:

  /// Default constructor.
  BasicLinearMatrixSolverImplementation(const Matrix& A);

  /// Destructor
  ~BasicLinearMatrixSolverImplementation(void);

protected:

  /// The linear solver instance used for this
  LinearSolver p_solver;

  /// Solve w/ the specified RHS Matrix (specialized)
  Matrix *p_solve(const Matrix& B) const;

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::Cursor *props);

};


} // namespace math
} // namespace gridpack
#endif
