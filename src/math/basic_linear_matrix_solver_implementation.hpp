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
 * @date   2015-03-06 12:17:58 d3g096
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
template <typename T, typename I = int>
class BasicLinearMatrixSolverImplementation 
  : public LinearMatrixSolverImplementation<T, I>
{
public:

  typedef typename BaseLinearMatrixSolverInterface<T, I>::MatrixType MatrixType;

  /// Default constructor.
  BasicLinearMatrixSolverImplementation(MatrixType& A)
    : LinearMatrixSolverImplementation<T, I>(A),
      p_solver(A)
  {
  }


  /// Destructor
  ~BasicLinearMatrixSolverImplementation(void)
  {}

protected:

  /// The linear solver instance used for this
  LinearSolverT<T, I> p_solver;

  /// Solve w/ the specified RHS Matrix (specialized)
  MatrixType *p_solve(const MatrixType& B) const
  {
    return p_solver.solve(B);
  }

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props)
  {
    p_solver.configure(props);
  }

};


} // namespace math
} // namespace gridpack
#endif
