// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   linear_matrix_solver.hpp
 * @author William A. Perkins
 * @date   2013-10-24 07:05:53 d3g096
 * 
 * @brief  Declaration of the LinearMatrixSolver class
 * 
 * 
 */
// -------------------------------------------------------------


#ifndef _linear_matrix_solver_hpp_
#define _linear_matrix_solver_hpp_

#include <boost/scoped_ptr.hpp>
#include "gridpack/math/linear_matrix_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearMatrixSolver
// -------------------------------------------------------------
/// A solver for multiple systems of linear equations
/**
 * This class implements a linear solver that solves systems with one
 * coefficient matrix and many right hand sides.  
 *
 * If the underlying math libray has a special facility for solving
 * multiple RHS, this class uses that. Otherwise, solve() is just an
 * interface to LinearSolver::solve(const Matrix& B).
 * 
 * 
 */
class LinearMatrixSolver
  : public parallel::WrappedDistributed,
    private utility::Uncopyable
{
public:

  /// Default constructor.
  LinearMatrixSolver(const Matrix& A);

  /// Destructor
  ~LinearMatrixSolver(void);

  /// Set the configuration key for this instance
  void configurationKey(const std::string& s)
  {
    p_impl->configurationKey(s);
  }

  /// Configure and do whatever is necessary to make this instance ready
  void configure(utility::Configuration::CursorPtr props)
  {
    p_impl->configure(props);
  }

  /// Solve multiple systems
  /** 
   * Each column in @c B represent a single right hand side. Each
   * column in the returned Matrix contains the solution corresponding
   * to the same column in @c B.  @c B must be \ref Matrix::Dense
   * "dense" and the returned Matrix will be \ref Matrix::Dense
   * "dense".
   * 
   * @param B RHS matrix
   * 
   * @return (dense) solution Matrix
   */
  Matrix *solve(const Matrix& B) const
  {
    return p_impl->solve(B);
  }

protected:

  /// Where the work is actually done
  boost::scoped_ptr<LinearMatrixSolverImplementation> p_impl;
};

} // namespace math
} // namespace gridpack


#endif
