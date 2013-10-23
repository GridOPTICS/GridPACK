// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   linear_matrix_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2013-10-21 10:22:12 d3g096
 * 
 * @brief  Declaration of LinearMatrixSolverImplementation class
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _linear_matrix_solver_implementation_hpp_
#define _linear_matrix_solver_implementation_hpp_

#include "gridpack/parallel/distributed.hpp"
#include "gridpack/configuration/configurable.hpp"
#include "gridpack/utilities/uncopyable.hpp"
#include "gridpack/math/matrix.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearMatrixSolverImplementation
// -------------------------------------------------------------
/// 
class LinearMatrixSolverImplementation 
  : public parallel::Distributed,
    public utility::Configurable,
    private utility::Uncopyable
{
public:

  /// Default constructor.
  LinearMatrixSolverImplementation(const Matrix& A);

  /// Destructor
  ~LinearMatrixSolverImplementation(void);

  /// Solve w/ the specified RHS Matrix, return (dense) Matrix
  Matrix *solve(const Matrix& B) const
  {
    return this->p_solve(B);
  }

protected:

  /// The coefficient matrix (may not need to remember)
  boost::scoped_ptr<Matrix> p_A;
  
  /// Solve w/ the specified RHS Matrix (specialized)
  virtual Matrix *p_solve(const Matrix& B) const = 0;

};


} // namespace math
} // namespace gridpack


#endif
