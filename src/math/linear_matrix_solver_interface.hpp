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
 * @file   linear_matrix_solver_interface.hpp
 * @author William A. Perkins
 * @date   2015-03-06 11:51:39 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _linear_matrix_solver_interface_hpp_
#define _linear_matrix_solver_interface_hpp_

#include "gridpack/parallel/distributed.hpp"
#include "gridpack/configuration/configurable.hpp"
#include "gridpack/utilities/uncopyable.hpp"
#include "gridpack/math/matrix.hpp"


namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class BaseLinearMatrixSolverInterface
// -------------------------------------------------------------
template <typename T, typename I = int>
class BaseLinearMatrixSolverInterface 
{
public:

  typedef MatrixT<T, I> MatrixType;

  /// Default constructor.
  BaseLinearMatrixSolverInterface(void)
  {}


  /// Destructor
  virtual ~BaseLinearMatrixSolverInterface(void)
  {}

  /// Solve w/ the specified RHS Matrix, return (dense) Matrix
  MatrixType *solve(const MatrixType& B) const
  {
    return this->p_solve(B);
  }

protected:

  /// Solve w/ the specified RHS Matrix (specialized)
  virtual MatrixType *p_solve(const MatrixType& B) const = 0;

};




} // namespace math
} // namespace gridpack

#endif
