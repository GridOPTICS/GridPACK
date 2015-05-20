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
 * @date   2015-03-06 12:17:02 d3g096
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
#include "gridpack/math/linear_matrix_solver_interface.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearMatrixSolverImplementation
// -------------------------------------------------------------
/// 
template <typename T, typename I>
class LinearMatrixSolverImplementation 
  : public BaseLinearMatrixSolverInterface<T, I>,
    public parallel::Distributed,
    public utility::Configurable,
    private utility::Uncopyable
{
public:

  typedef typename BaseLinearMatrixSolverInterface<T, I>::MatrixType MatrixType;

  /// Default constructor.
  LinearMatrixSolverImplementation(const MatrixType& A)
    : BaseLinearMatrixSolverInterface<T, I>(),
      parallel::Distributed(A.communicator()),
      utility::Configurable(),
      utility::Uncopyable(),
      p_A(A.clone())
  {
    configurationKey("LinearMatrixSolver");
  }


  /// Destructor
  ~LinearMatrixSolverImplementation(void)
  {}

  /// Solve w/ the specified RHS Matrix, return (dense) Matrix
  MatrixType *solve(const MatrixType& B) const
  {
    return this->p_solve(B);
  }

protected:

  /// The coefficient matrix (may not need to remember)
  boost::scoped_ptr<MatrixType> p_A;
  
  /// Solve w/ the specified RHS Matrix (specialized)
  virtual MatrixType *p_solve(const MatrixType& B) const = 0;

};


} // namespace math
} // namespace gridpack


#endif
