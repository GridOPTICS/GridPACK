// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   linear_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2014-10-22 09:13:34 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  1, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _linear_solver_implementation_hpp_
#define _linear_solver_implementation_hpp_

#include <boost/scoped_ptr.hpp>
#include <gridpack/math/linear_solver_interface.hpp>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/configuration/configurable.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearSolverImplementation
// -------------------------------------------------------------
class LinearSolverImplementation 
  : public parallel::Distributed,
    public utility::Configurable,
    private utility::Uncopyable,
    public BaseLinearSolverInterface
{
public:
  
  /// Default constructor.
  LinearSolverImplementation(const parallel::Communicator& comm);

  /// Destructor
  ~LinearSolverImplementation(void);
  
protected:

  /// The solution residual norm tolerance
  /**
   * This is the absolute solution tolerance. The linear system is
   * considered solved when the solution residual norm is below this value.
   * This is just handed to the underlying math library, so 
   * 
   */
  double p_solutionTolerance;

  /// The relative solution tolerance.
  /**
   * The linear system is considered solved when the solution residual
   * norm is @em reduced by the specified amount.
   * 
   */
  double p_relativeTolerance;

  /// The maximum number of iterations to perform
  int p_maxIterations;

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props);

  /// Get the solution tolerance (specialized)
  double p_tolerance(void) const;

  /// Set the solver tolerance (specialized)
  void p_tolerance(const double& tol);

  /// Get the maximum iterations (specialized)
  int p_maximumIterations(void) const;

  /// Set the maximum solution iterations
  void p_maximumIterations(const int& n);

  /// Solve multiple systems w/ each column of the Matrix a single RHS
  Matrix *p_solve(const Matrix& B) const;
};

} // namespace math
} // namespace gridpack



#endif
