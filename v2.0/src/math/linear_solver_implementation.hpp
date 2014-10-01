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
 * @date   2014-01-09 12:09:55 d3g096
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
#include <gridpack/math/matrix.hpp>
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
    private utility::Uncopyable
{
public:
  
  /// Default constructor.
  LinearSolverImplementation(const parallel::Communicator& comm);

  /// Destructor
  ~LinearSolverImplementation(void);
  
  /// Solve w/ the specified RHS, put result in specified vector
  void solve(const Vector& b, Vector& x) const
  {
    this->p_solve(b, x);
  }

  /// Solve multiple systems w/ each column of the Matrix a single RHS
  Matrix *solve(const Matrix& B) const;

  /// Get the solution tolerance
  /** 
   * 
   * 
   * 
   * @return current solution tolerance
   */
  double tolerance(void) const
  {
    return p_solutionTolerance;
  }

  /// Set the solver tolerance
  /** 
   * 
   * 
   * @param tol 
   */
  void tolerance(const double& tol) 
  {
    p_solutionTolerance = tol;
  }

  /// Get the maximum iterations
  /** 
   * 
   * 
   * 
   * @return 
   */
  int maximumIterations(void) const
  {
    return p_maxIterations;
  }

  /// Set the maximum iterations
  /** 
   * 
   * 
   * @param n current maximum number of iterations
   */
  void maximumIterations(const int& n) 
  {
    p_maxIterations = n;
  }

  /// Allow visits by implemetation visitor
  void accept(ImplementationVisitor& visitor)
  {
    this->p_accept(visitor);
  }

  /// Allow visits by const implemetation visitor
  void accept(ConstImplementationVisitor& visitor) const
  {
    this->p_accept(visitor);
  }

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

  /// Solve w/ the specified RHS, put result in specified vector
  /** 
   * Can be called repeatedly with different @c b and @c x vectors
   * 
   * @param b Vector containing right hand side of linear system
   * @param x Vector containing initial solution estimate; final
   * solution is put into this.
   */
  virtual void p_solve(const Vector& b, Vector& x) const = 0;

  /// Allow visits by implementation visitors
  virtual void p_accept(ImplementationVisitor& visitor) = 0;

  /// Allow visits by implementation visitors
  virtual void p_accept(ConstImplementationVisitor& visitor) const = 0;

};

} // namespace math
} // namespace gridpack



#endif
