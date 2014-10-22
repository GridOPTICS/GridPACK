// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------

// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   linear_solver_interface.hpp
 * @author William A. Perkins
 * @date   2014-10-22 09:12:07 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October 21, 2014 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _linear_solver_interface_hpp_
#define _linear_solver_interface_hpp_

#include "gridpack/math/matrix.hpp"
#include "gridpack/math/implementation_visitable.hpp"

namespace gridpack {
namespace math {


// -------------------------------------------------------------
//  class BaseLinearSolverInterface
// -------------------------------------------------------------
class BaseLinearSolverInterface 
  : public ImplementationVisitable
{
public:

  /// Default constructor.
  BaseLinearSolverInterface(void)
    : ImplementationVisitable()
  {}

  /// Destructor
  ~BaseLinearSolverInterface(void)
  {}

  /// Get the solution tolerance
  /** 
   * 
   * 
   * 
   * @return current solution tolerance
   */
  double tolerance(void) const
  {
    return this->p_tolerance();
  }

  /// Set the solver tolerance
  /** 
   * 
   * 
   * @param tol new solution tolerance
   */
  void tolerance(const double& tol)
  {
    this->p_tolerance(tol);
  }

  /// Get the maximum iterations
  /** 
   * 
   * 
   * 
   * @return current maximum number of solution iterations
   */
  int maximumIterations(void) const
  {
    return this->p_maximumIterations();
  }


  /// Set the maximum solution iterations
  /** 
   * 
   * 
   * @param n new maximum number of iterations
   */
  void maximumIterations(const int& n)
  {
    this->p_maximumIterations(n);
  }

  /// Solve w/ the specified RHS, put result in specified vector
  /** 
   * @e Collective.
   *
   * Solve a linear system of equations. When called, Vector @c x
   * should contain the initial solution estimate.  The final solution
   * is returned in @c x.  
   *
   * The \ref parallel::Communicator "communicator" @c x and @c b must
   * be the same and match that of the coefficient Matrix used for
   * construction. The length of both @c x and @c b must the number of
   * columns in the coeffienct Matrix used for constructor or passed
   * the last call to setMatrix().  If these conditions are not met,
   * an \ref Exception "exception" is thrown.
   *  
   * @param b Vector containing right hand side of linear system
   * @param x solution Vector
   */
  void solve(const Vector& b, Vector& x) const
  {
    this->p_solve(b, x);
  }
  /// Solve multiple systems w/ each column of the Matrix a single RHS
  /** 
   * 
   * 
   * @param B RHS matrix -- each column is used as a RHS Vector
   * 
   * @return @e dense solution Matrix -- each column is the solution for the corresponding column in @c B
   */
  Matrix *solve(const Matrix& B) const
  {
    return this->p_solve(B);
  }


protected:

  /// Get the solution tolerance (specialized)
  virtual double p_tolerance(void) const = 0;

  /// Set the solver tolerance (specialized)
  virtual void p_tolerance(const double& tol) = 0;

  /// Get the maximum iterations (specialized)
  virtual int p_maximumIterations(void) const = 0;

  /// Set the maximum solution iterations
  virtual void p_maximumIterations(const int& n) = 0;

  /// Solve w/ the specified RHS, put result in specified vector
  /** 
   * Can be called repeatedly with different @c b and @c x vectors
   * 
   * @param b Vector containing right hand side of linear system
   * @param x Vector containing initial solution estimate; final
   * solution is put into this.
   */
  virtual void p_solve(const Vector& b, Vector& x) const = 0;

  /// Solve multiple systems w/ each column of the Matrix a single RHS
  virtual Matrix *p_solve(const Matrix& B) const = 0;

};



} // namespace math
} // namespace gridpack


#endif
