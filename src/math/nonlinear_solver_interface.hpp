// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver_interface.hpp
 * @author William A. Perkins
 * @date   2015-03-25 14:28:06 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _nonlinear_solver_interface_hpp_
#define _nonlinear_solver_interface_hpp_

#include <gridpack/math/vector.hpp>
#include <gridpack/math/matrix.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NonlinearSolverInterface
// -------------------------------------------------------------
/// Interface to a solve system of nonlinear equations in parallel
/**
 * This serves as a base for classes that solve a system of nonlinear
 * equations. While not strictly abstract, it has no function if
 * instantiated on its own.
 * 
 * It encapuslates the nonlinear system solver of some
 * underlying implementation. The Pimpl idiom is used for \ref
 * NonlinearSolverImplementation "implementation", so user code is
 * completely independent of the underlying library. This class simply
 * provides an interface to a specific \ref
 * NonlinearSolverImplementation "implementation".  Subclasses are
 * required to call p_set_impl() at construction to set the \ref
 * NonlinearSolverImplementation "implementation".
 */
template <typename T, typename I = int>
class NonlinearSolverInterface
{
public:

  typedef VectorT<T, I> VectorType;
  typedef MatrixT<T, I> MatrixType;

  /// Default constructor.
  NonlinearSolverInterface()
  {
  }

  /// Destructor
  virtual ~NonlinearSolverInterface(void)
  {
  }

  /// Get the solution tolerance
  /** 
   * 
   * 
   * 
   * @return current solution tolerance
   */
  double tolerance(void) const
  {
    return p_tolerance();
  }

  /// Set the solver tolerance
  /** 
   * 
   * 
   * @param tol new solution tolerance
   */
  void tolerance(const double& tol)
  {
    p_tolerance(tol);
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
    return p_maximumIterations();
  }


  /// Set the maximum solution iterations
  /** 
   * 
   * 
   * @param n new maximum number of iterations
   */
  void maximumIterations(const int& n)
  {
    p_maximumIterations(n);
  }

  /// Solve w/ the specified initial estimated, put result in same vector
  /** 
   * This solves the system of nonlinear equations using the contents
   * of @c x as an initial solution estimate.  The final result is
   * placed back in @c x upon completion.
   * 
   * @param x solution Vector
   */
  void solve(VectorType& x)
  {
    p_solve(x);
  }

protected:

  /// Get the solution tolerance (specialized)
  virtual double p_tolerance(void) const = 0;

  /// Set the solver tolerance (specialized)
  virtual void p_tolerance(const double& tol) = 0;

  /// Get the maximum iterations (specialized)
  virtual int p_maximumIterations(void) const = 0;

  /// Set the maximum solution iterations  (specialized)
  virtual void p_maximumIterations(const int& n) = 0;

  /// Solve w/ the specified initial estimated, put result in same vector
  virtual void p_solve(VectorType& x) = 0;
  
};


} // namespace math
} // namespace gridpack


#endif
