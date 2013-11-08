// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2013-11-08 08:49:10 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------

#ifndef _nonlinear_solver_implementation_hpp_
#define _nonlinear_solver_implementation_hpp_

#include <boost/shared_ptr.hpp>
#include <gridpack/math/nonlinear_solver_functions.hpp>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/configuration/configurable.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NonlinearSolverImplementation
// -------------------------------------------------------------
class NonlinearSolverImplementation 
  : public parallel::Distributed,
    public utility::Configurable,
    private utility::Uncopyable
{
public:

  /// Default constructor.
  NonlinearSolverImplementation(const parallel::Communicator& comm,
                                const int& local_size,
                                JacobianBuilder form_jacobian,
                                FunctionBuilder form_function);

  /// Destructor
  virtual ~NonlinearSolverImplementation(void);

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

  /// Solve w/ using the specified initial guess, put solution in same vector
  void solve(Vector& x);

protected:

  /// A place to compute and store the Jacobian
  boost::shared_ptr<Matrix> p_J;

  /// A place to compute and store the RHS
  boost::shared_ptr<Vector> p_F;

  /// A vector containing to current solution estimate
  boost::shared_ptr<Vector> p_X;

  /// A thing to build a Jacobian
  JacobianBuilder p_jacobian;

  /// A thing to build the RHS
  FunctionBuilder p_function;

  /// The solution residual norm tolerance
  double p_solutionTolerance;

  /// The RHS function norm tolerance
  double p_functionTolerance;

  /// The maximum number of iterations to perform
  double p_maxIterations;

  /// Solve w/ using the specified initial guess (specialized)
  virtual void p_solve(void) = 0;

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props);
};


} // namespace math
} // namespace gridpack


#endif
