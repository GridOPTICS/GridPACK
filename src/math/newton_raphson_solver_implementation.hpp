// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: newton_raphson_solver_implementation.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created September  9, 2013 by William A. Perkins
// Last Change: 2013-09-10 11:36:14 d3g096
// -------------------------------------------------------------


#ifndef _newton_raphson_solver_implementation_hpp_
#define _newton_raphson_solver_implementation_hpp_

#include <boost/scoped_ptr.hpp>
#include "nonlinear_solver_implementation.hpp"
#include "nonlinear_solver_functions.hpp"
#include "linear_solver.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NewtonRaphsonSolverImplementation
// -------------------------------------------------------------
/// Implementation of Newton-Raphson method to solve a system of nonlinear equations in parallel
/**
 * This class implements the Newton-Raphson method to solve a system of
 * nonlinear equations in the form
 * \f[
 * \left[ \mathbf{J}\left( \mathbf{x} \right) \right] \Delta \mathbf{x} ~ = ~ 
 *    -\mathbf{F}\left( \mathbf{x} \right) 
 * \f]
 * where \f$\mathbf{J}\left( \mathbf{x} \right)\f$ is the Jacobian
 * matrix, \f$\mathbf{x}\f$ is the solution Vector, and
 * \f$\mathbf{F}\left( \mathbf{x} \right)\f$ is some Vector function
 * of \f$\mathbf{x}\f$. 
 *
 * Each successive solution estimate is computed as 
 *
 * \f[
 * \mathbf{x}^{k+1} ~ = \mathbf{x}^{k} + \Delta \mathbf{x}^{k}
 * \f]
 * where \f$\Delta \mathbf{x}^{k}\f$ is determined by solving the linear system 
 * \f[
 * \left[ \mathbf{J} \left( \mathbf{x}^{k} \right) \right] \Delta \mathbf{x}^{k} ~ = ~ 
 *    -\mathbf{F}\left( \mathbf{x}^{k} \right) 
 * \f]
 * and \f$ k \f$ is the number of the previous iteration.  
 *
 * The interative process is ended when the L<sup>2</sup> \ref
 * Vector::norm2() "norm" of \f$ \Delta \mathbf{x}^{k} \f$ is less
 * then some specified small tolerance.
 */
class NewtonRaphsonSolverImplementation 
  : public NonlinearSolverImplementation
{
public:

  /// Default constructor.
  /** 
   * @e Collective.
   *
   * A NonlinearSolverImplementation must be constructed
   * simultaneously on all processes involved in @c comm.
   * 
   * @param comm communicator on which the instance is to exist
   * @param local_size number Jacobian rows and Vector entries to be owned by this process
   * @param form_jacobian function to fill the Jacobian Matrix, \f$\left[ \mathbf{J}\left( \mathbf{x} \right) \right]\f$
   * @param form_function function to fill the RHS function Vector, \f$\mathbf{F}\left( \mathbf{x} \right)\f$
   * 
   * @return new NonlinearSolverImplementation instance
   */
  NewtonRaphsonSolverImplementation(const parallel::Communicator& comm,
                                    const int& local_size,
                                    JacobianBuilder form_jacobian,
                                    FunctionBuilder form_function);

  /// Destructor
  /**
   * This must be called simultaneously by all processes involved in
   * the \ref parallel::Communicator "communicator" used for \ref
   * NewtonRaphsonSolverImplementation() "construction".
   */
  ~NewtonRaphsonSolverImplementation(void);

  /// Get the solution tolerance
  /** 
   * 
   * 
   * 
   * @return current solution tolerance
   */
  double tolerance(void) const
  {
    return p_tolerance;
  }

  /// Set the solver tolerance
  /** 
   * 
   * 
   * @param tol 
   */
  void tolerance(const double& tol) 
  {
    p_tolerance = tol;
  }

  /// Get the maximum iterations
  /** 
   * 
   * 
   * 
   * @return 
   */
  int maximum_iterations(void) const
  {
    return p_max_iterations;
  }

  /// Set the maximum iterations
  /** 
   * 
   * 
   * @param n current maximum number of iterations
   */
  void maximum_iterations(const int& n) 
  {
    p_max_iterations = n;
  }

protected:

  /// Solver tolerance goal
  /**
   * 
   * 
   */
  double p_tolerance;

  /// Maximum number of iterations to perform
  int p_max_iterations;

  /// The linear solver
  boost::scoped_ptr<LinearSolver> p_linear_solver;

  /// Solve w/ using the specified initial guess (specialized)
  void p_solve(void);

};



} // namespace math
} // namespace gridpack
#endif
