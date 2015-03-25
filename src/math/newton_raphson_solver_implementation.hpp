// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   newton_raphson_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2015-03-25 14:52:31 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _newton_raphson_solver_implementation_hpp_
#define _newton_raphson_solver_implementation_hpp_

#include <iostream>
#include <boost/scoped_ptr.hpp>
#include "nonlinear_solver_functions.hpp"
#include "nonlinear_solver_implementation.hpp"
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
                                    FunctionBuilder form_function)
    : NonlinearSolverImplementation(comm, local_size, form_jacobian, form_function),
      p_linear_solver()
  {
    this->configurationKey("NewtonRaphsonSolver");
  }

  
  /// Construct with an existing Jacobian Matrix
  NewtonRaphsonSolverImplementation(Matrix& J,
                                    JacobianBuilder form_jacobian,
                                    FunctionBuilder form_function)
    : NonlinearSolverImplementation(J, form_jacobian, form_function),
      p_linear_solver()
  {
    this->configurationKey("NewtonRaphsonSolver");
  }

  /// Destructor
  /**
   * This must be called simultaneously by all processes involved in
   * the \ref parallel::Communicator "communicator" used for \ref
   * NewtonRaphsonSolverImplementation() "construction".
   */
  ~NewtonRaphsonSolverImplementation(void)
  {
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
  void p_solve(Vector& x)
  {
    NonlinearSolverImplementation::p_solve(x);
    ComplexType stol(1.0e+30);
    ComplexType ftol(1.0e+30);
    int iter(0);

    boost::scoped_ptr<Vector> deltaX(p_X->clone());
    while (real(stol) > p_solutionTolerance && iter < p_maxIterations) {
      p_function(*p_X, *p_F);
      p_F->scale(-1.0);
      p_jacobian(*p_X, *p_J);
      if (!p_linear_solver) {
        p_linear_solver.reset(new LinearSolver(*p_J));
        p_linear_solver->configure(this->p_configCursor);
      } 
      deltaX->zero();
      p_linear_solver->solve(*p_F, *deltaX);
      stol = deltaX->norm2();
      ftol = p_F->norm2();
      p_X->add(*deltaX);
      iter += 1;
      if (this->processor_rank() == 0) {
        std::cout << "Newton-Raphson "
                  << "iteration " << iter << ": "
                  << "solution residual norm = " << real(stol) << ", "
                  << "function norm = " << real(ftol)
                  << std::endl;
      }
    }
  }

};



} // namespace math
} // namespace gridpack
#endif
