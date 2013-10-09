// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   newton_raphson_solver.hpp
 * @author William A. Perkins
 * @date   2013-10-09 12:24:41 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _newton_raphson_solver_hpp_
#define _newton_raphson_solver_hpp_

#include "gridpack/math/nonlinear_solver_interface.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NewtonRaphsonSolver
// -------------------------------------------------------------
/// Solve a system of nonlinear equations using the Newton-Raphson method in parallel
/**
 * This (fully functional) class is intended as an example of how to
 * implement a nonlinear solver.  It has no advantage over the general
 * NonlinearSolver class, which relies on the underlying math library.
 * It is, however, a drop-in replacement for NonlinearSolver.
 * 
 * Users of this class must specify functions or functors that build
 * the \ref JacobianBuilder "Jacobian" Matrix and the \ref
 * FunctionBuilder "right hand side" Vector. Typically, it's best to
 * use functor classes or structs, since extra required information
 * can be available to the Matrix/Vector construction.
 *
 * This class is simply an interface to the
 * NewtonRaphsonSolverImplementation class.
 */
class NewtonRaphsonSolver 
  : public NonlinearSolverInterface {
public:

  /// Default constructor.
  /** 
   * @e Collective.
   *
   * A NonlinearSolver must be constructed simultaneously on all
   * processes involved in @c comm.  
   * 
   * @param comm communicator on which the instance is to exist
   * @param local_size number Jacobian rows and Vector entries to be owned by this process
   * @param form_jacobian function to fill the Jacobian Matrix, \f$\left[ \mathbf{J}\left( \mathbf{x} \right) \right]\f$
   * @param form_function function to fill the RHS function Vector, \f$\mathbf{F}\left( \mathbf{x} \right)\f$
   * 
   * @return new NonlinearSolver instance
   */
  NewtonRaphsonSolver(const parallel::Communicator& comm,
                      const int& local_size,
                      JacobianBuilder form_jacobian,
                      FunctionBuilder form_function);
  
  /// Destructor
  /**
   * This must be called simultaneously by all processes involved in
   * the \ref parallel::Communicator "communicator" used for \ref
   * NewtonRaphsonSolver() "construction".
   */
  ~NewtonRaphsonSolver(void);

  /// Get the solution tolerance
  /** 
   * 
   * 
   * 
   * @return current solution tolerance
   */
  double tolerance(void) const;

  /// Set the solver tolerance
  /** 
   * 
   * 
   * @param tol new solution tolerance
   */
  void tolerance(const double& tol);

  /// Get the maximum iterations
  /** 
   * 
   * 
   * 
   * @return current maximum number of solution iterations
   */
  int maximum_iterations(void) const;


  /// Set the maximum solution iterations
  /** 
   * 
   * 
   * @param n new maximum number of iterations
   */
  void maximum_iterations(const int& n);
};


} // namespace math
} // namespace gridpack

#endif
