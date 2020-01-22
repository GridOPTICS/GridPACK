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
 * @file   dae_solver_interface.hpp
 * @author William A. Perkins
 * @date   2019-12-05 07:51:14 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _dae_solver_interface_hpp_
#define _dae_solver_interface_hpp_

#include <gridpack/math/vector.hpp>
#include <gridpack/math/matrix.hpp>
#include <gridpack/utilities/uncopyable.hpp>

#include <gridpack/math/dae_solver_functions.hpp>
#include <gridpack/math/dae_event.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class DAESolverInterface
// -------------------------------------------------------------
template <typename T, typename I = int>
class DAESolverInterface 
{
public:

  typedef typename DAEBuilder<T, I>::VectorType VectorType;
  typedef typename DAEBuilder<T, I>::MatrixType MatrixType;
  typedef typename DAEBuilder<T, I>::Jacobian JacobianBuilder;
  typedef typename DAEBuilder<T, I>::Function FunctionBuilder;
  typedef typename DAEBuilder<T, I>::StepFunction StepFunction;
  typedef typename gridpack::math::DAEEventManagerT<T, I> EventManager;
  typedef typename boost::shared_ptr<EventManager> EventManagerPtr;
  typedef typename EventManager::Event Event;
  typedef typename boost::shared_ptr<Event> EventPtr;

  /// Default constructor.
  DAESolverInterface(void)
  {
  }

  /// Destructor
  virtual ~DAESolverInterface(void)
  {
  }

  /// Initialize the solver
  /** 
   * 
   * 
   * @param t0 start time, corresponding to @c x0
   * @param deltat0 initial time step
   * @param x0 initial solution corresponding to @c t0
   */
  void initialize(const double& t0,
                  const double& deltat0,
                  VectorType& x0)
  {
    this->p_initialize(t0, deltat0, x0);
  }

  /// Solve the system to when @c end_time or @c maxsteps is exceeded
  /** 
   * This solves the system from the initial time specified by
   * ::initialize() or continues the solution from the end time of the
   * previous call to ::solve().  After the solution completes, the
   * solution is placed back in the Vector passed to ::initialize()
   * 
   * @param maxtime time at which solution is to end (in); actual time
   * at which solution ended (out)
   *
   * @param maxsteps maximum number of time steps to perform (in);
   * actual number of time steps performed
   *
   * @param solution system solution at @c maxtime (out)
   */
  void solve(double& maxtime, int& maxsteps)
  {
    this->p_solve(maxtime, maxsteps);
  }

  /// Set a function to call before each time step
  void preStep(StepFunction& f)
  {
    this->p_preStep(f);
  }

  /// Set a function to call before each time step
  void postStep(StepFunction& f)
  {
    this->p_postStep(f);
  }

  /// Has the solver been terminated by an event
  bool terminated(void) const
  {
    return this->p_terminated();
  }

  /// Reset solver if it has been terminated by an event, maybe
  void terminated(const bool& flag)
  {
    this->p_terminated(flag);
  }

protected:

  /// Initialize the system (specialized)
  virtual void p_initialize(const double& t0,
                            const double& deltat0,
                            VectorType& x0) = 0;
                       

  /// Solve the system (specialized)
  virtual void p_solve(double& maxtime, int& maxsteps) = 0;

  /// Set a function to call before each time step (specialized)
  virtual void p_preStep(StepFunction& f) = 0;

  /// Set a function to call after each time step (specialized)
  virtual void p_postStep(StepFunction& f) = 0;

  /// Has the solver been terminated by an event (specialized)
  virtual bool p_terminated(void) const = 0;

  /// Reset solver if it has been terminated by an event, maybe (specialized)
  virtual void p_terminated(const bool& flag) = 0;
};



} // namespace math
} // namespace gridpack

#endif
