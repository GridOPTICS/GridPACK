// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dae_solver.hpp
 * @author William A. Perkins
 * @date   2013-11-14 13:57:55 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _dae_solver_hpp_
#define _dae_solver_hpp_

#include <gridpack/math/dae_solver_implementation.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class DAESolver
// -------------------------------------------------------------
/// A solver for systems differential algebraic equations (DAE)
/**
 * 
 * 
 */
class DAESolver 
  : public parallel::WrappedDistributed,
    public utility::WrappedConfigurable,
    private utility::Uncopyable
{
public:

  /// Default constructor.
  /** 
   * 
   * 
   * @param comm parallel environment for this instance
   * @param local_size size of problem owned by this process (rows in Jacobian, elements in solution Vector)
   * @param jbuilder function/functor to build Jacobian
   * @param fbuilder function/functor to build function
   * 
   * @return 
   */
  DAESolver(const parallel::Communicator& comm, 
            const int local_size,
            DAEJacobianBuilder& jbuilder,
            DAEFunctionBuilder& fbuilder);

  /// Destructor
  ~DAESolver(void);

  /// Initialize the solver
  /** 
   * 
   * 
   * @param t0 start time
   * @param deltat0 initial time step
   * @param x0 Vector to use for the solution, initially filled with
   * values corresponding to @c t0
   */
  void initialize(const double& t0,
                  const double& deltat0,
                  Vector& x0)
  {
    p_impl->initialize(t0, deltat0, x0);
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
    p_impl->solve(maxtime, maxsteps);
  }

  
protected:
  
  /// Where things really happen
  boost::scoped_ptr<DAESolverImplementation> p_impl;

  /// Set the implementation
  void p_setImpl(DAESolverImplementation *impl);
};



} // namespace math
} // namespace gridpack

#endif
