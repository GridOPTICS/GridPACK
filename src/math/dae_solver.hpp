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
 * @date   2015-05-07 13:15:37 d3g096
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
//  class DAESolverT
// -------------------------------------------------------------
/// A solver for systems differential algebraic equations (DAE)
/**
 * 
 * 
 */
template <typename T, typename I = int>
class DAESolverT 
  : public DAESolverInterface<T, I>,
    public parallel::WrappedDistributed,
    public utility::WrappedConfigurable,
    private utility::Uncopyable
{
public:

  typedef typename DAESolverInterface<T, I>::VectorType VectorType;
  typedef typename DAESolverInterface<T, I>::MatrixType MatrixType;
  typedef typename DAESolverInterface<T, I>::JacobianBuilder JacobianBuilder;
  typedef typename DAESolverInterface<T, I>::FunctionBuilder FunctionBuilder;
  typedef typename DAESolverInterface<T, I>::StepFunction StepFunction;


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
  DAESolverT(const parallel::Communicator& comm, 
             const int local_size,
             JacobianBuilder& jbuilder,
             FunctionBuilder& fbuilder);

  /// Destructor
  ~DAESolverT(void)
  {}

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
                  VectorType& x0)
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
  boost::scoped_ptr<DAESolverImplementation<T, I> > p_impl;

  /// Set the implementation
  void p_setImpl(DAESolverImplementation<T, I> *impl)
  {
    p_impl.reset(impl);
    p_setDistributed(p_impl.get());
    p_setConfigurable(p_impl.get());
  }


  /// Initialize the system (specialized)
  void p_initialize(const double& t0,
                    const double& deltat0,
                    VectorType& x0)
  {
    p_impl->initialize(t0, deltat0, x0);
  }
                       

  /// Solve the system (specialized)
  void p_solve(double& maxtime, int& maxsteps)
  {
    p_impl->solve(maxtime, maxsteps);
  }

  /// Set a function to call before each time step (specialized)
  void p_preStep(StepFunction& f)
  {
    p_impl->preStep(f);
  }

  /// Set a function to call after each time step (specialized)
  void p_postStep(StepFunction& f)
  {
    p_impl->postStep(f);
  }
};

typedef DAESolverT<ComplexType> ComplexDAESolver;
typedef DAESolverT<RealType> RealDAESolver;
typedef ComplexDAESolver DAESolver;

} // namespace math
} // namespace gridpack

#endif
