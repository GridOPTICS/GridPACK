// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dae_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2023-09-13 08:34:39 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _dae_solver_implementation_hpp_
#define _dae_solver_implementation_hpp_

#include <boost/shared_ptr.hpp>
#include <gridpack/math/dae_solver_interface.hpp>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/configuration/configurable.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class DAESolverImplementation
// -------------------------------------------------------------
template <typename T, typename I = int>
class DAESolverImplementation 
  : public DAESolverInterface<T, I>,
    public parallel::Distributed,
    public utility::Configurable,
    private utility::Uncopyable
{
public:

  typedef typename DAESolverInterface<T, I>::VectorType VectorType;
  typedef typename DAESolverInterface<T, I>::MatrixType MatrixType;
  typedef typename DAESolverInterface<T, I>::JacobianBuilder JacobianBuilder;
  typedef typename DAESolverInterface<T, I>::FunctionBuilder FunctionBuilder;
  typedef typename DAESolverInterface<T, I>::StepFunction StepFunction;
  typedef typename DAESolverInterface<T, I>::EventManagerPtr EventManagerPtr;

  /// Default constructor.
  DAESolverImplementation(const parallel::Communicator& comm, 
                          const int local_size,
                          JacobianBuilder& jbuilder,
                          FunctionBuilder& fbuilder,
                          EventManagerPtr eman)
    : parallel::Distributed(comm),
      utility::Configurable("DAESolver"),
      utility::Uncopyable(),
      p_J(new MatrixType(comm,local_size,local_size)),
      p_J_allocated(true),
      p_Fbuilder(fbuilder), p_Jbuilder(jbuilder),
      p_eventManager(eman),
      p_doAdaptive(true)
  { }

  DAESolverImplementation(const parallel::Communicator& comm, 
                          const int local_size,
			  MatrixType* J,
                          JacobianBuilder& jbuilder,
                          FunctionBuilder& fbuilder,
                          EventManagerPtr eman)
    : parallel::Distributed(comm),
      utility::Configurable("DAESolver"),
      utility::Uncopyable(),
      p_J(J),
      p_J_allocated(false),
      p_Fbuilder(fbuilder), p_Jbuilder(jbuilder),
      p_eventManager(eman),
      p_doAdaptive(true)
  {
  }


  /// Destructor
  ~DAESolverImplementation(void)
  {
    if (p_J_allocated) delete p_J;
  }

protected:

  /// A matrix to hold the Jacobian
  MatrixType *p_J;

  /// Is p_J locally allocated?
  bool p_J_allocated;

  /// A function to build the RHS vector
  FunctionBuilder p_Fbuilder;

  /// A function to build the Jacobian
  JacobianBuilder p_Jbuilder;

  /// An optional function to call before each time step
  StepFunction p_preStepFunc;

  /// An optional function to call after each time step
  StepFunction p_postStepFunc;

  /// An optional event manager
  EventManagerPtr p_eventManager;

  /// Is the time stepper adaptive?
  bool p_doAdaptive;

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props)
  {
    if (props) {
      p_doAdaptive = props->get("Adaptive", p_doAdaptive);
    }
  }

  /// Set a function to call before each time step (specialized)
  void p_preStep(StepFunction& f)
  {
    p_preStepFunc = f;
  }

  /// Set a function to call after each time step (specialized)
  void p_postStep(StepFunction& f)
  {
    p_postStepFunc = f;
  }

  /// Has the solver been terminated by an event (specialized)
  /** 
   * The event manager records terminating events. Library-specific
   * implementations should probably query the underlying solver
   * library.
   * 
   * 
   * @return true if an event has terminated further solution
   */
  bool p_terminated(void) const
  {
    bool result(false);
    if (p_eventManager) {
      result = (this->communicator()).any(p_eventManager->terminated());
    }
    
    return result;
  }

  /// Reset solver if it has been terminated by an event, maybe (specialized)
  void p_terminated(const bool& flag)
  {
    if (p_eventManager) {
      p_eventManager->terminated(flag);
    }
  }
};



} // namespace math
} // namespace gridpack

#endif
