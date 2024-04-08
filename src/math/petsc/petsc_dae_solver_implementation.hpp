// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

// -------------------------------------------------------------
/**
 * @file   petsc_dae_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2023-09-13 07:43:15 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_dae_solver_implementation_hpp_
#define _petsc_dae_solver_implementation_hpp_

#include <boost/format.hpp>
#include <petscts.h>

#include "petsc_exception.hpp"
#include "dae_solver_implementation.hpp"
#include "petsc_vector_implementation.hpp"
#include "petsc_vector_extractor.hpp"
#include "petsc_matrix_extractor.hpp"
#include "petsc_configurable.hpp"
#include "complex_operators.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScDAESolverImplementation
// -------------------------------------------------------------
template <typename T, typename I = int>
class PETScDAESolverImplementation 
  : public DAESolverImplementation<T, I>,
    private PETScConfigurable
{
public:

  typedef typename DAESolverInterface<T, I>::VectorType VectorType;
  typedef typename DAESolverInterface<T, I>::MatrixType MatrixType;
  typedef typename DAESolverInterface<T, I>::JacobianBuilder JacobianBuilder;
  typedef typename DAESolverInterface<T, I>::FunctionBuilder FunctionBuilder;
  typedef typename DAESolverInterface<T, I>::RHSFunctionBuilder RHSFunctionBuilder;
  typedef typename DAESolverInterface<T, I>::PreStepFunction PreStepFunction;
  typedef typename DAESolverInterface<T, I>::PostStepFunction PostStepFunction;
  typedef typename DAESolverInterface<T, I>::EventManagerPtr EventManagerPtr;

  /// Default constructor.
  PETScDAESolverImplementation(const parallel::Communicator& comm, 
                               const int local_size,
                               JacobianBuilder& jbuilder,
                               FunctionBuilder& fbuilder,
                               EventManagerPtr eman)
    : DAESolverImplementation<T, I>(comm, local_size, jbuilder, fbuilder, eman),
      PETScConfigurable(this->communicator()),
      p_ts(),
      p_petsc_J(NULL),
      p_eventv(),
      p_termFlag(false)
  {
    if (eman) p_eventv.resize(eman->size());
  }

  PETScDAESolverImplementation(const parallel::Communicator& comm, 
                               const int local_size,
			       MatrixType* J,
                               JacobianBuilder& jbuilder,
                               FunctionBuilder& fbuilder,
                               EventManagerPtr eman)
    : DAESolverImplementation<T, I>(comm, local_size, J, jbuilder, fbuilder, eman),
      PETScConfigurable(this->communicator()),
      p_ts(),
      p_petsc_J(NULL),
      p_eventv(),
      p_termFlag(false)
  {
    if (eman) p_eventv.resize(eman->size());
  }

  PETScDAESolverImplementation(const parallel::Communicator& comm, 
                               const int local_size,
			       MatrixType* J,
                               JacobianBuilder& jbuilder,
                               FunctionBuilder& fbuilder,
			       RHSFunctionBuilder& rbuilder,
                               EventManagerPtr eman)
    : DAESolverImplementation<T, I>(comm, local_size, J, jbuilder, fbuilder, rbuilder,eman),
      PETScConfigurable(this->communicator()),
      p_ts(),
      p_petsc_J(NULL),
      p_eventv(),
      p_termFlag(false)
  {
    if (eman) p_eventv.resize(eman->size());
  }

  /// Destructor
  ~PETScDAESolverImplementation(void)
  {
    PetscErrorCode ierr(0);
    try  {
      PetscBool ok;
      ierr = PetscInitialized(&ok); CHKERRXX(ierr);
      if (ok) {
        ierr = TSDestroy(&p_ts); CHKERRXX(ierr);
      }
    } catch (...) {
      // just eat it
    }
  }


protected:

  /// The actual solver
  TS p_ts;

  /// The Jacobian matrix
  Mat *p_petsc_J;

  /// An array to store event values
  std::vector<PetscScalar> p_eventv;

  /// Has an event terminated integration?
  bool p_termFlag;

  /// Has the solver been terminated by an event (specialized)
  bool p_terminated(void) const
  {
    return this->p_termFlag;
  }

  /// Reset solver if it has been terminated by an event, maybe (specialized)
  void p_terminated(const bool& flag)
  {
    p_termFlag = flag;
    DAESolverImplementation<T, I>::p_terminated(flag);
  }
  
  /// Do what is necessary to build this instance
  void p_build(const std::string& option_prefix)
  {
    PetscErrorCode ierr(0);
    try {
      ierr = TSCreate(this->communicator(), &p_ts); CHKERRXX(ierr);

      ierr = TSSetOptionsPrefix(p_ts, option_prefix.c_str()); CHKERRXX(ierr);

      p_petsc_J = PETScMatrix(*(this->p_J));
      ierr = TSSetIFunction(p_ts, NULL, FormIFunction, this); CHKERRXX(ierr);
      ierr = TSSetIJacobian(p_ts, *p_petsc_J, *p_petsc_J, FormIJacobian, this); CHKERRXX(ierr);
      ierr = TSSetApplicationContext(p_ts, this); CHKERRXX(ierr);
      ierr = TSSetPreStep(p_ts, PreTimeStep); CHKERRXX(ierr);
      ierr = TSSetPostStep(p_ts, PostTimeStep); CHKERRXX(ierr);

      if (this->p_eventManager) {
        int ne(this->p_eventManager->size());
        boost::scoped_array<PetscInt> pdir(new PetscInt[ne]);
        boost::scoped_array<PetscBool> pterm(new PetscBool[ne]);
        const DAEEventDirection *gdir(this->p_eventManager->directions());
        const bool *gterm(this->p_eventManager->terminateFlags());
        
        for (int i = 0; i < ne; ++i) {
          pdir[i] = gdir[i];
          pterm[i] = (gterm[i] ? PETSC_TRUE : PETSC_FALSE);
        }

        ierr = TSSetEventHandler(p_ts, ne, &(pdir[0]), &(pterm[0]),
                                 &EventHandler, &PostEventHandler,
                                 NULL);
      }

      ierr = TSSetProblemType(p_ts, TS_NONLINEAR); CHKERRXX(ierr);
      // ierr = TSSetExactFinalTime(p_ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRXX(ierr);

      TSAdapt adapt;

      ierr = TSGetAdapt(p_ts, &adapt); CHKERRXX(ierr);
      if (this->p_doAdaptive) {
        ierr = TSAdaptSetType(adapt, TSADAPTBASIC); CHKERRXX(ierr);
      } else {
        ierr = TSAdaptSetType(adapt, TSADAPTNONE); CHKERRXX(ierr);
      }

      ierr = TSSetFromOptions(p_ts); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
    
  }

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props)
  {
    DAESolverImplementation<T, I>::p_configure(props);
    this->build(props);
  }


  /// Initialize the system (specialized)
  void p_initialize(const double& t0,
                    const double& deltat0,
                    VectorType& x0)
  {
    PetscErrorCode ierr(0);
    try {
      ierr = TSSetTime(p_ts, t0); CHKERRXX(ierr);
      ierr = TSSetTimeStep(p_ts, deltat0); CHKERRXX(ierr);
      ierr = TSSetPostEventIntervalStep(p_ts, deltat0); CHKERRXX(ierr);
      Vec *xvec(PETScVector(x0));
      ierr = TSSetSolution(p_ts, *xvec);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  
  }

  /// Restart step
  void p_restartstep()
  {
    PetscErrorCode ierr(0);
    try {
      ierr = TSRestartStep(p_ts);CHKERRXX(ierr);
    } catch(const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Reuse preconditioner
  void p_reusepreconditioner(int niter)
  {
    PetscErrorCode ierr(0);
    try {
      SNES snes;
      ierr = TSGetSNES(p_ts,&snes);CHKERRXX(ierr);
      ierr = SNESSetLagPreconditioner(snes,niter);CHKERRXX(ierr);
    } catch(const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Get time step
  double p_gettimestep()
  {
    PetscErrorCode ierr(0);
    try {
      double dt;
      ierr = TSGetTimeStep(p_ts,&dt);CHKERRXX(ierr);
      return dt;
    } catch(const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }
  
  /// Get number of steps
  int p_getstepnumber()
  {
    PetscErrorCode ierr(0);
    try {
      int nsteps;
      ierr = TSGetStepNumber(p_ts,&nsteps);CHKERRXX(ierr);
      return nsteps;
    } catch(const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }

  /// Solve the system
  void p_solve(double& maxtime, int& maxsteps)
  {
    PetscErrorCode ierr(0);
    try {
      ierr = TSSetMaxSteps(p_ts, maxsteps); CHKERRXX(ierr);
      ierr = TSSetMaxTime(p_ts, maxtime); CHKERRXX(ierr);
      ierr = TSSetExactFinalTime(p_ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRXX(ierr); 
      ierr = TSSolve(p_ts, PETSC_NULL);
      // std::cout << this->processor_rank() << ": "
      //           << "DAESolver::solve() returned " << ierr 
      //           << std::endl;
      CHKERRXX(ierr);

      TSConvergedReason reason;
      ierr = TSGetConvergedReason(p_ts, &reason); CHKERRXX(ierr);

      PetscInt nstep;
      ierr = TSGetStepNumber(p_ts,&nstep);CHKERRXX(ierr);
      maxsteps = nstep;

      if (reason >= 0) {

        PetscReal tlast;
        ierr = TSGetSolveTime(p_ts,&tlast);CHKERRXX(ierr);
        maxtime = tlast;

        if (reason == TS_CONVERGED_EVENT) {

          // Note: the PETSc reason is already reduced across all
          // ranks. There's no need to do it a gain.
          
          this->p_termFlag = true;

          std::cout << this->processor_rank() << ": "
                    << "A DAE solver termination event occurred after "
                    << maxsteps << " steps, "
                    << "actual time = " << maxtime
                    << std::endl;
        } else {
          std::cout << this->processor_rank() << ": "
                    << "PETSc DAE Solver converged after " << maxsteps << " steps, "
                    << "actual time = " << maxtime
                    << std::endl;
        }
      } else {
        boost::format f("%d: PETSc DAE Solver diverged after %d steps, reason : %d");
        std::string msg = 
          boost::str(f % this->processor_rank() % maxsteps % reason );
        std::cerr << msg << std::endl;
        throw gridpack::Exception(msg);
      }
    
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    } catch (const gridpack::Exception& e) {
      throw;
    }
  }

  /// Routine to assemble Jacobian that is sent to PETSc
  static PetscErrorCode FormIJacobian(TS ts, PetscReal t, Vec x, Vec xdot, 
                                      PetscReal a, Mat jac, Mat B, 
                                      void *dummy)
  {
    PetscErrorCode ierr(0);

    // Necessary C cast
    PETScDAESolverImplementation *solver =
      (PETScDAESolverImplementation *)dummy;

    // Copy PETSc's current estimate into 

    // Should be the case, but just make sure
    BOOST_ASSERT(jac == *(solver->p_petsc_J));
    BOOST_ASSERT(B == *(solver->p_petsc_J));

    boost::scoped_ptr<VectorType> 
      xtmp(new VectorType(new PETScVectorImplementation<T, I>(x, false))),
      xdottmp(new VectorType(new PETScVectorImplementation<T, I>(xdot, false)));

    // Call the user-specified function (object) to form the Jacobian
    (solver->p_Jbuilder)(t, *xtmp, *xdottmp, a, *(solver->p_J));

    return ierr;
  
  }


  /// Routine to assemble IFunction that is sent to PETSc
  static PetscErrorCode FormIFunction(TS ts, PetscReal t, Vec x, Vec xdot, 
                                      Vec F, void *dummy)
  {
    PetscErrorCode ierr(0);
    // Necessary C cast
    PETScDAESolverImplementation *solver =
      (PETScDAESolverImplementation *)dummy;

    boost::scoped_ptr<VectorType> 
      xtmp(new VectorType(new PETScVectorImplementation<T, I>(x, false))),
      xdottmp(new VectorType(new PETScVectorImplementation<T, I>(xdot, false))),
      ftmp(new VectorType(new PETScVectorImplementation<T, I>(F, false)));

    (solver->p_Fbuilder)(t, *xtmp, *xdottmp, *ftmp);
    return ierr;
  }

  /// Routine to assemble RHSFunction that is sent to PETSc
  static PetscErrorCode FormRHSFunction(TS ts, PetscReal t, Vec x, 
                                      Vec F, void *dummy)
  {
    PetscErrorCode ierr(0);
    // Necessary C cast
    PETScDAESolverImplementation *solver =
      (PETScDAESolverImplementation *)dummy;

    boost::scoped_ptr<VectorType> 
      xtmp(new VectorType(new PETScVectorImplementation<T, I>(x, false))),
      ftmp(new VectorType(new PETScVectorImplementation<T, I>(F, false)));

    (solver->p_RHSbuilder)(t, *xtmp, *ftmp);
    return ierr;
  }

  /// Routine called after each time step
  static PetscErrorCode PostTimeStep(TS ts)
  {
    PetscErrorCode ierr(0);

    void *dummy;
    ierr = TSGetApplicationContext(ts, &dummy); CHKERRXX(ierr);

    // Necessary C cast
    PETScDAESolverImplementation *solver =
      (PETScDAESolverImplementation *)dummy;

    if (solver->p_postStepFunc) {
      PetscReal thetime;
      Vec       x;
      ierr = TSGetTime(ts, &thetime); CHKERRXX(ierr);
      ierr = TSGetSolution(ts,&x);CHKERRXX(ierr);

      boost::scoped_ptr<VectorType> 
	xtmp(new VectorType(new PETScVectorImplementation<T, I>(x, false)));
      solver->p_postStepFunc(thetime, *xtmp);
    }
    return ierr;
  }

  /// Routine called before each time step
  static PetscErrorCode PreTimeStep(TS ts)
  {
    PetscErrorCode ierr(0);

    void *dummy;
    ierr = TSGetApplicationContext(ts, &dummy); CHKERRXX(ierr);

    // Necessary C cast
    PETScDAESolverImplementation *solver =
      (PETScDAESolverImplementation *)dummy;

    if (solver->p_preStepFunc) {
      PetscReal thetime, thetimestep;
      Vec       x;
      ierr = TSGetTime(ts, &thetime); CHKERRXX(ierr);
      ierr = TSGetTimeStep(ts, &thetimestep); CHKERRXX(ierr);
      ierr = TSGetSolution(ts,&x);CHKERRXX(ierr);

      boost::scoped_ptr<VectorType> 
	xtmp(new VectorType(new PETScVectorImplementation<T, I>(x, false)));
      solver->p_preStepFunc(thetime, thetimestep, *xtmp);
    }
    return ierr;
  }

  /// Routine called every step if an event manager is specified
  static PetscErrorCode
  EventHandler(TS ts, PetscReal t, Vec U, PetscScalar fvalue[], void* ctx)
  {
    PetscErrorCode ierr(0);
    void *dummy;
    ierr = TSGetApplicationContext(ts, &dummy); CHKERRXX(ierr);

    // Necessary C cast
    PETScDAESolverImplementation *solver =
      (PETScDAESolverImplementation *)dummy;

    

    boost::scoped_ptr<VectorType>
      state(new VectorType(new PETScVectorImplementation<T, I>(U, false)));

    // This gets a little tricky.  If PETSc is built w/ complex,
    // fvalue will be complex, and evalues, whether real or complex,
    // can be assigned directly.  If T is complex and PetscScalar is
    // real, then equate<> will return the real part of T.  This is
    // what PETSc does internally to check a complex event value. So,
    // this *should* work for PETSc.

    const T *evalues = solver->p_eventManager->values(t, *state);
    for (int i = 0; i < solver->p_eventManager->size(); ++i) {
      fvalue[i] = equate<PetscScalar, T>(evalues[i]);
    }

    ierr = VecAssemblyBegin(U);
    ierr = VecAssemblyEnd(U);
    return ierr;
  }

  /// Routine called if events are triggered
  static PetscErrorCode
  PostEventHandler(TS ts, PetscInt nevents_zero, PetscInt events_zero[],
                   PetscReal t, Vec U, PetscBool forwardsolve, void* ctx)
  {
    PetscErrorCode ierr(0);
    void *dummy;
    ierr = TSGetApplicationContext(ts, &dummy); CHKERRXX(ierr);

    // Necessary C cast
    PETScDAESolverImplementation *solver =
      (PETScDAESolverImplementation *)dummy;

    boost::scoped_ptr<VectorType>
      state(new VectorType(new PETScVectorImplementation<T, I>(U, false)));

    solver->p_eventManager->handle(nevents_zero, events_zero, t, *state);
    return ierr;
  }
};



} // namespace math
} // namespace gridpack

#endif
