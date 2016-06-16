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
 * @date   2016-06-16 13:03:07 d3g096
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

  typedef typename DAESolverImplementation<T, I>::VectorType VectorType;
  typedef typename DAESolverImplementation<T, I>::MatrixType MatrixType;
  typedef typename DAESolverImplementation<T, I>::JacobianBuilder JacobianBuilder;
  typedef typename DAESolverImplementation<T, I>::FunctionBuilder FunctionBuilder;
  typedef typename DAESolverImplementation<T, I>::StepFunction StepFunction;


  /// Default constructor.
  PETScDAESolverImplementation(const parallel::Communicator& comm, 
                               const int local_size,
                               JacobianBuilder& jbuilder,
                               FunctionBuilder& fbuilder)
    : DAESolverImplementation<T, I>(comm, local_size, jbuilder, fbuilder),
      PETScConfigurable(this->communicator()),
      p_ts(),
      p_petsc_J(NULL)
  {
    
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

  /// Do what is necessary to build this instance
  void p_build(const std::string& option_prefix)
  {
    PetscErrorCode ierr(0);
    try {
      ierr = TSCreate(this->communicator(), &p_ts); CHKERRXX(ierr);
      ierr = TSSetOptionsPrefix(p_ts, option_prefix.c_str()); CHKERRXX(ierr);

      SNES snes;
      ierr = TSGetSNES(p_ts, &snes); CHKERRXX(ierr);
      ierr = SNESSetOptionsPrefix(snes, option_prefix.c_str()); CHKERRXX(ierr);

      KSP ksp;
      ierr = SNESGetKSP(snes, &ksp); CHKERRXX(ierr);
      ierr = KSPSetOptionsPrefix(ksp, option_prefix.c_str()); CHKERRXX(ierr);
    
      PC pc;
      ierr = KSPGetPC(ksp, &pc); CHKERRXX(ierr);
      ierr = PCSetOptionsPrefix(pc, option_prefix.c_str()); CHKERRXX(ierr);


      p_petsc_J = PETScMatrix(this->p_J);
      ierr = TSSetIFunction(p_ts, NULL, FormIFunction, this); CHKERRXX(ierr);
      ierr = TSSetIJacobian(p_ts, *p_petsc_J, *p_petsc_J, FormIJacobian, this); CHKERRXX(ierr);
      ierr = TSSetApplicationContext(p_ts, this); CHKERRXX(ierr);
      ierr = TSSetPreStep(p_ts, PreTimeStep); CHKERRXX(ierr);
      ierr = TSSetPostStep(p_ts, PostTimeStep); CHKERRXX(ierr);

      ierr = TSSetProblemType(p_ts, TS_NONLINEAR); CHKERRXX(ierr);
      // ierr = TSSetExactFinalTime(p_ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRXX(ierr);

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
      ierr = TSSetInitialTimeStep(p_ts, t0, deltat0); CHKERRXX(ierr); 
      Vec *xvec(PETScVector(x0));
      ierr = TSSetSolution(p_ts, *xvec);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  
  }                                      


  /// Solve the system
  void p_solve(double& maxtime, int& maxsteps)
  {
    PetscErrorCode ierr(0);
    try {
      ierr = TSSetDuration(p_ts, maxsteps, maxtime); CHKERRXX(ierr);
      ierr = TSSetExactFinalTime(p_ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRXX(ierr); 
      ierr = TSSolve(p_ts, PETSC_NULL);
      // std::cout << this->processor_rank() << ": "
      //           << "DAESolver::solve() returned " << ierr 
      //           << std::endl;
      CHKERRXX(ierr);

      TSConvergedReason reason;
      ierr = TSGetConvergedReason(p_ts, &reason); CHKERRXX(ierr);

      PetscInt nstep;
      ierr = TSGetTimeStepNumber(p_ts,&nstep);CHKERRXX(ierr);
      maxsteps = nstep;

      if (reason >= 0) {

        PetscReal tlast;
        ierr = TSGetSolveTime(p_ts,&tlast);CHKERRXX(ierr);
        maxtime = tlast;
      
        std::cout << this->processor_rank() << ": "
                  << "PETSc DAE Solver converged after " << maxsteps << " steps, "
                  << "actual time = " << maxtime
                  << std::endl;

      } else {
        boost::format f("%d: PETSc DAE Solver diverged after %d steps, reason : %d");
        std::string msg = 
          boost::str(f % this->processor_rank() % maxsteps % reason );
        std::cerr << msg << std::cerr;
        throw gridpack::Exception(msg);
      }
    
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    } catch (const gridpack::Exception& e) {
      throw;
    }
  }


#if PETSC_VERSION_LT(3,5,0)

  /// Routine to assemble Jacobian that is sent to PETSc
  static PetscErrorCode FormIJacobian(TS ts, PetscReal t, Vec x, Vec xdot, 
                                      PetscReal a, Mat *jac, Mat *B, 
                                      MatStructure *flag, void *dummy)
  {
    PetscErrorCode ierr(0);

    // Necessary C cast
    PETScDAESolverImplementation *solver =
      (PETScDAESolverImplementation *)dummy;

    // Copy PETSc's current estimate into 

    // Should be the case, but just make sure
    BOOST_ASSERT(*jac == *solver->p_petsc_J);
    BOOST_ASSERT(*B == *solver->p_petsc_J);

    boost::scoped_ptr<VectorType> 
      xtmp(new VectorType(new PETScVectorImplementation<T, I>(x, false))),
      xdottmp(new VectorType(new PETScVectorImplementation<T, I>(xdot, false)));

    // Call the user-specified function (object) to form the Jacobian
    (solver->p_Jbuilder)(t, *xtmp, *xdottmp, a, solver->p_J);

    *flag = SAME_NONZERO_PATTERN;

    return ierr;
  
  }


#else

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
    (solver->p_Jbuilder)(t, *xtmp, *xdottmp, a, solver->p_J);

    return ierr;
  
  }


#endif

  /// Routine to assemble RHS that is sent to PETSc
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
      ierr = TSGetTime(ts, &thetime); CHKERRXX(ierr);
      solver->p_postStepFunc(thetime);
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
      PetscReal thetime;
      ierr = TSGetTime(ts, &thetime); CHKERRXX(ierr);
      solver->p_preStepFunc(thetime);
    }
    return ierr;
  }


};



} // namespace math
} // namespace gridpack

#endif
