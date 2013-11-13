// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_dae_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2013-11-13 13:58:44 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "petsc/petsc_dae_solver_implementation.hpp"
#include "petsc/petsc_matrix_implementation.hpp"
#include "petsc/petsc_vector_implementation.hpp"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_vector_extractor.hpp"
#include "petsc/petsc_matrix_extractor.hpp"
#include "petsc/petsc_vector_implementation.hpp"
#include "petsc/petsc_configuration.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScDAESolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// PETScDAESolverImplementation::FormIJacobian
// -------------------------------------------------------------
PetscErrorCode 
PETScDAESolverImplementation::FormIJacobian(TS ts, PetscReal t, Vec x, Vec xdot, 
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

  boost::scoped_ptr<Vector> 
    xtmp(new Vector(new PETScVectorImplementation(x, false))),
    xdottmp(new Vector(new PETScVectorImplementation(xdot, false)));

  // Call the user-specified function (object) to form the Jacobian
  (solver->p_Jbuilder)(t, *xtmp, *xdottmp, a, solver->p_J);

  *flag = SAME_NONZERO_PATTERN;

  return ierr;
  
}

// -------------------------------------------------------------
// PETScDAESolverImplementation::FormIFunction
// -------------------------------------------------------------
PetscErrorCode 
PETScDAESolverImplementation::FormIFunction(TS ts, PetscReal t, Vec x, Vec xdot, 
                                            Vec F, void *dummy)
{
  PetscErrorCode ierr(0);
  // Necessary C cast
  PETScDAESolverImplementation *solver =
    (PETScDAESolverImplementation *)dummy;

  boost::scoped_ptr<Vector> 
    xtmp(new Vector(new PETScVectorImplementation(x, false))),
    xdottmp(new Vector(new PETScVectorImplementation(xdot, false))),
    ftmp(new Vector(new PETScVectorImplementation(F, false)));

  (solver->p_Fbuilder)(t, *xtmp, *xdottmp, *ftmp);
  return ierr;
}


// -------------------------------------------------------------
// PETScDAESolverImplementation:: constructors / destructor
// -------------------------------------------------------------
PETScDAESolverImplementation::PETScDAESolverImplementation(const parallel::Communicator& comm, 
                                                           const int local_size,
                                                           DAEJacobianBuilder& jbuilder,
                                                           DAEFunctionBuilder& fbuilder)
  : DAESolverImplementation(comm, local_size, jbuilder, fbuilder),
    p_ts(),
    p_petsc_J(NULL)
{
  
}

PETScDAESolverImplementation::~PETScDAESolverImplementation(void)
{
  PetscErrorCode ierr;
  try  {
    PetscBool ok;
    ierr = PetscInitialized(&ok);
    if (ok) {
      ierr = TSDestroy(&p_ts);
    }
  } catch (...) {
    // just eat it
  }
}

// -------------------------------------------------------------
// PETScDAESolverImplementation::p_build
// -------------------------------------------------------------
void 
PETScDAESolverImplementation::p_build(const std::string& option_prefix)
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


    p_petsc_J = PETScMatrix(p_J);
    ierr = TSSetIFunction(p_ts, NULL, FormIFunction, this); CHKERRXX(ierr);
    ierr = TSSetIJacobian(p_ts, *p_petsc_J, *p_petsc_J, FormIJacobian, this); CHKERRXX(ierr);

    ierr = TSSetProblemType(p_ts, TS_NONLINEAR); CHKERRXX(ierr);
    ierr = TSSetFromOptions(p_ts); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
    
}

// -------------------------------------------------------------
// PETScDAESolverImplementation::p_configure
// -------------------------------------------------------------
void 
PETScDAESolverImplementation::p_configure(utility::Configuration::CursorPtr props)
{
  DAESolverImplementation::p_configure(props);
  std::string prefix(petscProcessOptions(this->communicator(), props));
  p_build(prefix);
}

// -------------------------------------------------------------
// PETScDAESolverImplementation::p_solve
// -------------------------------------------------------------
void PETScDAESolverImplementation::p_solve(const double& time,
                                           const double& deltat0,
                                           double& maxtime,
                                           int& maxsteps,
                                           Vector& solution)
{
  PetscErrorCode ierr(0);
  try {
    ierr = TSSetInitialTimeStep(p_ts, time, deltat0); CHKERRXX(ierr); 
    ierr = TSSetDuration(p_ts, maxsteps, maxtime); CHKERRXX(ierr);

    Vec *x(PETScVector(solution));

    ierr = TSSolve(p_ts,*x);
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

      std::cout << this->processor_rank() << ": " 
                << "PETSc DAE Solver diverged after " << maxsteps << " steps, " 
                << "reason: " << reason 
                << std::endl;
    }

  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  
}


} // namespace math
} // namespace gridpack
