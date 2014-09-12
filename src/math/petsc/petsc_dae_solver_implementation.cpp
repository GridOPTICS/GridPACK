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
 * @date   2014-09-12 14:44:40 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/format.hpp>
#include "petsc/petsc_dae_solver_implementation.hpp"
#include "petsc/petsc_matrix_implementation.hpp"
#include "petsc/petsc_vector_implementation.hpp"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_vector_extractor.hpp"
#include "petsc/petsc_matrix_extractor.hpp"
#include "petsc/petsc_vector_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScDAESolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// PETScDAESolverImplementation::FormIJacobian
// -------------------------------------------------------------

#if PETSC_VERSION_LT(3,5,0)

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

#else

PetscErrorCode 
PETScDAESolverImplementation::FormIJacobian(TS ts, PetscReal t, Vec x, Vec xdot, 
                                            PetscReal a, Mat jac, Mat B, 
                                            void *dummy)
{
  PetscErrorCode ierr(0);

  // Necessary C cast
  PETScDAESolverImplementation *solver =
    (PETScDAESolverImplementation *)dummy;

  // Copy PETSc's current estimate into 

  // Should be the case, but just make sure
  BOOST_ASSERT(jac == solver->p_petsc_J);
  BOOST_ASSERT(B == solver->p_petsc_J);

  boost::scoped_ptr<Vector> 
    xtmp(new Vector(new PETScVectorImplementation(x, false))),
    xdottmp(new Vector(new PETScVectorImplementation(xdot, false)));

  // Call the user-specified function (object) to form the Jacobian
  (solver->p_Jbuilder)(t, *xtmp, *xdottmp, a, solver->p_J);

  return ierr;
  
}

#endif

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
    PETScConfigurable(this->communicator()),
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
    // ierr = TSSetExactFinalTime(p_ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRXX(ierr);

    ierr = TSSetFromOptions(p_ts); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
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
  this->build(props);
}

// -------------------------------------------------------------
// PETScDAESolverImplementation::p_initialize
// -------------------------------------------------------------
void
PETScDAESolverImplementation::p_initialize(const double& t0,
                                           const double& deltat0,
                                           Vector& x0)
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

// -------------------------------------------------------------
// PETScDAESolverImplementation::p_solve
// -------------------------------------------------------------
void 
PETScDAESolverImplementation::p_solve(double& maxtime,
                                      int& maxsteps)
{
  PetscErrorCode ierr(0);
  try {
    ierr = TSSetDuration(p_ts, maxsteps, maxtime); CHKERRXX(ierr);
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


} // namespace math
} // namespace gridpack
