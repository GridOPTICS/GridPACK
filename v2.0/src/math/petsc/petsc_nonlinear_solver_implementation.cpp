// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_nonlinear_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2014-09-18 09:27:16 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/format.hpp>

#include "petscsys.h"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_nonlinear_solver_implementation.hpp"
#include "petsc/petsc_matrix_implementation.hpp"
#include "petsc/petsc_matrix_extractor.hpp"
#include "petsc/petsc_vector_extractor.hpp"
#include "petsc/petsc_vector_implementation.hpp"

static PetscErrorCode  
MonitorNorms(SNES snes, PetscInt its, PetscReal fgnorm, void *dummy)
{
  PetscViewer    viewer = PETSC_VIEWER_STDOUT_WORLD;

  Vec dx;
  SNESGetSolutionUpdate(snes, &dx);
  PetscReal dxnorm;
  VecNorm(dx, NORM_2, &dxnorm);
  PetscInt tablevel;
  PetscObjectGetTabLevel(((PetscObject)snes), &tablevel);
  PetscViewerASCIIAddTab(viewer, tablevel);
  PetscViewerASCIIPrintf(viewer,"%3D SNES Function norm %14.12e, Solution residual norm %14.12e\n",
                         its, (double)fgnorm, (double)dxnorm);
  PetscViewerASCIISubtractTab(viewer, tablevel);
  return(0);
}

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PetscNonlinearSolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// PetscNonlinearSolverImplementation:: constructors / destructor
// -------------------------------------------------------------
PetscNonlinearSolverImplementation::PetscNonlinearSolverImplementation(const parallel::Communicator& comm,
                                                                       const int& local_size,
                                                                       JacobianBuilder form_jacobian,
                                                                       FunctionBuilder form_function)
  : NonlinearSolverImplementation(comm, local_size, form_jacobian, form_function),
    PETScConfigurable(this->communicator()),
    p_snes(), 
    p_petsc_J(), p_petsc_F(),
    p_petsc_X()                 // set by p_solve()
{
  
}

PetscNonlinearSolverImplementation::PetscNonlinearSolverImplementation(Matrix& J,
                                                                       JacobianBuilder form_jacobian,
                                                                       FunctionBuilder form_function)
  : NonlinearSolverImplementation(J, form_jacobian, form_function),
    PETScConfigurable(this->communicator()),
    p_snes(), 
    p_petsc_J(), p_petsc_F(),
    p_petsc_X()                 // set by p_solve()
{
  
}

PetscNonlinearSolverImplementation::~PetscNonlinearSolverImplementation(void)
{
  PetscErrorCode ierr;
  try  {
    PetscBool ok;
    ierr = PetscInitialized(&ok);
    if (ok) {
      ierr = SNESDestroy(&p_snes);
    }
  } catch (...) {
    // just eat it
  }
}

// -------------------------------------------------------------
// PetscNonlinearSolverImplementation::p_build
// -------------------------------------------------------------
void
PetscNonlinearSolverImplementation::p_build(const std::string& option_prefix)
{
  PetscErrorCode ierr(0);
  try {
    ierr  = SNESCreate(this->communicator(), &p_snes); CHKERRXX(ierr);
    p_petsc_F = PETScVector(*p_F);

    if (!p_function.empty()) {
      ierr = SNESSetFunction(p_snes, *p_petsc_F, FormFunction, 
                             static_cast<void *>(this)); CHKERRXX(ierr);
    }

    p_petsc_J = PETScMatrix(*p_J);
    
    if (!p_jacobian.empty()) {
      ierr = SNESSetJacobian(p_snes, *p_petsc_J, *p_petsc_J, FormJacobian, 
                             static_cast<void *>(this)); CHKERRXX(ierr);
    }

    // set the 
    ierr = SNESSetOptionsPrefix(p_snes, option_prefix.c_str()); CHKERRXX(ierr);
    KSP ksp;
    ierr = SNESGetKSP(p_snes, &ksp); CHKERRXX(ierr);
    ierr = KSPSetOptionsPrefix(ksp, option_prefix.c_str()); CHKERRXX(ierr);
    
    PC pc;
    ierr = KSPGetPC(ksp, &pc); CHKERRXX(ierr);
    ierr = PCSetOptionsPrefix(pc, option_prefix.c_str()); CHKERRXX(ierr);

    ierr = SNESMonitorSet(p_snes, MonitorNorms, PETSC_NULL, PETSC_NULL); CHKERRXX(ierr);

    ierr = SNESSetTolerances(p_snes, 
                             p_functionTolerance, 
                             PETSC_DEFAULT,
                             p_solutionTolerance,
                             p_maxIterations, 
                             PETSC_DEFAULT);
                             
    ierr = SNESSetFromOptions(p_snes); CHKERRXX(ierr);
    
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscNonlinearSolverImplementation::p_configure
// -------------------------------------------------------------
void
PetscNonlinearSolverImplementation::p_configure(utility::Configuration::CursorPtr props)
{
  NonlinearSolverImplementation::p_configure(props);
  this->build(props);
}

// -------------------------------------------------------------
// PetscNonlinearSolverImplementation::FormJacobian
// -------------------------------------------------------------

#if PETSC_VERSION_LT(3,5,0)

PetscErrorCode 
PetscNonlinearSolverImplementation::FormJacobian(SNES snes, Vec x, 
                                                 Mat *jac, Mat *B, 
                                                 MatStructure *flag, 
                                                 void *dummy)
{
  PetscErrorCode ierr(0);

  // Necessary C cast
  PetscNonlinearSolverImplementation *solver =
    (PetscNonlinearSolverImplementation *)dummy;

  // Copy PETSc's current estimate into 

  // Should be the case, but just make sure
  BOOST_ASSERT(*jac == *solver->p_petsc_J);
  BOOST_ASSERT(*B == *solver->p_petsc_J);

  // Not sure about this
  BOOST_ASSERT(x == *(solver->p_petsc_X));

  // May need to do this, which seems slow.
  // ierr = VecCopy(x, *(solver->p_petsc_X)); CHKERRQ(ierr);

  // Call the user-specified function (object) to form the Jacobian
  (solver->p_jacobian)(*(solver->p_X), *(solver->p_J));

  *flag = SAME_NONZERO_PATTERN;

  return ierr;
}

#else

PetscErrorCode 
PetscNonlinearSolverImplementation::FormJacobian(SNES snes, Vec x, 
                                                 Mat jac, Mat B, 
                                                 void *dummy)
{
  PetscErrorCode ierr(0);

  // Necessary C cast
  PetscNonlinearSolverImplementation *solver =
    (PetscNonlinearSolverImplementation *)dummy;

  // Copy PETSc's current estimate into 

  // Should be the case, but just make sure
  BOOST_ASSERT(jac == *(solver->p_petsc_J));
  BOOST_ASSERT(B == *(solver->p_petsc_J));

  // Not sure about this
  BOOST_ASSERT(x == *(solver->p_petsc_X));

  // May need to do this, which seems slow.
  // ierr = VecCopy(x, *(solver->p_petsc_X)); CHKERRQ(ierr);

  // Call the user-specified function (object) to form the Jacobian
  (solver->p_jacobian)(*(solver->p_X), *(solver->p_J));

  return ierr;
}

#endif

// -------------------------------------------------------------
// PetscNonlinearSolverImplementation::FormFunction
// -------------------------------------------------------------
PetscErrorCode
PetscNonlinearSolverImplementation::FormFunction(SNES snes, Vec x, Vec f, void *dummy)
{
  PetscErrorCode ierr(0);

  // Necessary C cast
  PetscNonlinearSolverImplementation *solver =
    (PetscNonlinearSolverImplementation *)dummy;

  // Apparently, you cannot count on x and f (like in FormJacobian())
  // being the same as those used to set up SNES, so the PETSc
  // solution and function vectors are wrapped temporarily for the
  // user function

  boost::scoped_ptr<Vector> 
    xtmp(new Vector(new PETScVectorImplementation(x, false)));

  boost::scoped_ptr<Vector> 
    ftmp(new Vector(new PETScVectorImplementation(f, false)));

  // Call the user-specified function (object) to form the RHS
  (solver->p_function)(*xtmp, *ftmp);

  xtmp.reset();
  ftmp.reset();

  return ierr;
}
  

// -------------------------------------------------------------
// PetscNonlinearSolverImplementation::p_solve
// -------------------------------------------------------------
void
PetscNonlinearSolverImplementation::p_solve(void)
{
  PetscErrorCode ierr(0);
  p_petsc_X = PETScVector(*p_X);
  int me(this->processor_rank());

  try {
    ierr = SNESSolve(p_snes, NULL, *p_petsc_X); CHKERRXX(ierr);
    SNESConvergedReason reason;
    PetscInt iter;
    ierr = SNESGetConvergedReason(p_snes, &reason); CHKERRXX(ierr);
    ierr = SNESGetIterationNumber(p_snes, &iter); CHKERRXX(ierr);

    std::string msg;
    if (reason < 0) {
      msg = 
        boost::str(boost::format("%d: PETSc SNES diverged after %d iterations, reason: %d") % 
                   me % iter % reason);
      throw Exception(msg);
    } else if (me == 0) {
      msg = 
        boost::str(boost::format("%d: PETSc SNES converged after %d iterations, reason: %d") % 
                   me % iter % reason);
      std::cerr << msg << std::endl;
    }
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  } catch (const Exception& e) {
    throw e;
  }
}

} // namespace math
} // namespace gridpack
