/**
 * @file   petsc_nonlinear_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2013-10-01 14:41:07 d3g096
 * 
 * @brief  
 * 
 * 
 */

#include <boost/format.hpp>

#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_nonlinear_solver_implementation.hpp"
#include "petsc/petsc_matrix_implementation.hpp"
#include "petsc/petsc_matrix_extractor.hpp"
#include "petsc/petsc_vector_extractor.hpp"

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
    p_snes(), 
    p_petsc_J(), p_petsc_F(),
    p_petsc_X()                 // set by p_solve()
{
  
  PetscErrorCode ierr(0);
  try {
    ierr  = SNESCreate(this->communicator(), &p_snes); CHKERRXX(ierr);
    p_petsc_F = PETScVector(*p_F);

    if (!form_function.empty()) {
      ierr = SNESSetFunction(p_snes, *p_petsc_F, FormFunction, this); CHKERRXX(ierr);
    }

    p_petsc_J = PETScMatrix(*p_J);
    
    if (!form_jacobian.empty()) {
      ierr = SNESSetJacobian(p_snes, *p_petsc_J, *p_petsc_J, FormJacobian, this); CHKERRXX(ierr);
    }

    ierr = SNESSetFromOptions(p_snes); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
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
// PetscNonlinearSolverImplementation::FormJacobian
// -------------------------------------------------------------
PetscErrorCode 
PetscNonlinearSolverImplementation::FormJacobian(SNES snes, Vec x, Mat *jac, Mat *B, 
                                                 MatStructure *flag, void *dummy)
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
  // being the same as those used to set up SNES, so copying is
  // necessary
  ierr = VecCopy(x, *(solver->p_petsc_X)); CHKERRQ(ierr);

  // Call the user-specified function (object) to form the RHS
  (solver->p_function)(*(solver->p_X), *(solver->p_F));

  ierr = VecCopy(*(solver->p_petsc_F), f); CHKERRQ(ierr);

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

  try {
    ierr = SNESSolve(p_snes, NULL, *p_petsc_X); CHKERRXX(ierr);
    SNESConvergedReason reason;
    PetscInt iter;
    ierr = SNESGetConvergedReason(p_snes, &reason); CHKERRXX(ierr);
    ierr = SNESGetIterationNumber(p_snes, &iter); CHKERRXX(ierr);
    std::string msg;
    if (reason < 0) {
      msg = 
        boost::str(boost::format("PETSc SNES diverged after %d iterations, reason: %d") % 
                   iter % reason);
      throw Exception(msg);
    } else {
      msg = 
        boost::str(boost::format("PETSc SNES converged after %d iterations, reason: %d") % 
                   iter % reason);
      std::cerr << msg << std::endl;
    }
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  } catch (const Exception& e) {
    throw e;
  }
}

// -------------------------------------------------------------
// PetscNonlinearSolverImplementation::p_configure
// -------------------------------------------------------------
void
PetscNonlinearSolverImplementation::p_configure(utility::Configuration::Cursor *props)
{

}

} // namespace math
} // namespace gridpack
