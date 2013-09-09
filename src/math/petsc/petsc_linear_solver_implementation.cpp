/**
 * @file   petsc_linear_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2013-09-09 12:24:46 d3g096
 * 
 * @brief  
 * 
 * 
 */


#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_matrix_implementation.hpp"
#include "petsc/petsc_linear_solver_implementation.hpp"
#include "petsc/petsc_matrix_extractor.hpp"
#include "petsc/petsc_vector_extractor.hpp"

namespace gridpack {
namespace math {


// -------------------------------------------------------------
//  class PETScLinearSolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// PETScLinearSolverImplementation:: constructors / destructor
// -------------------------------------------------------------
PETScLinearSolverImplementation::PETScLinearSolverImplementation(const Matrix& A)
  : LinearSolverImplementation(A)
{
  PetscErrorCode ierr;
  try  {
    ierr = KSPCreate(this->communicator(), &p_KSP); CHKERRXX(ierr);
    ierr = KSPSetInitialGuessNonzero(p_KSP,PETSC_TRUE); CHKERRXX(ierr); 
    ierr = KSPSetFromOptions(p_KSP);CHKERRXX(ierr);
    p_set_matrix(*p_A);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

PETScLinearSolverImplementation::~PETScLinearSolverImplementation(void)
{
  PetscErrorCode ierr;
  try  {
    PetscBool ok;
    ierr = PetscInitialized(&ok);
    if (ok) {
      ierr = KSPDestroy(&p_KSP);
    }
  } catch (...) {
    // just eat it
  }
}

// -------------------------------------------------------------
// PETScLinearSolverImplementation::p_accept
// -------------------------------------------------------------
void 
PETScLinearSolverImplementation::p_accept(ImplementationVisitor& visitor)
{
  visitor.visit(*this);
}

void 
PETScLinearSolverImplementation::p_accept(ConstImplementationVisitor& visitor) const
{
  visitor.visit(*this);
}

// -------------------------------------------------------------
// PETScLinearSolverImplementation::p_solve
// -------------------------------------------------------------
void
PETScLinearSolverImplementation::p_solve(const Vector& b, Vector& x) const
{
  PetscErrorCode ierr(0);
  try {
    const Vec *bvec(PETScVector(b));
    Vec *xvec(PETScVector(x));
    ierr = KSPSolve(p_KSP, *bvec, *xvec); CHKERRXX(ierr);
    int its;
    KSPConvergedReason reason;
    PetscReal rnorm;
    ierr = KSPGetIterationNumber(p_KSP, &its); CHKERRXX(ierr);
    ierr = KSPGetConvergedReason(p_KSP, &reason); CHKERRXX(ierr);
    ierr = KSPGetResidualNorm(p_KSP, &rnorm); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScLinearSolverImplementation::p_set_matrix
// -------------------------------------------------------------
void
PETScLinearSolverImplementation::p_set_matrix(const Matrix& A)
{
  PetscErrorCode ierr(0);
  try  {
    Mat *Amat(PETScMatrix(*p_A));
    ierr = KSPSetOperators(p_KSP, *Amat, *Amat, SAME_NONZERO_PATTERN); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}  


} // namespace math
} // namespace gridpack
