/**
 * @file   petsc_linear_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2013-10-03 13:51:29 d3g096
 * 
 * @brief  
 * 
 * 
 */

#include <boost/format.hpp>

#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_matrix_implementation.hpp"
#include "petsc/petsc_linear_solver_implementation.hpp"
#include "petsc/petsc_matrix_extractor.hpp"
#include "petsc/petsc_vector_extractor.hpp"
#include "petsc/petsc_configuration.hpp"

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
// PETScLinearSolverImplementation::p_build
// -------------------------------------------------------------
void
PETScLinearSolverImplementation::p_build(const std::string& option_prefix)
{
  PetscErrorCode ierr;
  try  {
    ierr = KSPCreate(this->communicator(), &p_KSP); CHKERRXX(ierr);
    ierr = KSPSetInitialGuessNonzero(p_KSP,PETSC_TRUE); CHKERRXX(ierr); 
    ierr = KSPSetOptionsPrefix(p_KSP, option_prefix.c_str()); CHKERRXX(ierr);
    PC pc;
    ierr = KSPGetPC(p_KSP, &pc); CHKERRXX(ierr);
    ierr = KSPSetFromOptions(p_KSP);CHKERRXX(ierr);
    ierr = PCSetOptionsPrefix(pc, option_prefix.c_str()); CHKERRXX(ierr);
    p_set_matrix(*p_A);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}  

// -------------------------------------------------------------
// PETScLinearSolverImplementation::p_configure
// -------------------------------------------------------------
void
PETScLinearSolverImplementation::p_configure(utility::Configuration::Cursor *props)
{
  std::string prefix(petscProcessOptions(this->communicator(), props));
  p_build(prefix);
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
  int me(this->processor_rank());
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
    std::string msg;
    if (reason < 0) {
      msg = 
        boost::str(boost::format("%d: PETSc KSP diverged after %d iterations, reason: %d") % 
                   me % its % reason);
      throw Exception(msg);
    } else {
      msg = 
        boost::str(boost::format("%d: PETSc KSP converged after %d iterations, reason: %d") % 
                   me % its % reason);
      std::cerr << msg << std::endl;
    }
    
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  } catch (const Exception& e) {
    throw e;
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
