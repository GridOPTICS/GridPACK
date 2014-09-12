// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_linear_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2014-09-12 13:51:08 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/format.hpp>

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
PETScLinearSolverImplementation::PETScLinearSolverImplementation(Matrix& A)
  : LinearSolverImplementation(A.communicator()),
    PETScConfigurable(this->communicator()),
    p_A(PETScMatrix(A))
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
    ierr = PCSetOptionsPrefix(pc, option_prefix.c_str()); CHKERRXX(ierr);

    ierr = KSPSetTolerances(p_KSP, 
                            p_relativeTolerance, 
                            p_solutionTolerance, 
                            PETSC_DEFAULT,
                            p_maxIterations); CHKERRXX(ierr);

    ierr = KSPSetFromOptions(p_KSP);CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}  

// -------------------------------------------------------------
// PETScLinearSolverImplementation::p_configure
// -------------------------------------------------------------
void
PETScLinearSolverImplementation::p_configure(utility::Configuration::CursorPtr props)
{
  LinearSolverImplementation::p_configure(props);
  this->build(props);
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

#if PETSC_VERSION_LT(3,5,0)
    ierr = KSPSetOperators(p_KSP, *p_A, *p_A, SAME_NONZERO_PATTERN); CHKERRXX(ierr);
#else
    ierr = KSPSetOperators(p_KSP, *p_A, *p_A); CHKERRXX(ierr);
#endif

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
    } else if (me == 0) {
      msg = 
        boost::str(boost::format("%d: PETSc KSP converged after %d iterations, reason: %d") % 
                   me % its % reason);
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
