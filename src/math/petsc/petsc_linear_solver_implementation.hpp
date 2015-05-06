// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_linear_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2015-03-05 12:54:57 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_linear_solver_implementation_hpp_
#define _petsc_linear_solver_implementation_hpp_

#include <boost/format.hpp>

#include <petscksp.h>
#include "petsc/petsc_exception.hpp"
#include "linear_solver_implementation.hpp"
#include "petsc_configurable.hpp"
#include "petsc/petsc_matrix_extractor.hpp"
#include "petsc/petsc_vector_extractor.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScLinearSolverImplementation
// -------------------------------------------------------------
template <typename T, typename I>
class PETScLinearSolverImplementation 
  : public LinearSolverImplementation<T, I>,
    private PETScConfigurable

{
public:

  typedef typename LinearSolverImplementation<T, I>::MatrixType MatrixType;
  typedef typename LinearSolverImplementation<T, I>::VectorType VectorType;

  /// Default constructor.
  PETScLinearSolverImplementation(MatrixType& A)
    : LinearSolverImplementation<T, I>(A.communicator()),
      PETScConfigurable(this->communicator()),
      p_A(PETScMatrix(A))
  {
  }

  /// Destructor
  ~PETScLinearSolverImplementation(void)
  {
    PetscErrorCode ierr(0);
    try  {
      PetscBool ok;
      ierr = PetscInitialized(&ok);
      if (ok) {
        ierr = KSPDestroy(&p_KSP); CHKERRXX(ierr);
      }
    } catch (...) {
      // just eat it
    }
  }


protected:

  /// The coefficient Matrix
  Mat *p_A;

  /// The PETSc linear solver 
  KSP p_KSP;

  /// Do what is necessary to build this instance
  void p_build(const std::string& option_prefix)
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
                              LinearSolverImplementation<T, I>::p_relativeTolerance, 
                              LinearSolverImplementation<T, I>::p_solutionTolerance, 
                              PETSC_DEFAULT,
                              LinearSolverImplementation<T, I>::p_maxIterations); CHKERRXX(ierr);

      ierr = KSPSetFromOptions(p_KSP);CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }  

  /// Solve w/ the specified RHS and estimate (result in x)
  void p_solve(const VectorType& b, VectorType& x) const
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

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props)
  {
    LinearSolverImplementation<T, I>::p_configure(props);
    this->build(props);
  }

};

} // namespace math
} // namespace gridpack

#endif
