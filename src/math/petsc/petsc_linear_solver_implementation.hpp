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
 * @date   2015-08-18 13:40:01 d3g096
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
    : LinearSolverImplementation<T, I>(A),
      PETScConfigurable(this->communicator()),
      p_matrixSet(false)
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

  /// The PETSc linear solver 
  KSP p_KSP;

  /// For constant matrices, has the coefficient matrix been set
  mutable bool p_matrixSet;

  /// Do what is necessary to build this instance
  void p_build(const std::string& option_prefix)
  {
    PetscErrorCode ierr;
    try  {
      parallel::Communicator comm(this->communicator());
      if (this->p_doSerial) {
        comm = this->communicator().self();
      }
      ierr = KSPCreate(comm, &p_KSP); CHKERRXX(ierr);
      if (!this->p_guessZero) {
        ierr = KSPSetInitialGuessNonzero(p_KSP,PETSC_TRUE); CHKERRXX(ierr); 
      } else {
        ierr = KSPSetInitialGuessNonzero(p_KSP,PETSC_FALSE); CHKERRXX(ierr); 
      }
      ierr = KSPSetOptionsPrefix(p_KSP, option_prefix.c_str()); CHKERRXX(ierr);

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
  void p_solveImpl(MatrixType& A, const VectorType& b, VectorType& x) const
  {
    PetscErrorCode ierr(0);
    try {
      Mat *Amat(PETScMatrix(A));

      if (p_matrixSet && this->p_constSerialMatrix) {
        // KSPSetOperators can be skipped
      } else {
#if PETSC_VERSION_LT(3,5,0)
        ierr = KSPSetOperators(p_KSP, *Amat, *Amat, SAME_NONZERO_PATTERN); CHKERRXX(ierr);
#else
        ierr = KSPSetOperators(p_KSP, *Amat, *Amat); CHKERRXX(ierr);
#endif
        p_matrixSet = true;
      }

      this->p_resolveImpl(b, x);
          
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    } catch (const Exception& e) {
      throw e;
    }
  }

  /// Solve again w/ the specified RHS, put result in specified vector (specialized)
  void p_resolveImpl(const VectorType& b, VectorType& x) const
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
