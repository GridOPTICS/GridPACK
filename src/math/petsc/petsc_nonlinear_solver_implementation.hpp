// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_nonlinear_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2015-03-26 14:25:10 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/format.hpp>
#include <petscsys.h>
#include <petscsnes.h>
#include "petsc/petsc_exception.hpp"
#include "nonlinear_solver_implementation.hpp"
#include "petsc/petsc_matrix_extractor.hpp"
#include "petsc/petsc_vector_extractor.hpp"
#include "petsc/petsc_vector_implementation.hpp"
#include "petsc/petsc_configurable.hpp"

namespace gridpack {
namespace math {

extern PetscErrorCode  
MonitorNorms(SNES snes, PetscInt its, PetscReal fgnorm, void *dummy);

// -------------------------------------------------------------
//  class PetscNonlinearSolverImplementation
// -------------------------------------------------------------
template <typename T, typename I = int>
class PetscNonlinearSolverImplementation 
  : public NonlinearSolverImplementation<T, I>,
    private PETScConfigurable
{
public:

  typedef typename NonlinearSolverImplementation<T, I>::VectorType VectorType;
  typedef typename NonlinearSolverImplementation<T, I>::MatrixType MatrixType;
  typedef typename NonlinearSolverImplementation<T, I>::JacobianBuilder JacobianBuilder;
  typedef typename NonlinearSolverImplementation<T, I>::FunctionBuilder FunctionBuilder;

  /// Default constructor.
  PetscNonlinearSolverImplementation(const parallel::Communicator& comm,
                                     const int& local_size,
                                     JacobianBuilder form_jacobian,
                                     FunctionBuilder form_function)
    : NonlinearSolverImplementation<T, I>(comm, local_size, form_jacobian, form_function),
      PETScConfigurable(this->communicator()),
      p_snes(), 
      p_petsc_J(), p_petsc_F(),
      p_petsc_X()                 // set by p_solve()
  {
    
  }


  /// Construct with an existing Jacobian Matrix
  PetscNonlinearSolverImplementation(MatrixType& J,
                                     JacobianBuilder form_jacobian,
                                     FunctionBuilder form_function)
    : NonlinearSolverImplementation<T, I>(J, form_jacobian, form_function),
      PETScConfigurable(this->communicator()),
      p_snes(), 
      p_petsc_J(), p_petsc_F(),
      p_petsc_X()                 // set by p_solve()
  {
    
  }

  /// Destructor
  ~PetscNonlinearSolverImplementation(void)
  {
    PetscErrorCode ierr(0);
    try  {
      PetscBool ok;
      ierr = PetscInitialized(&ok); CHKERRXX(ierr);
      if (ok) {
        ierr = SNESDestroy(&p_snes); CHKERRXX(ierr);
      }
    } catch (...) {
      // just eat it
    }
  }

protected:

  /// The PETSc nonlinear solver instance
  SNES p_snes;

  /// A pointer to PETSc matrix part of ::p_J
  Mat *p_petsc_J;

  /// A pointer to the PETSc vector part of ::p_F
  Vec *p_petsc_F;

  /// A pointer to the PETSc vector part of ::p_X
  Vec *p_petsc_X;

  /// Do what is necessary to build this instance
  void p_build(const std::string& option_prefix)
  {
    PetscErrorCode ierr(0);
    try {
      ierr  = SNESCreate(this->communicator(), &p_snes); CHKERRXX(ierr);
      p_petsc_F = PETScVector(*(this->p_F));

      if (!this->p_function.empty()) {
        ierr = SNESSetFunction(p_snes, *p_petsc_F, FormFunction, 
                               static_cast<void *>(this)); CHKERRXX(ierr);
      }

      p_petsc_J = PETScMatrix(*(this->p_J));
    
      if (!this->p_jacobian.empty()) {
        ierr = SNESSetJacobian(p_snes, *p_petsc_J, *p_petsc_J, FormJacobian, 
                               static_cast<void *>(this)); CHKERRXX(ierr);
      }

      // set the 
      ierr = SNESSetOptionsPrefix(p_snes, option_prefix.c_str()); CHKERRXX(ierr);

      ierr = SNESMonitorSet(p_snes, MonitorNorms, PETSC_NULL, PETSC_NULL); CHKERRXX(ierr);

      ierr = SNESSetTolerances(p_snes, 
                               this->p_functionTolerance, 
                               PETSC_DEFAULT,
                               this->p_solutionTolerance,
                               this->p_maxIterations, 
                               PETSC_DEFAULT);
                             
      ierr = SNESSetFromOptions(p_snes); CHKERRXX(ierr);
    
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }


  /// Solve w/ using the specified initial guess (specialized)
  void p_solve(VectorType& x)
  {
    NonlinearSolverImplementation<T, I>::p_solve(x);

    PetscErrorCode ierr(0);
    p_petsc_X = PETScVector(*(this->p_X));
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


  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props)
  {
    NonlinearSolverImplementation<T, I>::p_configure(props);
    this->build(props);
  }


#if PETSC_VERSION_LT(3,5,0)
  /// Routine to assemble Jacobian that is sent to PETSc
  static PetscErrorCode FormJacobian(SNES snes, Vec x, Mat *jac, Mat *B, 
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

#else
  /// Routine to assemble Jacobian that is sent to PETSc
  static PetscErrorCode FormJacobian(SNES snes, Vec x, Mat jac, Mat B, 
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

  /// Routine to assemble RHS that is sent to PETSc
  static PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *dummy)
  {
    PetscErrorCode ierr(0);

    // Necessary C cast
    PetscNonlinearSolverImplementation *solver =
      (PetscNonlinearSolverImplementation *)dummy;

    // Apparently, you cannot count on x and f (like in FormJacobian())
    // being the same as those used to set up SNES, so the PETSc
    // solution and function vectors are wrapped temporarily for the
    // user function

    boost::scoped_ptr< VectorType > 
      xtmp(new VectorType(new PETScVectorImplementation<T, I>(x, false)));

    boost::scoped_ptr< VectorType > 
      ftmp(new VectorType(new PETScVectorImplementation<T, I>(f, false)));

    // Call the user-specified function (object) to form the RHS
    (solver->p_function)(*xtmp, *ftmp);

    xtmp.reset();
    ftmp.reset();

    return ierr;
  }

};



} // namespace math
} // namespace gridpack
