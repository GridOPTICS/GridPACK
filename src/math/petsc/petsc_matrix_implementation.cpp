/**
 * @file   petsc_matrix_implementation.cpp
 * @author William A. Perkins
 * @date   2013-05-20 08:47:27 d3g096
 * 
 * @brief  PETSc-specific matrix implementation
 * 
 * 
 */

#include "gridpack/math/implementation_visitor.hpp"
#include "gridpack/math/petsc/petsc_matrix_implementation.hpp"
#include "gridpack/math/petsc/petsc_exception.hpp"
#include "gridpack/math/petsc/petsc_matrix_extractor.hpp"


namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScMatrixImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// PETScMatrixImplementation:: constructors / destructor
// -------------------------------------------------------------
/** 
 * 
 * 
 * @param comm 
 * @param local_rows rows to make local to this processor
 * @param cols 
 * @param dense 
 */
PETScMatrixImplementation::PETScMatrixImplementation(const parallel::Communicator& comm,
                                                     const int& local_rows, const int& cols,
                                                     const bool& dense)
  : MatrixImplementation(comm)
{
  PetscErrorCode ierr(0);
  static const PetscInt diagonal_non_zero_guess(10);
  static const PetscInt offdiagonal_non_zero_guess(2);
  try {
    ierr = MatCreate(this->communicator(), &p_matrix); CHKERRXX(ierr);
    MatType the_type;
    ierr = MatSetSizes(p_matrix,
                       local_rows, PETSC_DECIDE,
                       PETSC_DETERMINE, cols); CHKERRXX(ierr);
    if (dense) {
      if (this->communicator().size() == 1) {
        ierr = MatSetType(p_matrix, MATSEQDENSE); CHKERRXX(ierr);
        ierr = MatSeqDenseSetPreallocation(p_matrix, PETSC_NULL); CHKERRXX(ierr);
      } else {
        ierr = MatSetType(p_matrix, MATDENSE); CHKERRXX(ierr);
        ierr = MatMPIDenseSetPreallocation(p_matrix, PETSC_NULL); CHKERRXX(ierr);
      }
    } else {
      if (this->communicator().size() == 1) {
        ierr = MatSetType(p_matrix, MATSEQAIJ); CHKERRXX(ierr);
        ierr = MatSeqAIJSetPreallocation(p_matrix, 
                                         diagonal_non_zero_guess + offdiagonal_non_zero_guess,
                                         PETSC_NULL); CHKERRXX(ierr);
      } else {
        ierr = MatSetType(p_matrix, MATMPIAIJ); CHKERRXX(ierr);
        ierr = MatMPIAIJSetPreallocation(p_matrix, 
                                         diagonal_non_zero_guess,
                                         PETSC_NULL,
                                         offdiagonal_non_zero_guess, 
                                         PETSC_NULL); CHKERRXX(ierr);
      }
    }
    ierr = MatSetFromOptions(p_matrix); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

PETScMatrixImplementation::~PETScMatrixImplementation(void)
{
  PetscErrorCode ierr;
  try  {
    PetscBool ok;
    ierr = PetscInitialized(&ok);
    if (ok) {
      ierr = MatDestroy(&p_matrix);
    }
  } catch (...) {
    // just eat it
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_local_row_range
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_local_row_range(int& lo, int& hi) const
{
  PetscErrorCode ierr(0);
  try {
    PetscInt plo, phi;
    ierr = MatGetOwnershipRange(p_matrix, &plo, &phi); CHKERRXX(ierr);
    lo = plo;
    hi = phi;
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_rows
// -------------------------------------------------------------
int 
PETScMatrixImplementation::p_rows(void) const
{
  PetscErrorCode ierr(0);
  int result(0);
  try {
    PetscInt rows;
    ierr = MatGetSize(p_matrix, &rows, PETSC_NULL); CHKERRXX(ierr);
    result = rows;
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  return result;
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_local_rows
// -------------------------------------------------------------
int
PETScMatrixImplementation::p_local_rows(void) const
{
  PetscErrorCode ierr(0);
  int result(0);
  try {
    PetscInt rows;
    ierr = MatGetLocalSize(p_matrix, &rows, PETSC_NULL); CHKERRXX(ierr);
    result = rows;
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  return result;
}  

// -------------------------------------------------------------
// PETScMatrixImplementation::p_cols
// -------------------------------------------------------------
int
PETScMatrixImplementation::p_cols(void) const
{
  PetscErrorCode ierr(0);
  int result(0);
  try {
    PetscInt cols;
    ierr = MatGetSize(p_matrix, PETSC_NULL, &cols); CHKERRXX(ierr);
    result = cols;
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  return result;
}  

// -------------------------------------------------------------
// PETScMatrixImplementation::p_ready
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_ready(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = MatAssemblyBegin(p_matrix, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
    ierr = MatAssemblyEnd(p_matrix, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_accept
// -------------------------------------------------------------
void 
PETScMatrixImplementation::p_accept(ImplementationVisitor& visitor)
{
  visitor.visit(*this);
}

void 
PETScMatrixImplementation::p_accept(ConstImplementationVisitor& visitor) const
{
  visitor.visit(*this);
}

} // namespace math
} // namespace gridpack
