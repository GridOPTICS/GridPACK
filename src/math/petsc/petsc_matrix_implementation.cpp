/**
 * @file   petsc_matrix_implementation.cpp
 * @author William A. Perkins
 * @date   2013-05-16 10:52:29 d3g096
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
PETScMatrixImplementation::PETScMatrixImplementation(const parallel::Communicator& comm,
                                                     const int& rows, const int& cols,
                                                     const bool& dense)
  : MatrixImplementation(comm)
{
  PetscErrorCode ierr(0);
  try {
    ierr = MatCreate(this->communicator(), &p_matrix); CHKERRXX(ierr);
    MatType the_type;
    if (dense) {
      if (this->communicator().size() == 1) {
        the_type = MATSEQDENSE;
      } else {
        the_type = MATDENSE;
      }
    } else {
      if (this->communicator().size() == 1) {
        the_type = MATSEQAIJ;
      } else {
        the_type = MATMPIAIJ;
      }
    }
    ierr = MatSetType(p_matrix, the_type); CHKERRXX(ierr);
    ierr = MatSetSizes(p_matrix,
                       PETSC_DECIDE, PETSC_DECIDE,
                       rows, cols); CHKERRXX(ierr);
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
// PETScMatrixImplementation::p_rows
// -------------------------------------------------------------
int 
PETScMatrixImplementation::p_rows(void) const
{
  PetscErrorCode ierr(0);
  int result(0);
  try {
    PetscInt rows, cols;
    ierr = MatGetSize(p_matrix, &rows, &cols);
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
    PetscInt rows, cols;
    ierr = MatGetLocalSize(p_matrix, &rows, &cols);
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
    PetscInt rows, cols;
    ierr = MatGetLocalSize(p_matrix, &rows, &cols);
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
    ierr = MatAssemblyBegin(p_matrix, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(p_matrix, MAT_FINAL_ASSEMBLY);
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
