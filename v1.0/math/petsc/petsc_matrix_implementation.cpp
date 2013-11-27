// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_implementation.cpp
 * @author William A. Perkins
 * @date   2013-11-12 10:27:53 d3g096
 * 
 * @brief  PETSc-specific matrix implementation
 * 
 * 
 */
// -------------------------------------------------------------

#include "matrix.hpp"
#include "implementation_visitor.hpp"
#include "petsc/petsc_matrix_implementation.hpp"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_matrix_extractor.hpp"


namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScMatrixImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// PETScMatrixImplementation::p_getCommunicator
// -------------------------------------------------------------
parallel::Communicator
PETScMatrixImplementation::p_getCommunicator(const Mat& m)
{
  MPI_Comm comm(PetscObjectComm((PetscObject)m));
  parallel::Communicator result(comm, boost::mpi::comm_attach);
  return result;
}


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
                                                     const int& local_rows, const int& global_cols,
                                                     const bool& dense)
  : MatrixImplementation(comm),
    p_matrixWrapped(false)
{
  PetscErrorCode ierr(0);
  static const PetscInt diagonal_non_zero_guess(10);
  static const PetscInt offdiagonal_non_zero_guess(10);
  try {

    // If any ownership arguments are specifed, *all* ownership
    // arguments need to be specified.  AND, if the matrix is square,
    // the local rows and cols need to be the same.

    PetscInt lrows(local_rows), grows(PETSC_DECIDE);
    ierr = PetscSplitOwnership(comm, &lrows, &grows); CHKERRXX(ierr);

    PetscInt lcols(PETSC_DECIDE), gcols(global_cols);
    if (grows == global_cols) {
      lcols = lrows;
    } else {
      ierr = PetscSplitOwnership(comm, &lcols, &gcols); CHKERRXX(ierr);
    }
    
    ierr = MatCreate(this->communicator(), &p_matrix); CHKERRXX(ierr);
    MatType the_type;
    ierr = MatSetSizes(p_matrix, lrows, lcols, grows, gcols); CHKERRXX(ierr);
    if (dense) {
      if (this->communicator().size() == 1) {
        ierr = MatSetType(p_matrix, MATSEQDENSE); CHKERRXX(ierr);
        // ierr = MatSeqDenseSetPreallocation(p_matrix, PETSC_NULL); CHKERRXX(ierr);
      } else {
        ierr = MatSetType(p_matrix, MATDENSE); CHKERRXX(ierr);
        // ierr = MatMPIDenseSetPreallocation(p_matrix, PETSC_NULL); CHKERRXX(ierr);
      }
    } else {
      if (this->communicator().size() == 1) {
        ierr = MatSetType(p_matrix, MATSEQAIJ); CHKERRXX(ierr);
        // ierr = MatSeqAIJSetPreallocation(p_matrix, 
        //                                  diagonal_non_zero_guess + offdiagonal_non_zero_guess,
        //                                  PETSC_NULL); CHKERRXX(ierr);
      } else {
        ierr = MatSetType(p_matrix, MATMPIAIJ); CHKERRXX(ierr);
        // ierr = MatMPIAIJSetPreallocation(p_matrix, 
        //                                  diagonal_non_zero_guess,
        //                                  PETSC_NULL,
        //                                  offdiagonal_non_zero_guess, 
        //                                  PETSC_NULL); CHKERRXX(ierr);
      }
      // By default, new elements that are not pre-allocated cause an
      // error. Let's disable that. With the preallocation above, it
      // should not happen very often.
      // ierr = MatSetOption(p_matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); CHKERRXX(ierr);
      // ierr = MatSetOption(p_matrix, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE); CHKERRXX(ierr)
                                                                                  ;
    }
    ierr = MatSetFromOptions(p_matrix); CHKERRXX(ierr);
    ierr = MatSetUp(p_matrix); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

PETScMatrixImplementation::PETScMatrixImplementation(Mat& m, const bool& copyMat)
  : MatrixImplementation(p_getCommunicator(m)),
    p_matrixWrapped(false)
{
  PetscErrorCode ierr;
  try {

    if (copyMat) {
      ierr = MatDuplicate(m, MAT_COPY_VALUES, &p_matrix); CHKERRXX(ierr);
    } else {
      p_matrix = m;
      p_matrixWrapped = true;
    }

  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

PETScMatrixImplementation::~PETScMatrixImplementation(void)
{
  PetscErrorCode ierr;
  if (!p_matrixWrapped) {
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
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_local_row_range
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_localRowRange(int& lo, int& hi) const
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
PETScMatrixImplementation::p_localRows(void) const
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
// PETScMatrixImplementation::p_setElement
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_setElement(const int& i, const int& j, 
                                         const ComplexType& x)
{
  PetscErrorCode ierr(0);
  try {
    ierr = MatSetValue(p_matrix, i, j, x, INSERT_VALUES); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_set_elements
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_setElements(const int& n, 
                                          const int *i, const int *j, 
                                          const ComplexType *x)
{
  // FIXME: There's probably a better way
  for (int k = 0; k < n; k++) {
    this->p_setElement(i[k], j[k], x[k]);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_add_element
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_addElement(const int& i, const int& j, 
                                         const ComplexType& x)
{
  PetscErrorCode ierr(0);
  try {
    ierr = MatSetValue(p_matrix, i, j, x, ADD_VALUES); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_add_elements
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_addElements(const int& n, 
                                          const int *i, const int *j, 
                                          const ComplexType *x)
{
  // FIXME: There's probably a better way
  for (int k = 0; k < n; k++) {
    this->p_addElement(i[k], j[k], x[k]);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_get_element
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_getElement(const int& i, const int& j, 
                                         ComplexType& x) const
{
  PetscErrorCode ierr(0);
  try {
    static const int one(1);
    PetscInt iidx[one] = { i };
    PetscInt jidx[one] = { j };
    ierr = MatGetValues(p_matrix, one, &iidx[0], one, &jidx[0], &x); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_get_elements
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_getElements(const int& n,
                                          const int *i, const int *j, 
                                          ComplexType *x) const
{
  // FIXME: There is a better way
  for (int k = 0; k < n; k++) {
    this->p_getElement(i[k], j[k], x[k]);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_real
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_real(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = MatRealPart(p_matrix); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_imaginary
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_imaginary(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = MatImaginaryPart(p_matrix); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_conjugate
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_conjugate(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = MatConjugate(p_matrix); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_norm2
// -------------------------------------------------------------
ComplexType
PETScMatrixImplementation::p_norm2(void) const
{
  PetscErrorCode ierr(0);
  ComplexType result(0.0);
  try {
    PetscReal norm;
    ierr = MatNorm(p_matrix, NORM_FROBENIUS, &norm); CHKERRXX(ierr);
    result = norm;
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

// -------------------------------------------------------------
// PETScMatrixImplementation::p_clone
// -------------------------------------------------------------
MatrixImplementation *
PETScMatrixImplementation::p_clone(void) const
{
  PETScMatrixImplementation *result =
    new PETScMatrixImplementation(const_cast<Mat&>(this->p_matrix), true);
  return result;
}

} // namespace math
} // namespace gridpack
