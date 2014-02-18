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
 * @date   2014-02-17 15:53:32 d3g096
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
  parallel::Communicator result(comm);
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
  p_build_matrix(comm, local_rows, global_cols);
  if (dense) {
    p_set_dense_matrix();
  } else {
    p_set_sparse_matrix();
  }
}

PETScMatrixImplementation::PETScMatrixImplementation(const parallel::Communicator& comm,
                                                     const int& local_rows, const int& cols,
                                                     const int& max_nonzero_per_row)
  : MatrixImplementation(comm),
    p_matrixWrapped(false)
{
  p_build_matrix(comm, local_rows, cols);
  p_set_sparse_matrix(max_nonzero_per_row);
}

PETScMatrixImplementation::PETScMatrixImplementation(const parallel::Communicator& comm,
                                                     const int& local_rows, const int& cols,
                                                     const int *nonzero_by_row)
  : MatrixImplementation(comm),
    p_matrixWrapped(false)
{
  p_build_matrix(comm, local_rows, cols);
  p_set_sparse_matrix(nonzero_by_row);
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
// PETScMatrixImplementation::p_build_matrix
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_build_matrix(const parallel::Communicator& comm,
                                          const int& local_rows, 
                                          const int& global_cols)
{
  PetscErrorCode ierr(0);
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
    // FIXME: ierr = MatSetSizes(p_matrix, lrows, lcols, grows, gcols); CHKERRXX(ierr);
    ierr = MatSetSizes(p_matrix, PETSC_DECIDE, PETSC_DECIDE, grows, gcols); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_set_dense_matrix
// -------------------------------------------------------------
void 
PETScMatrixImplementation::p_set_dense_matrix(void)
{
  PetscErrorCode ierr(0);
  try {
    if (this->communicator().size() == 1) {
      ierr = MatSetType(p_matrix, MATSEQDENSE); CHKERRXX(ierr);
      ierr = MatSeqDenseSetPreallocation(p_matrix, PETSC_NULL); CHKERRXX(ierr);
    } else {
      ierr = MatSetType(p_matrix, MATDENSE); CHKERRXX(ierr);
      ierr = MatMPIDenseSetPreallocation(p_matrix, PETSC_NULL); CHKERRXX(ierr);
    }
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_set_sparse_matrix
// -------------------------------------------------------------
void 
PETScMatrixImplementation::p_set_sparse_matrix(void)
{
  PetscErrorCode ierr(0);
  try {
    if (this->communicator().size() == 1) {
      ierr = MatSetType(p_matrix, MATSEQAIJ); CHKERRXX(ierr);
    } else {
      ierr = MatSetType(p_matrix, MATMPIAIJ); CHKERRXX(ierr);
    }
    ierr = MatSetFromOptions(p_matrix); CHKERRXX(ierr);
    ierr = MatSetUp(p_matrix); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

void 
PETScMatrixImplementation::p_set_sparse_matrix(const int& max_nz_per_row)
{
  PetscErrorCode ierr(0);
  const PetscInt diagonal_non_zero_guess(max_nz_per_row);
  PetscInt offdiagonal_non_zero_guess(static_cast<PetscInt>(0.5*diagonal_non_zero_guess));
  offdiagonal_non_zero_guess = std::max(offdiagonal_non_zero_guess, 10);

  try {
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
    ierr = MatSetFromOptions(p_matrix); CHKERRXX(ierr);
    ierr = MatSetUp(p_matrix); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

void 
PETScMatrixImplementation::p_set_sparse_matrix(const int *nz_by_row)
{
  std::vector<PetscInt> diagnz;
  int lrows(this->localRows());
  diagnz.reserve(lrows);
  std::copy(nz_by_row, nz_by_row+lrows, 
            std::back_inserter(diagnz));

  PetscErrorCode ierr(0);
  try {
    if (this->communicator().size() == 1) {
      ierr = MatSetType(p_matrix, MATSEQAIJ); CHKERRXX(ierr);
      ierr = MatSeqAIJSetPreallocation(p_matrix, 
                                       PETSC_DECIDE,
                                       &diagnz[0]); CHKERRXX(ierr);
    } else {
      std::vector<PetscInt> offdiagnz;
      offdiagnz.reserve(lrows);
      std::copy(nz_by_row, nz_by_row+lrows, 
                std::back_inserter(offdiagnz));
      ierr = MatSetType(p_matrix, MATMPIAIJ); CHKERRXX(ierr);
      ierr = MatMPIAIJSetPreallocation(p_matrix, 
                                       PETSC_DECIDE,
                                       &diagnz[0],
                                       PETSC_DECIDE, 
                                       &offdiagnz[0]); CHKERRXX(ierr);
    }
    ierr = MatSetFromOptions(p_matrix); CHKERRXX(ierr);
    ierr = MatSetUp(p_matrix); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
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
    // MatInfo info;
    // ierr = MatGetInfo(p_matrix,MAT_LOCAL,&info);
    // std::cerr << this->processor_rank() << ": Matrix::ready(): "
    //           << "size = (" << this->p_rows() << "x" << this->p_cols() << "), "
    //           << "#assembled = " << info.assemblies << ", "
    //           << "#mallocs = " << info.mallocs << ", "
    //           << "#nz allocated = " << info.nz_allocated << ", "
    //           << "#nz used = " << info.nz_used << ", "
    //           << std::endl;
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
