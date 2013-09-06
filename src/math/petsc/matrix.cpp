/**
 * @file   matrix.cpp
 * @author William A. Perkins
 * @date   2013-09-06 11:09:19 d3g096
 * 
 * @brief  PETSc specific part of Matrix
 * 
 * 
 */

#include "boost/assert.hpp"
#include "matrix.hpp"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_matrix_implementation.hpp"
#include "petsc/petsc_matrix_extractor.hpp"


namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class Matrix
// -------------------------------------------------------------

// -------------------------------------------------------------
// Matrix:: constructors / destructor
// -------------------------------------------------------------
Matrix::Matrix(const parallel::Communicator& comm,
               const int& local_rows,
               const int& cols,
               const StorageType& storage_type)
  : parallel::WrappedDistributed(), utility::Uncopyable(),
    p_matrix_impl()
{
  switch (storage_type) {
  case Sparse:
    p_matrix_impl.reset(new PETScMatrixImplementation(comm,
                                                      local_rows, cols, false));
    break;
  case Dense:
    p_matrix_impl.reset(new PETScMatrixImplementation(comm,
                                                      local_rows, cols, true));
    break;
  default:
    BOOST_ASSERT(false);
  }
  BOOST_ASSERT(p_matrix_impl);
  p_set_distributed(p_matrix_impl.get());
}

// -------------------------------------------------------------
// Matrix::equate
// -------------------------------------------------------------
void
Matrix::equate(const Matrix& B)
{
  this->p_check_compatible(B);
  Mat *pA(PETScMatrix(*this));
  const Mat *pB(PETScMatrix(B));

  PetscErrorCode ierr(0);
  try {
    PetscScalar one(1.0);
    ierr = MatCopy(*pB, *pA, DIFFERENT_NONZERO_PATTERN); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// Matrix::scale
// -------------------------------------------------------------
void
Matrix::scale(const ComplexType& xin)
{
  Mat *pA(PETScMatrix(*this));

  PetscErrorCode ierr(0);
  
  try {
    PetscScalar x(xin);
    ierr = MatScale(*pA, x); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}
    

// -------------------------------------------------------------
// Matrix::add
// -------------------------------------------------------------
void
Matrix::add(const Matrix& B)
{
  this->p_check_compatible(B);
  Mat *pA(PETScMatrix(*this));
  const Mat *pB(PETScMatrix(B));

  PetscErrorCode ierr(0);
  try {
    PetscScalar one(1.0);
    ierr = MatAXPY(*pA, one, *pB, DIFFERENT_NONZERO_PATTERN); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// Matrix::identity
// -------------------------------------------------------------
void
Matrix::identity(void)
{
  this->zero();
  Mat *pA(PETScMatrix(*this));

  PetscErrorCode ierr(0);
  try {
    ierr = MatZeroEntries(*pA); CHKERRXX(ierr);
    PetscScalar one(1.0);
    ierr = MatShift(*pA, one);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// Matrix::zero
// -------------------------------------------------------------
void
Matrix::zero(void)
{
  Mat *pA(PETScMatrix(*this));

  PetscErrorCode ierr(0);
  try {
    ierr = MatZeroEntries(*pA); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// petsc_make_viewer
// -------------------------------------------------------------
static void
petsc_make_viewer(const char* filename, PetscViewer *viewer)
{
  PetscErrorCode ierr;
  if (filename != NULL) {
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, viewer); ; CHKERRXX(ierr);
  } else {
    *viewer = PETSC_VIEWER_STDOUT_(PETSC_COMM_WORLD);
  }
}

// -------------------------------------------------------------
// petsc_print_matrix
// -------------------------------------------------------------
static void
petsc_print_matrix(const Mat mat, const char* filename, PetscViewerFormat format)
{
  PetscErrorCode ierr;
  try {
    PetscViewer viewer;
    petsc_make_viewer(filename, &viewer);
    ierr = PetscViewerSetFormat(viewer, format); ; CHKERRXX(ierr);
    ierr = MatView(mat, viewer); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}


// -------------------------------------------------------------
// Matrix::print
// -------------------------------------------------------------
void
Matrix::print(const char* filename) const
{
  const Mat *mat(PETScMatrix(*this));
  petsc_print_matrix(*mat, filename, PETSC_VIEWER_DEFAULT);
}

// -------------------------------------------------------------
// Matrix::save
// -------------------------------------------------------------
void
Matrix::save(const char* filename) const
{
  const Mat *mat(PETScMatrix(*this));
  petsc_print_matrix(*mat, filename, PETSC_VIEWER_ASCII_MATLAB);
}


} // namespace math
} // namespace gridpack
