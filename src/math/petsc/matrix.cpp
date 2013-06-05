/**
 * @file   matrix.cpp
 * @author William A. Perkins
 * @date   2013-06-05 14:35:29 d3g096
 * 
 * @brief  PETSc specific part of Matrix
 * 
 * 
 */

#include "boost/assert.hpp"
#include "gridpack/math/matrix.hpp"
#include "gridpack/math/petsc/petsc_exception.hpp"
#include "gridpack/math/petsc/petsc_matrix_implementation.hpp"
#include "gridpack/math/petsc/petsc_matrix_extractor.hpp"


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
  : parallel::Distributed(comm), utility::Uncopyable(),
    p_matrix_impl()
{
  switch (storage_type) {
  case Sparse:
    p_matrix_impl.reset(new PETScMatrixImplementation(this->communicator(),
                                                      local_rows, cols, false));
    break;
  case Dense:
    p_matrix_impl.reset(new PETScMatrixImplementation(this->communicator(),
                                                      local_rows, cols, true));
    break;
  default:
    BOOST_ASSERT(false);
  }
  BOOST_ASSERT(p_matrix_impl);
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
Matrix::scale(const complex_type& xin)
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

} // namespace math
} // namespace gridpack
