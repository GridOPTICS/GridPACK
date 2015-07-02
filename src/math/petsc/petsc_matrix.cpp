// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   matrix.cpp
 * @author William A. Perkins
 * @date   2015-06-16 08:34:11 d3g096
 * 
 * @brief  PETSc specific part of Matrix
 * 
 * 
 */
// -------------------------------------------------------------

#include "boost/assert.hpp"
#include "boost/format.hpp"
#include "matrix.hpp"
#include "fallback_matrix_operations.hpp"
#include "complex_operators.hpp"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_matrix_implementation.hpp"
#include "petsc/petsc_matrix_extractor.hpp"
#include "petsc/petsc_vector_extractor.hpp"


namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class Matrix
// -------------------------------------------------------------

// -------------------------------------------------------------
// Matrix:: constructors / destructor
// -------------------------------------------------------------
template <typename T, typename I>
MatrixT<T, I>::MatrixT(const parallel::Communicator& comm,
                       const int& local_rows,
                       const int& cols,
                       const MatrixStorageType& storage_type)
  : parallel::WrappedDistributed(), utility::Uncopyable(),
    p_matrix_impl()
{
  switch (storage_type) {
  case Sparse:
    p_matrix_impl.reset(new PETScMatrixImplementation<T, I>(comm,
                                                            local_rows, cols, false));
    break;
  case Dense:
    p_matrix_impl.reset(new PETScMatrixImplementation<T, I>(comm,
                                                            local_rows, cols, true));
    break;
  default:
    BOOST_ASSERT(false);
  }
  BOOST_ASSERT(p_matrix_impl);
  p_setDistributed(p_matrix_impl.get());
}

template 
MatrixT<ComplexType>::MatrixT(const parallel::Communicator& comm,
                              const int& local_rows,
                              const int& cols,
                              const MatrixStorageType& storage_type);

template 
MatrixT<RealType>::MatrixT(const parallel::Communicator& comm,
                           const int& local_rows,
                           const int& cols,
                           const MatrixStorageType& storage_type);

template <typename T, typename I>
MatrixT<T, I>::MatrixT(const parallel::Communicator& comm,
                       const int& local_rows,
                       const int& cols,
                       const int& max_nz_per_row)
  : parallel::WrappedDistributed(), utility::Uncopyable(),
    p_matrix_impl()
{
  p_matrix_impl.reset(new PETScMatrixImplementation<T, I>(comm,
                                                          local_rows, cols, 
                                                          max_nz_per_row));
  BOOST_ASSERT(p_matrix_impl);
  p_setDistributed(p_matrix_impl.get());
}

template 
MatrixT<ComplexType>::MatrixT(const parallel::Communicator& comm,
                              const int& local_rows,
                              const int& cols,
                              const int& max_nz_per_row);

template 
MatrixT<RealType>::MatrixT(const parallel::Communicator& comm,
                           const int& local_rows,
                           const int& cols,
                           const int& max_nz_per_row);

template <typename T, typename I>
MatrixT<T, I>::MatrixT(const parallel::Communicator& comm,
                       const int& local_rows,
                       const int& cols,
                       const int *nz_by_row)
  : parallel::WrappedDistributed(), utility::Uncopyable(),
    p_matrix_impl()
{
  p_matrix_impl.reset(new PETScMatrixImplementation<T, I>(comm,
                                                          local_rows, cols, 
                                                          nz_by_row));
  BOOST_ASSERT(p_matrix_impl);
  p_setDistributed(p_matrix_impl.get());
}

template 
MatrixT<ComplexType>::MatrixT(const parallel::Communicator& comm,
                              const int& local_rows,
                              const int& cols,
                              const int *nz_by_row);

template 
MatrixT<RealType>::MatrixT(const parallel::Communicator& comm,
                              const int& local_rows,
                              const int& cols,
                              const int *nz_by_row);


// -------------------------------------------------------------
// Matrix::createDense
// -------------------------------------------------------------
template <typename T, typename I>
MatrixT<T, I> *
MatrixT<T, I>::createDense(const parallel::Communicator& comm,
                           const int& global_rows,
                           const int& global_cols,
                           const int& local_rows,
                           const int& local_cols)
{
  PETScMatrixImplementation<T, I> *impl = 
    PETScMatrixImplementation<T, I>::createDense(comm, 
                                                 global_rows, global_cols,
                                                 local_rows, local_cols);
  MatrixT<T, I> *result = new MatrixT<T, I>(impl);
  return result;
}

template 
MatrixT<ComplexType> *
MatrixT<ComplexType>::createDense(const parallel::Communicator& comm,
                                  const int& global_rows,
                                  const int& global_cols,
                                  const int& local_rows,
                                  const int& local_cols);

template 
MatrixT<RealType> *
MatrixT<RealType>::createDense(const parallel::Communicator& comm,
                               const int& global_rows,
                               const int& global_cols,
                               const int& local_rows,
                               const int& local_cols);


// -------------------------------------------------------------
// Matrix::equate
// -------------------------------------------------------------
template <typename T, typename I>
void
MatrixT<T, I>::equate(const MatrixT<T, I>& B)
{
  this->p_check_compatible(B);
  Mat *pA(PETScMatrix(*this));
  const Mat *pB(PETScMatrix(B));

  PetscErrorCode ierr(0);
  try {
    ierr = MatCopy(*pB, *pA, DIFFERENT_NONZERO_PATTERN); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

template 
void
MatrixT<ComplexType>::equate(const MatrixT<ComplexType>& B);

template 
void
MatrixT<RealType>::equate(const MatrixT<RealType>& B);

// -------------------------------------------------------------
// Matrix::add
// -------------------------------------------------------------
template <typename T, typename I>
void
MatrixT<T, I>::add(const MatrixT<T, I>& B)
{
  this->p_check_compatible(B);
  Mat *pA(PETScMatrix(*this));
  const Mat *pB(PETScMatrix(B));

  PetscErrorCode ierr(0);
  try {
    PetscScalar one(1.0);
    ierr = MatAXPY(*pA, one, *pB, DIFFERENT_NONZERO_PATTERN); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

template 
void
MatrixT<ComplexType>::add(const MatrixT<ComplexType>& B);

template 
void
MatrixT<RealType>::add(const MatrixT<RealType>& B);

// -------------------------------------------------------------
// Matrix::addDiagonal
// -------------------------------------------------------------
template <typename T, typename I>
void
MatrixT<T, I>::addDiagonalVector(const VectorT<T, I>& x)
{
  if (PETScMatrixImplementation<T, I>::useLibrary) {
    const Vec *pX(PETScVector(x));
    Mat *pA(PETScMatrix(*this));
    PetscErrorCode ierr(0);
    try {
      ierr = MatDiagonalSet(*pA, *pX, ADD_VALUES); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  } else {
    fallback::addDiagonal<T, I>(*this, x);
  }
}

template 
void
MatrixT<ComplexType>::addDiagonalVector(const VectorT<ComplexType>& x);

template 
void
MatrixT<RealType>::addDiagonalVector(const VectorT<RealType>& x);

// -------------------------------------------------------------
// Matrix::multiplyDiagonal
// -------------------------------------------------------------
template <typename T, typename I>
void
MatrixT<T, I>::multiplyDiagonal(const VectorT<T, I>& x)
{
  if (PETScMatrixImplementation<T, I>::useLibrary) {
    const Vec *pscale(PETScVector(x));
    Mat *pA(PETScMatrix(*this));
    PetscErrorCode ierr(0);
    try {
      Vec diagorig, diagnew;
      ierr = VecDuplicate(*pscale, &diagorig);  CHKERRXX(ierr);
      ierr = VecDuplicate(*pscale, &diagnew);  CHKERRXX(ierr);
      ierr = MatGetDiagonal(*pA, diagorig); CHKERRXX(ierr);
      ierr = VecPointwiseMult(diagnew, diagorig, *pscale); CHKERRXX(ierr);
      ierr = MatDiagonalSet(*pA, diagnew, INSERT_VALUES); CHKERRXX(ierr);
      ierr = VecDestroy(&diagorig); CHKERRXX(ierr); 
      ierr = VecDestroy(&diagnew); CHKERRXX(ierr); 
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  } else {
    fallback::multiplyDiagonal<T, I>(*this, x);
  }
}

template 
void
MatrixT<ComplexType>::multiplyDiagonal(const VectorT<ComplexType>& x);

template 
void
MatrixT<RealType>::multiplyDiagonal(const VectorT<RealType>& x);

// -------------------------------------------------------------
// Matrix::storageType
// -------------------------------------------------------------
} // namespace math
} // namespace gridpack
