// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_operations.cpp
 * @author William A. Perkins
 * @date   2019-12-03 08:14:38 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#include <boost/assert.hpp>
#include <boost/format.hpp>
#include "matrix.hpp"
#include "fallback_matrix_operations.hpp"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_matrix_implementation.hpp"
#include "petsc/petsc_matrix_extractor.hpp"
#include "petsc/petsc_vector_implementation.hpp"
#include "petsc/petsc_vector_extractor.hpp"
#include "petsc/ga_matrix.hpp"

namespace gridpack {
namespace math {



// -------------------------------------------------------------
// transpose
// -------------------------------------------------------------
template <typename T, typename I>
void 
transpose(const MatrixT<T, I>& A, MatrixT<T, I>& result)
{
  // This only works if the matrixes are the correct size
  if ((A.cols() == result.rows() &&
       A.rows() == result.cols())) {
    const Mat *pA(PETScMatrix(A));
    Mat *pAtrans(PETScMatrix(result));
    PetscErrorCode ierr(0);
    try {
      ierr = MatTranspose(*pA, MAT_REUSE_MATRIX, pAtrans); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
    if (!PETScMatrixImplementation<T, I>::useLibrary) {
      result.conjugate();
    }
  } else {
    std::string msg = 
      boost::str(boost::format("transpose(Matrix, Matrix): Matrix size mismatch: (%dx%d) and (%dx%d)") %
                 A.rows() % A.cols() % result.rows() % result.cols());
    throw Exception(msg);
  }
}

template
void 
transpose<ComplexType, int>(const MatrixT<ComplexType, int>& A, 
                            MatrixT<ComplexType, int>& result);

template
void 
transpose<RealType, int>(const MatrixT<RealType, int>& A, 
                         MatrixT<RealType, int>& result);

// -------------------------------------------------------------
// transpose
// -------------------------------------------------------------
template <typename T, typename I>
MatrixT<T, I> *
transpose(const MatrixT<T, I>& A)
{
  PetscErrorCode ierr(0);
  const Mat *pA(PETScMatrix(A));
  MatrixT<T, I> *result;

  try {

    // There is some indication that dense complex transposes using a
    // real PETSc library are distributed incorrectly, perhaps
    // splitting the rows between the rows representing a single
    // complex row.  So, let's be explicit about the distribution of
    // the transpose matrix before we get PETSc to compute it.

    // Doing this is harder with sparse matrix, so let's just handle
    // the dense case (for which we have a actual problem.

    // if (!PETScMatrixImplementation<T, I>::useLibrary &&
    //     A.storageType() == Dense) {
      
    //   PetscInt N;
    //   PetscInt n;
      
    //   N = A.cols();
    //   n = PETSC_DECIDE;
    //   ierr = PetscSplitOwnership(A.communicator(), &n, &N);
    //   I trows(n);
      
    //   N = A.rows();
    //   n = PETSC_DECIDE;
    //   ierr = PetscSplitOwnership(A.communicator(), &n, &N);
    //   I tcols(n);

    //   result = new MatrixT<T, I>(A.communicator(), trows, tcols, A.storageType());

    //   transpose(A, *result);

    //   return result;

    // } else {

      Mat pAtrans;
      ierr = MatTranspose(*pA, MAT_INITIAL_MATRIX, &pAtrans); CHKERRXX(ierr);
      
      PETScMatrixImplementation<T, I> *result_impl = 
        new PETScMatrixImplementation<T, I>(pAtrans, false, true);
      result = new MatrixT<T, I>(result_impl);
    // }
  
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
  if (!PETScMatrixImplementation<T, I>::useLibrary) {
    result->conjugate();
  }
  return result;
}

template
MatrixT<ComplexType, int> *
transpose(const MatrixT<ComplexType, int>& A);

template
MatrixT<RealType, int> *
transpose(const MatrixT<RealType, int>& A);

// -------------------------------------------------------------
// transposeMultiply
// -------------------------------------------------------------
/** 
 * Works for complex regardless of underlying PETSc element type.
 * 
 * @param T 
 * @param A 
 * @param T 
 * @param x 
 * @param T 
 * @param result 
 */
template <typename T, typename I>
void
transposeMultiply(const MatrixT<T, I>& A, const VectorT<T, I>& x, VectorT<T, I>& result)
{
  const Mat *Amat(PETScMatrix(A));
  const Vec *Xvec(PETScVector(x));
  Vec *Yvec(PETScVector(result));
  
  PetscErrorCode ierr(0);
  try {
    ierr = MatMultTranspose(*Amat, *Xvec, *Yvec); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

template
void
transposeMultiply(const MatrixT<ComplexType, int>& A, 
                  const VectorT<ComplexType, int>& x, 
                  VectorT<ComplexType, int>& result);

template
void
transposeMultiply(const MatrixT<RealType, int>& A, 
                  const VectorT<RealType, int>& x, 
                  VectorT<RealType, int>& result);

// -------------------------------------------------------------
// column
// -------------------------------------------------------------
/** 
 * Works for complex regardless of underlying PETSc element type
 * 
 * @param T 
 * @param A 
 * @param cidx 
 * @param T 
 * @param result 
 */
template <typename T, typename I>
void
column(const MatrixT<T, I>& A, const int& cidx, VectorT<T, I>& result)
{
  if (result.communicator() != A.communicator()) {
    throw gridpack::Exception("incompatible: communicators do not match");
  }

  // this is a requirement of PETSc
  if (result.localSize() != A.localRows()) {
    throw gridpack::Exception("incompatible: sizes do not match");
  }

  const Mat *pA(PETScMatrix(A));
  Vec *pX(PETScVector(result));
  PetscErrorCode ierr(0);
  try {
    PetscInt c(cidx*PetscElementSize<T>::value);
    ierr = MatGetColumnVector(*pA, *pX, c); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}  

template
void
column(const MatrixT<ComplexType, int>& A, 
       const int& cidx, VectorT<ComplexType, int>& result);

template
void
column(const MatrixT<RealType, int>& A, 
       const int& cidx, VectorT<RealType, int>& result);

// -------------------------------------------------------------
// diagonal
// -------------------------------------------------------------
/** 
 * For complex matrix/vector on top of a real PETSc, a fallback has to
 * be used.
 * 
 * @param A 
 * @param result 
 */
template <typename T, typename I>
void
diagonal(const MatrixT<T, I>& A, VectorT<T, I>& result)
{
  if (result.communicator() != A.communicator()) {
    throw gridpack::Exception("incompatible: communicators do not match");
  }

  // only try this on square matrices
  if (A.rows() != A.cols()) {
    throw gridpack::Exception("can only get diagonal from square matrices");
  }

  if (result.size() != A.rows()) {
    char buf[128];
    sprintf(buf,"Matrix::diagonal incompatible: sizes do not match."
            " Matrix rows: %d Vector length: %d",A.rows(),result.size());
    throw gridpack::Exception(buf);
  }

  if (PETScMatrixImplementation<T, I>::useLibrary) {
    const Mat *pA(PETScMatrix(A));
    Vec *pX(PETScVector(result));
    PetscErrorCode ierr(0);
    try {
      ierr = MatGetDiagonal(*pA, *pX); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  } else {
    fallback::diagonal(A, result);
  }
}  

template
void
diagonal<ComplexType, int>(const MatrixT<ComplexType, int>& A, 
                           VectorT<ComplexType, int>& result);

template
void
diagonal<RealType, int>(const MatrixT<RealType, int>& A, 
                        VectorT<RealType, int>& result);

template <typename T, typename I>
MatrixT<T, I> *
diagonal(const VectorT<T, I>& x, const MatrixStorageType& stype)
{
  MatrixT<T, I> *result;

  if (PETScMatrixImplementation<T, I>::useLibrary) {
    result = new MatrixT<T,I>(x.communicator(), 
                              x.localSize(), x.localSize(), stype);
    const Vec *pX(PETScVector(x));
    Mat *pA(PETScMatrix(*result));
    PetscErrorCode ierr(0);
    try {
      ierr = MatDiagonalSet(*pA, *pX, INSERT_VALUES); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  } else {
    result = fallback::diagonal(x, stype);
  }
  return result;
}

template
MatrixT<ComplexType, int> *
diagonal(const VectorT<ComplexType, int>& x, const MatrixStorageType& stype);

template
MatrixT<RealType, int> *
diagonal(const VectorT<RealType, int>& x, const MatrixStorageType& stype);


// -------------------------------------------------------------
// (Vector) multiply
// -------------------------------------------------------------
template <typename T, typename I>
void
multiply(const MatrixT<T, I>& A, const VectorT<T, I>& x, VectorT<T, I>& result)
{
  const Mat *Amat(PETScMatrix(A));
  const Vec *Xvec(PETScVector(x));
  Vec *Yvec(PETScVector(result));

  PetscErrorCode ierr(0);
  try {
    PetscInt plo, phi;
    MatType type;
    ierr = MatGetType(*Amat,&type);
    ierr = MatGetOwnershipRange(*Amat,&plo,&phi);
    ierr = MatGetSize(*Amat,&plo,&phi);
    ierr = VecGetOwnershipRange(*Xvec,&plo,&phi);
    ierr = VecGetSize(*Xvec,&plo);
    ierr = VecGetOwnershipRange(*Yvec,&plo,&phi);
    ierr = VecGetSize(*Yvec,&plo);
    ierr = MatMult(*Amat, *Xvec, *Yvec); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

template
void
multiply(const MatrixT<ComplexType, int>& A, 
         const VectorT<ComplexType, int>& x, 
         VectorT<ComplexType, int>& result);

template
void
multiply(const MatrixT<RealType, int>& A, 
         const VectorT<RealType, int>& x, 
         VectorT<RealType, int>& result);


// -------------------------------------------------------------
// multiply_dense
// -------------------------------------------------------------
static 
PetscErrorCode 
multiply_dense(const Mat& A, const Mat& B, Mat& C)
{
  PetscErrorCode ierr(0);
  Mat Aga, Bga, Cga;

  ierr = MatConvertToDenseGA(A, &Aga); CHKERRQ(ierr);
  ierr = MatConvertToDenseGA(B, &Bga); CHKERRQ(ierr);
  ierr = MatMatMult(Aga, Bga, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Cga); CHKERRQ(ierr);
  ierr = MatConvertGAToDense(Cga, &C); CHKERRXX(ierr);
  
  ierr = MatDestroy(&Aga); CHKERRQ(ierr);
  ierr = MatDestroy(&Bga); CHKERRQ(ierr);
  ierr = MatDestroy(&Cga); CHKERRQ(ierr);
  return ierr;
}

// -------------------------------------------------------------
// (Matrix) multiply
// -------------------------------------------------------------
template <typename T, typename I>
void
multiply(const MatrixT<T, I>& A, const MatrixT<T, I>& B, MatrixT<T, I>& result)
{
  PetscErrorCode ierr(0);

  // special method required for parallel dense*dense
  if (A.communicator().size() > 1 &&
      A.storageType() == Dense &&
      B.storageType() == Dense) {
    const Mat *Amat(PETScMatrix(A));
    const Mat *Bmat(PETScMatrix(B));
    Mat *Cmat(PETScMatrix(result));
    ierr = MatMultbyGA(*Amat, *Bmat, *Cmat); CHKERRXX(ierr);
    // ierr = multiply_dense(*Amat, *Bmat, *Cmat); CHKERRXX(ierr);
  } else {
    const Mat *Amat(PETScMatrix(A));
    const Mat *Bmat(PETScMatrix(B));
    Mat *Cmat(PETScMatrix(result));
    
    try {
      ierr = MatDestroy(Cmat); CHKERRXX(ierr);
      ierr = MatMatMult(*Amat, *Bmat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, Cmat); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
  }
}

template
void
multiply(const MatrixT<ComplexType, int>& A, 
         const MatrixT<ComplexType, int>& B, 
         MatrixT<ComplexType, int>& result);

template
void
multiply(const MatrixT<RealType, int>& A, 
         const MatrixT<RealType, int>& B, 
         MatrixT<RealType, int>& result);

template <typename T, typename I>
MatrixT<T, I> *
multiply(const MatrixT<T, I>& A, const MatrixT<T, I>& B)
{
  PetscErrorCode ierr(0);
  MatrixT<T, I> *result;
  // special method required for parallel dense*dense
  if (A.communicator().size() > 1 &&
      A.storageType() == Dense &&
      B.storageType() == Dense) {
    const Mat *Amat(PETScMatrix(A));
    const Mat *Bmat(PETScMatrix(B));
    Mat Cmat;
    ierr = MatMultbyGA_new(*Amat, *Bmat, Cmat);

    PETScMatrixImplementation<T, I> *result_impl = 
      new PETScMatrixImplementation<T, I>(Cmat, false, true);
    result = new MatrixT<T, I>(result_impl);
  } else {
    const Mat *Amat(PETScMatrix(A));
    const Mat *Bmat(PETScMatrix(B));
    Mat Cmat;

    try {
      ierr = MatMatMult(*Amat, *Bmat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Cmat); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }

    PETScMatrixImplementation<T, I> *result_impl = 
      new PETScMatrixImplementation<T, I>(Cmat, false, true);
    result = new MatrixT<T, I>(result_impl);
  }
  return result;
}

template
MatrixT<ComplexType, int> *
multiply(const MatrixT<ComplexType, int>& A, const MatrixT<ComplexType, int>& B);

template
MatrixT<RealType, int> *
multiply(const MatrixT<RealType, int>& A, const MatrixT<RealType, int>& B);

// -------------------------------------------------------------
// storageType
// -------------------------------------------------------------
template <typename T, typename I>
MatrixStorageType
MatrixT<T, I>::storageType(void) const
{
  MatrixStorageType result;
  PetscErrorCode ierr;
  const Mat *mat(PETScMatrix(*this));
  try {
    MatType thetype;
    ierr = MatGetType(*mat, &thetype); CHKERRXX(ierr);

    // this relies on the fact that MatType is actual a pointer to a
    // char string -- need to test string contents, not pointers
    std::string stype(thetype);

    if (stype == MATSEQDENSE ||
        stype == MATDENSE ||
        stype == MATMPIDENSE) {
      result = Dense;
    } else if (stype == MATAIJ || 
               stype == MATSEQAIJ || 
               stype == MATMPIAIJ) {
      result = Sparse;
    } else {
      std::string msg("Matrix: unexpected PETSc storage type: ");
      msg += "\"";
      msg += stype;
      msg += "\"";
      throw Exception(msg);
    }
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  } catch (const Exception& e) {
    throw;
  }
  return result;
}

template 
MatrixStorageType
MatrixT<ComplexType>::storageType(void) const;

template 
MatrixStorageType
MatrixT<RealType>::storageType(void) const;

template <typename T, typename I>
MatrixT<T, I> *
storageType(const MatrixT<T, I>& A, const MatrixStorageType& new_type)
{
  int nproc(A.processor_size());

  MatrixT<T, I> *result;
  MatType new_mat_type(MATSEQAIJ);

  if (A.storageType() != new_type) {
    switch (new_type) {
    case (Dense):
      if (nproc > 1) {
        new_mat_type = MATDENSE;
      } else {
        new_mat_type = MATSEQDENSE;
      } 
      break;
    case (Sparse):
      if (nproc > 1) {
        new_mat_type = MATMPIAIJ;
      } else {
        new_mat_type = MATSEQAIJ;
      } 
      break;
    }
  
    const Mat *Amat(PETScMatrix(A));
    Mat B;
    PetscErrorCode ierr(0);
    try {
      ierr = MatConvert(*Amat, new_mat_type, MAT_INITIAL_MATRIX, &B); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }

    PETScMatrixImplementation<T, I> *result_impl = 
      new PETScMatrixImplementation<T, I>(B, true, true);
    result = new MatrixT<T, I>(result_impl);
    ierr = MatDestroy(&B); CHKERRXX(ierr);
    
  } else {
    result = A.clone();
  }

  return result;
}

template
MatrixT<ComplexType, int> *
storageType(const MatrixT<ComplexType, int>& A, 
            const MatrixStorageType& new_type);

template
MatrixT<RealType, int> *
storageType(const MatrixT<RealType, int>& A, 
            const MatrixStorageType& new_type);

// -------------------------------------------------------------
// matrixLoadBinary
// -------------------------------------------------------------
template <typename T, typename I>
MatrixT<T, I> *
matrixLoadBinary(const parallel::Communicator& comm, const char* filename)
{
  PetscErrorCode ierr;
  MatrixT<T, I> *result = NULL;
  try {
    Mat mat;
    ierr = MatCreate(comm, &mat); CHKERRXX(ierr);
    if (comm.size() == 1) {
      ierr = MatSetType(mat, MATSEQAIJ); CHKERRXX(ierr);
    } else {
      ierr = MatSetType(mat, MATMPIAIJ); CHKERRXX(ierr);
    }
    ierr = MatSetFromOptions(mat); CHKERRXX(ierr);

    PetscViewer viewer;
    ierr = PetscViewerBinaryOpen(comm,
                                 filename,
                                 FILE_MODE_READ,
                                 &viewer); CHKERRXX(ierr);
#if PETSC_VERSION_GE(3,7,0)
    ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_NATIVE); CHKERRXX(ierr);
#else
    ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_NATIVE); CHKERRXX(ierr);
#endif
    ierr = MatLoad(mat, viewer); CHKERRXX(ierr);

#if PETSC_VERSION_GE(3,7,0)
    ierr = PetscViewerPopFormat(viewer); CHKERRXX(ierr);
#endif
    ierr = PetscViewerDestroy(&viewer); CHKERRXX(ierr);

    PETScMatrixImplementation<T, I> *result_impl = 
      new PETScMatrixImplementation<T, I>(mat, true);
    result = new MatrixT<T, I>(result_impl);
    ierr = MatDestroy(&mat); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
  return result;
}  

template
MatrixT<ComplexType, int> *
matrixLoadBinary(const parallel::Communicator& comm, const char* filename);

template
MatrixT<RealType, int> *
matrixLoadBinary(const parallel::Communicator& comm, const char* filename);


} // namespace math
} // namespace gridpack
