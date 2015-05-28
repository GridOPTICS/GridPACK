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
 * @date   2015-05-28 10:13:16 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#define GA_DENSE_TRANSPOSE 0

#if GA_DENSE_TRANSPOSE
#include <ga.h>
#endif

#include <boost/assert.hpp>
#include <boost/format.hpp>
#include "matrix.hpp"
#include "fallback_matrix_operations.hpp"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_matrix_implementation.hpp"
#include "petsc/petsc_matrix_extractor.hpp"
#include "petsc/petsc_vector_implementation.hpp"
#include "petsc/petsc_vector_extractor.hpp"

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
  const Mat *pA(PETScMatrix(A));
  Mat pAtrans;
  
  PetscErrorCode ierr(0);
  try {
#if GA_DENSE_TRANSPOSE
    MatType type;
    MatGetType(*pA,&type);
    if (strcmp(type,MATMPIDENSE) == 0 && A.processor_size() > 1) {
      // This is a hack to get around problems with changing distributions in
      // PETSc. Copy the matrix to a global array, do the transpose in GA, and
      // then copy to a new PETSc matrix while preserving the distribution of the
      // original PETSc matrix.
      // Get dimensions of matrix A and allocate local arrays
      int lrows(A.localRows());
      int grows(A.rows());
      int gcols(A.cols());
      ComplexType vals[lrows*gcols];
      int nprocs = A.processor_size();
      int me = A.processor_rank();
      int tmapc[nprocs+1];
      int mapc[nprocs+1];
      //printf("p[%d] (transpose) Got to 1\n",me);
      // Set up global arrays with same row partition as A
      int i, j;
      for (i=0; i<nprocs+1; i++) tmapc[i] = 0;
      tmapc[me] = lrows;
      GA_Pgroup_igop(A.communicator().getGroup(),tmapc,nprocs+1,"+");
      mapc[0] = 0;
      for (i=1; i<nprocs; i++) mapc[i] = mapc[i-1]+tmapc[i-1];
      mapc[nprocs] = 0;
      int dims[2],blocks[2];
      dims[0] = grows;
      dims[1] = gcols;
      blocks[0] = nprocs;
      blocks[1] = 1;
      //printf("p[%d] (transpose) Got to 2\n",me);
      int g_tmp = GA_Create_handle();
      GA_Set_data(g_tmp,2,dims,C_DCPL);
      GA_Set_irreg_distr(g_tmp,mapc,blocks);
      GA_Set_pgroup(g_tmp,A.communicator().getGroup());
      if (!GA_Allocate(g_tmp)) {
        //TODO: some kind of error
      }
      int g_trns = GA_Create_handle();
      GA_Set_data(g_trns,2,dims,C_DCPL);
      GA_Set_irreg_distr(g_trns,mapc,blocks);
      GA_Set_pgroup(g_trns,A.communicator().getGroup());
      if (!GA_Allocate(g_trns)) {
        //TODO: some kind of error
      }
      //printf("p[%d] (transpose) Got to 3\n",me);
      // Get low and high row indices for matrix A and copy local row block of A
      // to g_tmp
      int llo,lhi, ncount;
      A.localRowRange(llo,lhi);
      ncount = 0;
      for (i=llo; i<lhi; i++) {
        for (j=0; j<gcols; j++) {
          A.getElement(i,j,vals[ncount]);
          ncount++;
        }
      }
      int lo[2],hi[2];
      lo[0] = llo;
      lo[1] = 0;
      hi[0] = lhi-1;
      hi[1] = gcols-1;
      int ld = gcols;
      //printf("p[%d] (transpose) Got to 4\n",me);
      // Transpose g_tmp
      NGA_Put(g_tmp,lo,hi,vals,&ld);
      GA_Transpose(g_tmp, g_trns);
      NGA_Get(g_trns,lo,hi,vals,&ld);
      // Copy g_trns to new PETSc matrix
      GA_Pgroup_sync(A.communicator().getGroup());
      //printf("p[%d] (transpose) Got to 5\n",me);
      ierr = MatCreate(A.communicator(), &pAtrans); CHKERRXX(ierr);
      ierr = MatSetSizes(pAtrans, lrows,lrows,PETSC_DETERMINE,PETSC_DETERMINE); CHKERRXX(ierr);
      //printf("p[%d] (transpose) Got to 5a\n",me);
      ierr = MatSetType(pAtrans, type); CHKERRXX(ierr);
      //printf("p[%d] (transpose) Got to 5b\n",me);
      ierr = MatMPIDenseSetPreallocation(pAtrans, PETSC_NULL); CHKERRXX(ierr);
      //printf("p[%d] (transpose) Got to 5c\n",me);
      ncount = 0;
      for (i=llo; i<lhi; i++) {
        for (j=0; j<gcols; j++) {
          ierr = MatSetValue(pAtrans,i,j,vals[ncount],INSERT_VALUES); CHKERRXX(ierr);
          ncount++;
        }
      }
      //printf("p[%d] (transpose) Got to 5d\n",me);
      ierr = MatAssemblyBegin(pAtrans, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
      ierr = MatAssemblyEnd(pAtrans, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
      //printf("p[%d] (transpose) Got to 6\n",me);
      GA_Destroy(g_tmp);
      GA_Destroy(g_trns);
      GA_Pgroup_sync(A.communicator().getGroup());
      //printf("p[%d] (transpose) Got to 7\n",me);
    } else {
      ierr = MatTranspose(*pA, MAT_INITIAL_MATRIX, &pAtrans); CHKERRXX(ierr);
    }
#else
    ierr = MatTranspose(*pA, MAT_INITIAL_MATRIX, &pAtrans); CHKERRXX(ierr);
#endif
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
  PETScMatrixImplementation<T, I> *result_impl = 
    new PETScMatrixImplementation<T, I>(pAtrans, true);
  MatrixT<T, I> *result = new MatrixT<T, I>(result_impl);

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
// (Matrix) multiply
// -------------------------------------------------------------
template <typename T, typename I>
void
multiply(const MatrixT<T, I>& A, const MatrixT<T, I>& B, MatrixT<T, I>& result)
{
  // special method required for parallel dense*dense
  if (A.communicator().size() > 1 &&
      A.storageType() == Dense &&
      B.storageType() == Dense) {
    fallback::denseMatrixMultiply(A, B, result);
  } else {
    const Mat *Amat(PETScMatrix(A));
    const Mat *Bmat(PETScMatrix(B));
    Mat *Cmat(PETScMatrix(result));
    
    PetscErrorCode ierr(0);
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
  MatrixT<T, I> *result;
  // special method required for parallel dense*dense
  if (A.communicator().size() > 1 &&
      A.storageType() == Dense &&
      B.storageType() == Dense) {
    result = fallback::denseMatrixMultiply(A, B);
  } else {
    const Mat *Amat(PETScMatrix(A));
    const Mat *Bmat(PETScMatrix(B));
    Mat Cmat;

    PetscErrorCode ierr(0);
    try {
      ierr = MatMatMult(*Amat, *Bmat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Cmat); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }

    PETScMatrixImplementation<T, I> *result_impl = 
      new PETScMatrixImplementation<T, I>(Cmat, true);
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
  MatrixT<T, I> *result(A.clone());
  int nproc(result->processor_size());

  MatType new_mat_type(MATSEQAIJ);

  if (result->storageType() != new_type) {
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
  
    Mat *mat(PETScMatrix(*result));
    PetscErrorCode ierr(0);
    try {
      ierr = MatConvert(*mat, new_mat_type, MAT_REUSE_MATRIX, mat); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
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
    ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_NATIVE); CHKERRXX(ierr);
    ierr = MatLoad(mat, viewer); CHKERRXX(ierr);
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
