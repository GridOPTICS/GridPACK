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
 * @date   2014-02-20 08:32:07 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#define GA_DENSE_TRANSPOSE 0

#if GA_DENSE_TRANSPOSE
#include <ga.h>
#endif

#include <boost/assert.hpp>
#include <boost/format.hpp>
#include "matrix.hpp"
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
void 
transpose(const Matrix& A, Matrix& result)
{
  // This only works if the matrixes are the correct size
  if ((A.cols() == result.rows() &&
       A.rows() == result.cols())) {
    const Mat *pA(PETScMatrix(A));
    Mat *pAtrans(PETScMatrix(result));
    PetscErrorCode ierr(0);
    try {
      ierr = MatTranspose(*pA, MAT_REUSE_MATRIX, pAtrans); CHKERRXX(ierr);
    } catch (const PETSc::Exception& e) {
      throw PETScException(ierr, e);
    }
  } else {
    std::string msg = 
      boost::str(boost::format("transpose(Matrix, Matrix): Matrix size mismatch: (%dx%d) and (%dx%d)") %
                 A.rows() % A.cols() % result.rows() % result.cols());
    throw Exception(msg);
  }
}

// -------------------------------------------------------------
// transpose
// -------------------------------------------------------------
Matrix *
transpose(const Matrix& A)
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
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  PETScMatrixImplementation *result_impl = 
    new PETScMatrixImplementation(pAtrans, true);
  Matrix *result = new Matrix(result_impl);

  return result;
}

// -------------------------------------------------------------
// transposeMultiply
// -------------------------------------------------------------
void
transposeMultiply(const Matrix& A, const Vector& x, Vector& result)
{
  const Mat *Amat(PETScMatrix(A));
  const Vec *Xvec(PETScVector(x));
  Vec *Yvec(PETScVector(result));
  
  PetscErrorCode ierr(0);
  try {
    ierr = MatMultTranspose(*Amat, *Xvec, *Yvec); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// column
// -------------------------------------------------------------
void
column(const Matrix& A, const int& cidx, Vector& result)
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
    ierr = MatGetColumnVector(*pA, *pX, cidx); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}  

// -------------------------------------------------------------
// diagonal
// -------------------------------------------------------------
void
diagonal(const Matrix& A, Vector& result)
{
  if (result.communicator() != A.communicator()) {
    throw gridpack::Exception("incompatible: communicators do not match");
  }

  // only try this on square matrices
  if (A.rows() != A.cols()) {
    throw gridpack::Exception("can only get diagonal from square matrices");
  }

  if (result.size() != A.rows()) {
    throw gridpack::Exception("incompatible: sizes do not match");
  }

  const Mat *pA(PETScMatrix(A));
  Vec *pX(PETScVector(result));
  PetscErrorCode ierr(0);
  try {
    ierr = MatGetDiagonal(*pA, *pX); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}  

Matrix *
diagonal(const Vector& x, const Matrix::StorageType& stype)
{
  Matrix *result(new Matrix(x.communicator(), x.localSize(), x.localSize(), stype));
  const Vec *pX(PETScVector(x));
  Mat *pA(PETScMatrix(*result));
  PetscErrorCode ierr(0);
  try {
    ierr = MatDiagonalSet(*pA, *pX, INSERT_VALUES); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  return result;
}


// -------------------------------------------------------------
// multiply
// -------------------------------------------------------------
void
multiply(const Matrix& A, const Vector& x, Vector& result)
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
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

void
multiply(const Matrix& A, const Matrix& B, Matrix& result)
{
  const Mat *Amat(PETScMatrix(A));
  const Mat *Bmat(PETScMatrix(B));
  Mat *Cmat(PETScMatrix(result));

  PetscErrorCode ierr(0);
  try {
    ierr = MatDestroy(Cmat); CHKERRXX(ierr);
    ierr = MatMatMult(*Amat, *Bmat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, Cmat); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

Matrix *
multiply(const Matrix& A, const Matrix& B)
{
  const Mat *Amat(PETScMatrix(A));
  const Mat *Bmat(PETScMatrix(B));
  Mat Cmat;

  PetscErrorCode ierr(0);
  try {
    ierr = MatMatMult(*Amat, *Bmat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Cmat); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }

  PETScMatrixImplementation *result_impl = 
    new PETScMatrixImplementation(Cmat, true);
  Matrix *result = new Matrix(result_impl);
  return result;
}

// -------------------------------------------------------------
// storageType
// -------------------------------------------------------------
Matrix *
storageType(const Matrix& A, const Matrix::StorageType& new_type)
{
  Matrix *result(A.clone());
  int nproc(result->processor_size());

  MatType new_mat_type;

  if (result->storageType() != new_type) {
    switch (new_type) {
    case (Matrix::Dense):
      if (nproc > 1) {
        new_mat_type = MATDENSE;
      } else {
        new_mat_type = MATSEQDENSE;
      } 
      break;
    case (Matrix::Sparse):
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
    } catch (const PETSc::Exception& e) {
      throw PETScException(ierr, e);
    }
  }

  return result;
}

} // namespace math
} // namespace gridpack
