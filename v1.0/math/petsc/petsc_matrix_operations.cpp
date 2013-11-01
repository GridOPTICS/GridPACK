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
 * @date   2013-10-31 09:15:36 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

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
    ierr = MatTranspose(*pA, MAT_INITIAL_MATRIX, &pAtrans); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  PETScMatrixImplementation *result_impl = 
    new PETScMatrixImplementation(A.communicator(), pAtrans);
  Matrix *result = new Matrix(result_impl);

  return result;
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
  Matrix *result(new Matrix(x.communicator(), x.localSize(), x.size(), stype));
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
    new PETScMatrixImplementation(A.communicator(), Cmat);
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
