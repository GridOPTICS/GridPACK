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
 * @date   2013-10-09 13:24:58 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include <boost/assert.hpp>
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
  result.equate(A);
  Mat *pA(PETScMatrix(result));
  PetscErrorCode ierr(0);
  try {
    PetscScalar one(1.0);
    ierr = MatTranspose(*pA, MAT_REUSE_MATRIX, pA); CHKERRXX(ierr);
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

  

} // namespace math
} // namespace gridpack
