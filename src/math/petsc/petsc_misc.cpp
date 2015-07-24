// -------------------------------------------------------------
// file: petsc_misc.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created May 27, 2015 by William A. Perkins
// Last Change: 2015-07-23 09:06:38 d3g096
// -------------------------------------------------------------


#include <petscsys.h>
#include <petscmath.h>
#include <petscmat.h>

#include "petsc_misc.hpp"

/// Scale a complex DENSE matrix
/** 
 * Only call this if A is real-valued but represents a complex-valued
 * matrix.  
 * 
 * @param A 
 * @param x 
 * 
 * @return 
 */
PetscErrorCode 
sillyMatScaleComplex(Mat A, const gridpack::ComplexType& px)
{
  PetscErrorCode ierr(0);
#ifndef PETSC_USE_COMPLEX
  PetscInt lo, hi;
  Mat B;

  ierr = MatDuplicate(A, MAT_SHARE_NONZERO_PATTERN, &B); CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A, &lo, &hi); CHKERRQ(ierr);
  
  for (PetscInt i = lo; i < hi; i += 2) {
    PetscInt ncols;
    const PetscInt *cols;
    const PetscScalar *r1vals;
    ierr = MatGetRow(A, i, &ncols, &cols, &r1vals); CHKERRQ(ierr);
    for (PetscInt k = 0; k < ncols; k += 2) {
      PetscInt j = cols[k];
      gridpack::ComplexType y(r1vals[k], r1vals[k+1]);
      y *= px;
      PetscInt iidx[2] = { i, i+1 };
      PetscInt jidx[2] = { j, j+1 };
      PetscScalar v[4] = { std::real(y), -std::imag(y),
                           std::imag(y), std::real(y) };
      ierr = MatSetValues(B, 2, iidx, 2, jidx, v, INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = MatRestoreRow(A, i, &ncols, &cols, &r1vals); CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatCopy(B, A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = MatDestroy(&B); CHKERRQ(ierr);
#else
  SETERRQ(PetscObjectComm((PetscObject)A), PETSC_ERR_SUP,
          "This should not have been called");
#endif
  return ierr;
}






