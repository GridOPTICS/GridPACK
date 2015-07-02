// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   densemat.cpp
 * @author William A. Perkins
 * @date   2015-06-04 13:47:41 d3g096
 * 
 * @brief  A test of dense PETSc matrix
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created June  3, 2015 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#include <iostream>
#include <cassert>
#include <boost/mpi.hpp>
#include <petsc.h>

// -------------------------------------------------------------
// fillIt
// -------------------------------------------------------------
PetscErrorCode
fillIt(Mat A)
{
  PetscErrorCode ierr(0);
  PetscInt lo, hi;
  ierr = MatGetOwnershipRange(A, &lo, &hi); CHKERRQ(ierr);

  PetscScalar v(0.0);
  for (PetscInt i = lo; i < hi; i += 1) {
    for (PetscInt j = lo; j < hi; j += 1) {
      ierr = MatSetValue(A, i, j, v, INSERT_VALUES);
      v += 1.0;
    }
  }
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return ierr;
}

// -------------------------------------------------------------
// getRows
// -------------------------------------------------------------
PetscErrorCode
getRows(Mat A)
{
  PetscErrorCode ierr(0);
  PetscInt lo, hi;

  ierr = MatGetOwnershipRange(A, &lo, &hi); CHKERRQ(ierr);

  for (PetscInt i = lo; i < hi; ++i) {
    PetscInt ncols;
    const PetscScalar *p;
    ierr = MatGetRow(A, i, &ncols, PETSC_NULL, &p); CHKERRQ(ierr);
    ierr = MatRestoreRow(A, i, &ncols, PETSC_NULL, &p); CHKERRQ(ierr);
  }
  return ierr;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  int rank(world.rank()), size(world.size());

  PetscErrorCode ierr(0);
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);

  PetscInt lrow(5), lcol(5);
  Mat A, B;

  ierr = MatCreate(world, &A); CHKERRXX(ierr);
  ierr = MatSetSizes(A, lrow, lcol, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRXX(ierr);
  ierr = MatSetType(A, MATDENSE); CHKERRXX(ierr);
  ierr = MatSetFromOptions(A); CHKERRXX(ierr);
  ierr = MatSetUp(A); CHKERRXX(ierr);
  ierr = fillIt(A); CHKERRXX(ierr);

  ierr = MatConvert(A, MATAIJ, MAT_INITIAL_MATRIX, &B); CHKERRXX(ierr);

  // ierr = MatCreate(world, &B); CHKERRXX(ierr);
  // ierr = MatSetSizes(B, lrow, lcol, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRXX(ierr);
  // ierr = MatSetType(B, MATAIJ); CHKERRXX(ierr);
  // ierr = MatSetFromOptions(B); CHKERRXX(ierr);
  // ierr = MatSetUp(B); CHKERRXX(ierr);

  ierr = fillIt(B); CHKERRXX(ierr);

  getRows(A);
  getRows(B);

  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = MatDestroy(&B); CHKERRQ(ierr);
  return 0;
}

