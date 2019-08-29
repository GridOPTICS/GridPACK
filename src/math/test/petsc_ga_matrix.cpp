// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   petsc_ga_matrix.cpp
 * @author William A. Perkins
 * @date   2019-08-01 08:51:08 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <ga.h>

#include "math.hpp"
#include "gridpack/parallel/parallel.hpp"
#include "petsc/ga_matrix.hpp"
#include "petsc/petsc_exception.hpp"

#include "test_main.cpp"

static const PetscInt local_size(5);

static
PetscErrorCode
fill_pattern(Mat A, InsertMode addv)
{
  PetscErrorCode ierr(0);
  PetscScalar x(0.0);
  PetscInt lo, hi;

  ierr = MatGetOwnershipRange(A, &lo, &hi);  CHKERRQ(ierr);
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      ierr = MatSetValues(A, 1, &i, 1, &j, &x, addv);  CHKERRQ(ierr);
      x += 1.0;
    }
  }

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);  CHKERRQ(ierr);
  return ierr;
}  

static
void
convert_and_check(Mat A)
{
  PetscErrorCode ierr(0);
  PetscInt lo, hi;
  Mat B;

  ierr = MatConvertToDenseGA(A, &B); CHKERRXX(ierr);
  ierr = MatGetOwnershipRange(B, &lo, &hi);  CHKERRXX(ierr);
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      PetscScalar x;
      PetscScalar y;
      ierr = MatGetValues(A, 1, &i, 1, &j, &x);  CHKERRXX(ierr);
      ierr = MatGetValues(B, 1, &i, 1, &j, &y);  CHKERRXX(ierr);
      BOOST_CHECK_EQUAL(x, y);
    }
  }
  ierr = MatDestroy(&B); CHKERRXX(ierr);
}  

BOOST_AUTO_TEST_SUITE(GAMatrixTest)

BOOST_AUTO_TEST_CASE( ConstructConvertDuplicate )
{
  PetscErrorCode ierr(0);
  gridpack::parallel::Communicator world;

  Mat A, B, C;
  PetscScalar x(0.0);
  PetscInt lrows, lcols;
  PetscInt lo, hi;


  BOOST_TEST_MESSAGE("Create Matrix");

  ierr = MatCreateDenseGA(world, local_size, local_size, PETSC_DETERMINE, PETSC_DETERMINE, &A); CHKERRXX(ierr);
  ierr = MatGetLocalSize(A, &lrows, &lcols);  CHKERRXX(ierr);

  BOOST_CHECK_EQUAL(lrows, 5);
  BOOST_CHECK_EQUAL(lcols, 5);

  BOOST_TEST_MESSAGE("Fill Matrix");

  ierr = fill_pattern(A, INSERT_VALUES); CHKERRXX(ierr);

  BOOST_TEST_MESSAGE("Get Matrix Values");

  ierr = MatGetOwnershipRange(A, &lo, &hi);  CHKERRXX(ierr);
  PetscScalar y(0.0);
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      ierr = MatGetValues(A, 1, &i, 1, &j, &x);  CHKERRXX(ierr);
      BOOST_CHECK_EQUAL(x, y);
      y += 1.0;
    }
  }

  BOOST_TEST_MESSAGE("Add Matrix Values");

  ierr = fill_pattern(A, ADD_VALUES); CHKERRXX(ierr);

  BOOST_TEST_MESSAGE("Get Matrix Values");

  y = 0.0;
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      ierr = MatGetValues(A, 1, &i, 1, &j, &x);  CHKERRXX(ierr);
      BOOST_CHECK_EQUAL(x, 2.0*y);
      y += 1.0;
    }
  }

  BOOST_TEST_MESSAGE("Convert Matrix");

  ierr = MatConvert(A, MATDENSE, MAT_INITIAL_MATRIX, &B);  CHKERRXX(ierr);

  BOOST_TEST_MESSAGE("View Converted Matrix");

  ierr = MatView(B, PETSC_VIEWER_STDOUT_WORLD); CHKERRXX(ierr);

  BOOST_TEST_MESSAGE("Get Matrix Values");

   for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      ierr = MatGetValues(A, 1, &i, 1, &j, &x);  CHKERRXX(ierr);
      ierr = MatGetValues(B, 1, &i, 1, &j, &y);  CHKERRXX(ierr);
      BOOST_CHECK_EQUAL(x, y);
    }
  }
  ierr = MatDestroy(&B); CHKERRXX(ierr);

  BOOST_TEST_MESSAGE("Duplicate Matrix");
  ierr = MatDuplicate(A, MAT_COPY_VALUES, &C);  CHKERRXX(ierr);

  BOOST_TEST_MESSAGE("Get Matrix Values");

  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      ierr = MatGetValues(A, 1, &i, 1, &j, &x);  CHKERRXX(ierr);
      ierr = MatGetValues(C, 1, &i, 1, &j, &y);  CHKERRXX(ierr);
      BOOST_CHECK_EQUAL(x, y);
    }
  }
  ierr = MatDestroy(&C); CHKERRXX(ierr);

  ierr = MatDestroy(&A); CHKERRXX(ierr);
}

BOOST_AUTO_TEST_CASE( Dense2GA )
{
  PetscErrorCode ierr(0);
  gridpack::parallel::Communicator comm;

  Mat A;
  PetscInt lrows = 5, lcols = 5;

  ierr = MatCreate(comm, &A); CHKERRXX(ierr);
  ierr = MatSetSizes(A, lrows, lcols, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRXX(ierr);
  ierr = MatSetType(A, MATDENSE); CHKERRXX(ierr);
  ierr = MatSetFromOptions(A); CHKERRXX(ierr);
  ierr = MatSetUp(A); CHKERRXX(ierr);
  
  ierr = fill_pattern(A, INSERT_VALUES); CHKERRXX(ierr);
  
  convert_and_check(A);

  ierr = MatDestroy(&A); CHKERRXX(ierr);
}

BOOST_AUTO_TEST_CASE( Sparse2GA )
{
  PetscErrorCode ierr(0);
  gridpack::parallel::Communicator comm;

  Mat A;
  PetscInt lrows = 5, lcols = 5;

  ierr = MatCreate(comm, &A); CHKERRXX(ierr);
  ierr = MatSetSizes(A, lrows, lcols, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRXX(ierr);
  ierr = MatSetType(A, MATAIJ); CHKERRXX(ierr);
  ierr = MatSetFromOptions(A); CHKERRXX(ierr);
  ierr = MatSetUp(A); CHKERRXX(ierr);
  
  ierr = fill_pattern(A, INSERT_VALUES); CHKERRXX(ierr);

  convert_and_check(A);

  ierr = MatDestroy(&A); CHKERRXX(ierr);
}


BOOST_AUTO_TEST_CASE( VectorMultiply )
{
  PetscErrorCode ierr(0);
  gridpack::parallel::Communicator world;

  Mat A;
  PetscInt lrows, lcols;
  PetscInt grows, gcols;
  PetscInt lo, hi;

  ierr = MatCreateDenseGA(world, local_size, local_size, PETSC_DETERMINE, PETSC_DETERMINE, &A); CHKERRXX(ierr);
  ierr = MatGetLocalSize(A, &lrows, &lcols);  CHKERRXX(ierr);
  ierr = MatGetSize(A, &grows, &gcols); CHKERRXX(ierr);

  BOOST_CHECK_EQUAL(lrows, local_size);
  BOOST_CHECK_EQUAL(lcols, local_size);

  PetscScalar v(2.0);
  ierr = MatGetOwnershipRange(A, &lo, &hi);  CHKERRXX(ierr);
  for (int i = lo; i < hi; ++i) {
    ierr = MatSetValues(A, 1, &i, 1, &i, &v, INSERT_VALUES);  CHKERRXX(ierr);
  }
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  CHKERRXX(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);  CHKERRXX(ierr);

  Vec x, y;
  ierr = VecCreate(world, &x); CHKERRXX(ierr);
  ierr = VecSetSizes(x, lcols, gcols); CHKERRXX(ierr);
  ierr = VecSetFromOptions(x); CHKERRXX(ierr);
  v = 1.0;
  ierr = VecSet(x, v); CHKERRXX(ierr);
  ierr = VecAssemblyBegin(x);
  ierr = VecAssemblyEnd(x);

  ierr = VecCreate(world, &y); CHKERRXX(ierr);
  ierr = VecSetSizes(y, lrows, grows); CHKERRXX(ierr);
  ierr = VecSetFromOptions(y); CHKERRXX(ierr);
  v = 0.0;
  ierr = VecSet(y, v); CHKERRXX(ierr);
  ierr = VecAssemblyBegin(y);
  ierr = VecAssemblyEnd(y);

  ierr = MatMult(A, x, y); CHKERRXX(ierr);

  ierr = VecGetOwnershipRange(y, &lo, &hi); CHKERRXX(ierr);
  v = 2.0;
  for (int i = lo; i < hi; ++i) {
    PetscScalar xtmp;
    ierr = VecGetValues(y, 1, &i, &xtmp); CHKERRXX(ierr);
    BOOST_CHECK_EQUAL(v, xtmp);
  }

  ierr = VecDestroy(&x); CHKERRXX(ierr);
  ierr = VecDestroy(&y); CHKERRXX(ierr);
  ierr = MatDestroy(&A); CHKERRXX(ierr);
}

BOOST_AUTO_TEST_CASE( SquareMatrixMultiply )
{
  Mat A, B, C, Cbase;
  PetscErrorCode ierr(0);
  gridpack::parallel::Communicator comm;

  // do a multiply with a pair of normal sparse matrices for a check
  // (dense multiply may not be available in parallel)

  ierr = MatCreate(comm, &A); CHKERRXX(ierr);
  ierr = MatSetSizes(A, local_size, local_size, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRXX(ierr);
  ierr = MatSetType(A, MATAIJ); CHKERRXX(ierr);
  ierr = MatSetFromOptions(A); CHKERRXX(ierr);
  ierr = MatSetUp(A); CHKERRXX(ierr);
  ierr = fill_pattern(A, INSERT_VALUES); CHKERRXX(ierr);
  ierr = MatDuplicate(A, MAT_COPY_VALUES, &B); CHKERRXX(ierr);
  ierr = MatMatMult(A, B, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Cbase); CHKERRXX(ierr);

  ierr = MatDestroy(&A); CHKERRXX(ierr);
  ierr = MatDestroy(&B); CHKERRXX(ierr);

  // Do the same multiply with GA-based matrices

  ierr = MatCreateDenseGA(comm, local_size, local_size, PETSC_DETERMINE, PETSC_DETERMINE, &A); 
  CHKERRXX(ierr);
  ierr = fill_pattern(A, INSERT_VALUES); CHKERRXX(ierr);
  ierr = MatDuplicate(A, MAT_COPY_VALUES, &B); CHKERRXX(ierr);

  ierr = MatMatMult(A, B, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);
  PetscInt lo, hi;
  ierr = MatGetOwnershipRange(Cbase, &lo, &hi); CHKERRXX(ierr);
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      PetscScalar x;
      PetscScalar y;
      ierr = MatGetValues(Cbase, 1, &i, 1, &j, &x);  CHKERRXX(ierr);
      ierr = MatGetValues(C, 1, &i, 1, &j, &y);  CHKERRXX(ierr);
      BOOST_CHECK_EQUAL(x, y);
    }
  }

  ierr = MatDestroy(&A); CHKERRXX(ierr);
  ierr = MatDestroy(&B); CHKERRXX(ierr);
  ierr = MatDestroy(&C); CHKERRXX(ierr);
  ierr = MatDestroy(&Cbase); CHKERRXX(ierr);

}

BOOST_AUTO_TEST_CASE(Transpose)
{
  PetscErrorCode ierr(0);
  gridpack::parallel::Communicator world;

  Mat A, Atrans;
  PetscInt lrows, lcols;
  PetscInt grows, gcols;
  PetscInt lo, hi;

  ierr = MatCreateDenseGA(world, local_size, local_size, PETSC_DETERMINE, PETSC_DETERMINE, &A); CHKERRXX(ierr);
  ierr = MatGetLocalSize(A, &lrows, &lcols);  CHKERRXX(ierr);
  ierr = MatGetSize(A, &grows, &gcols); CHKERRXX(ierr);

  BOOST_CHECK_EQUAL(lrows, local_size);
  BOOST_CHECK_EQUAL(lcols, local_size);

  PetscScalar v;
  ierr = MatGetOwnershipRange(A, &lo, &hi);  CHKERRXX(ierr);
  for (int i = lo; i < hi; ++i) {
    v = (PetscScalar)i;
    for (int j = lo; j < hi; ++j) {
      ierr = MatSetValues(A, 1, &i, 1, &j, &v, INSERT_VALUES);  CHKERRXX(ierr);
    }
  }
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  CHKERRXX(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);  CHKERRXX(ierr);

  ierr = MatTranspose(A, MAT_INITIAL_MATRIX, &Atrans); CHKERRXX(ierr);

  PetscScalar x;
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      v = (PetscScalar)j;
      ierr = MatGetValues(Atrans, 1, &i, 1, &j, &x);  CHKERRXX(ierr);
      BOOST_CHECK_EQUAL(v, x);
    }
  }

  ierr = MatDestroy(&Atrans); CHKERRXX(ierr);
  ierr = MatDestroy(&A); CHKERRXX(ierr);
}

BOOST_AUTO_TEST_SUITE_END()




