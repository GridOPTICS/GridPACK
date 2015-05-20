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
 * @date   2015-05-20 06:58:36 d3g096
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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

static const PetscInt local_size(5);

BOOST_AUTO_TEST_SUITE(GAMatrixTest)

BOOST_AUTO_TEST_CASE( ConstructConvertDuplicate )
{
  PetscErrorCode ierr(0);
  gridpack::parallel::Communicator world;

  Mat A, B, C;
  PetscInt lrows, lcols;
  PetscInt lo, hi;


  BOOST_TEST_MESSAGE("Create Matrix");

  ierr = MatCreateDenseGA(world, local_size, local_size, PETSC_DETERMINE, PETSC_DETERMINE, &A); CHKERRXX(ierr);
  ierr = MatGetLocalSize(A, &lrows, &lcols);  CHKERRXX(ierr);

  BOOST_CHECK_EQUAL(lrows, 5);
  BOOST_CHECK_EQUAL(lcols, 5);

  BOOST_TEST_MESSAGE("Fill Matrix");

  ierr = MatGetOwnershipRange(A, &lo, &hi);  CHKERRXX(ierr);
  PetscScalar x(0.0);
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      ierr = MatSetValues(A, 1, &i, 1, &j, &x, INSERT_VALUES);  CHKERRXX(ierr);
      x += 1.0;
    }
  }

  BOOST_TEST_MESSAGE("Assemble Matrix");
  
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  CHKERRXX(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);  CHKERRXX(ierr);


  BOOST_TEST_MESSAGE("Get Matrix Values");

  PetscScalar y(0.0);
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      ierr = MatGetValues(A, 1, &i, 1, &j, &x);  CHKERRXX(ierr);
      BOOST_CHECK_EQUAL(x, y);
      y += 1.0;
    }
  }

  BOOST_TEST_MESSAGE("Add Matrix Values");

  x = 0.0;
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      ierr = MatSetValues(A, 1, &i, 1, &j, &x, ADD_VALUES);  CHKERRXX(ierr);
      x += 1.0;
    }
  }

  BOOST_TEST_MESSAGE("Assemble Matrix");

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  CHKERRXX(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);  CHKERRXX(ierr);

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

// -------------------------------------------------------------
// init_function
// -------------------------------------------------------------
bool init_function()
{
  return true;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  gridpack::parallel::Communicator world;
  gridpack::math::Initialize();
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  gridpack::math::Finalize();
  return result;
}

