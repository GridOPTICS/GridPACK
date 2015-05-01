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
 * @date   2015-05-01 12:35:30 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "math.hpp"
#include "gridpack/parallel/parallel.hpp"
#include "petsc/ga_matrix.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(GAMatrixTest)

BOOST_AUTO_TEST_CASE( construction )
{
  PetscErrorCode ierr(0);
  gridpack::parallel::Communicator world;

  Mat A;
  PetscInt lrows, lcols;
  PetscInt lo, hi;

  ierr = MatCreateDenseGA(world, 5, 5, PETSC_DETERMINE, PETSC_DETERMINE, &A); CHKERRXX(ierr);
  ierr = MatGetLocalSize(A, &lrows, &lcols);  CHKERRXX(ierr);

  BOOST_CHECK_EQUAL(lrows, 5);
  BOOST_CHECK_EQUAL(lcols, 5);

  ierr = MatGetOwnershipRange(A, &lo, &hi);  CHKERRXX(ierr);
  PetscScalar x(0.0);
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      ierr = MatSetValues(A, 1, &i, 1, &j, &x, INSERT_VALUES);  CHKERRXX(ierr);
      x += 1.0;
    }
  }
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  CHKERRXX(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);  CHKERRXX(ierr);

  PetscScalar y(0.0);
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      ierr = MatGetValues(A, 1, &i, 1, &j, &x);  CHKERRXX(ierr);
      BOOST_CHECK_EQUAL(x, y);
      y += 1.0;
    }
  }

  x = 0.0;
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      ierr = MatSetValues(A, 1, &i, 1, &j, &x, ADD_VALUES);  CHKERRXX(ierr);
      x += 1.0;
    }
  }
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  CHKERRXX(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);  CHKERRXX(ierr);

  y = 0.0;
  for (int i = lo; i < hi; ++i) {
    for (int j = lo; j < hi; ++j) {
      ierr = MatGetValues(A, 1, &i, 1, &j, &x);  CHKERRXX(ierr);
      BOOST_CHECK_EQUAL(x, 2.0*y);
      y += 1.0;
    }
  }


  ierr = MatDestroy(&A);
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

