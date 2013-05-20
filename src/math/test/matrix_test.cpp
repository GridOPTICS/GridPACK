/**
 * @file   matrix_test.cpp
 * @author William A. Perkins
 * @date   2013-05-20 10:52:42 d3g096
 * 
 * @brief  Unit tests for Matrix
 * 
 * 
 */

#include <iostream>
#include <boost/mpi/collectives.hpp>
#include "gridpack/parallel/parallel.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/math/matrix.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

static const int local_size(5);
static const double delta(0.0001);

static const gridpack::math::Matrix::StorageType the_storage_type =
#ifdef TEST_DENSE
  gridpack::math::Matrix::Dense;
#else
  gridpack::math::Matrix::Sparse;
#endif

BOOST_AUTO_TEST_SUITE(Matrix)

BOOST_AUTO_TEST_CASE( construction )
{
  gridpack::parallel::Communicator world;
  int global_size;
  boost::mpi::all_reduce(world, local_size, global_size, std::plus<int>());

  switch (the_storage_type) {
  case gridpack::math::Matrix::Dense:
    BOOST_TEST_MESSAGE("Testing Dense Matrices ...");
    break;
  case gridpack::math::Matrix::Sparse:
    BOOST_TEST_MESSAGE("Testing Sparse Matrices ...");
    break;
  default:
    BOOST_FAIL("Undefined Matrix storage type");
  }

  gridpack::math::Matrix A(world, local_size, global_size, the_storage_type);

  int lo, hi;
  A.local_row_range(lo, hi);

  BOOST_CHECK_GT(hi, lo);

  BOOST_CHECK_EQUAL(A.local_rows(), local_size);
  BOOST_CHECK_EQUAL(A.rows(), global_size);
  BOOST_CHECK_EQUAL(A.cols(), global_size);
  
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
  gridpack::math::Initialize();
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  gridpack::math::Finalize();
}
