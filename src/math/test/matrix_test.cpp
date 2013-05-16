/**
 * @file   matrix_test.cpp
 * @author William A. Perkins
 * @date   2013-05-16 11:12:41 d3g096
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


BOOST_AUTO_TEST_SUITE(Matrix)

BOOST_AUTO_TEST_CASE( sparse_construction )
{
  gridpack::parallel::Communicator world;
  int global_size;
  boost::mpi::all_reduce(world, local_size, global_size, std::plus<int>());

  gridpack::math::Matrix A(world, global_size, global_size);
  
}

BOOST_AUTO_TEST_CASE( dense_construction )
{
  gridpack::parallel::Communicator world;
  int global_size;
  boost::mpi::all_reduce(world, local_size, global_size, std::plus<int>());

  gridpack::math::Matrix A(world, global_size, global_size, 
                           gridpack::math::Matrix::Dense);
  
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
