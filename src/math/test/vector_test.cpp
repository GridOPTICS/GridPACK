/**
 * @file   vector_construction_test.cpp
 * @author William A. Perkins
 * @date   2013-05-08 15:11:55 d3g096
 * 
 * @brief  Construction/clone unit testing for gridpack::math::Vector
 * 
 * 
 */


// #include <iostream>
#include "gridpack/parallel/parallel.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/math/vector.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Vector)

BOOST_AUTO_TEST_CASE( construction )
{
  static int local_size(5);
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v(world, local_size);

  BOOST_CHECK_EQUAL( v.local_size(), local_size);
  BOOST_CHECK_EQUAL( v.size(), world.size()*local_size);
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
