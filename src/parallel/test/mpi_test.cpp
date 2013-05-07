/**
 * @file   mpi_test.cpp
 * @author William A. Perkins
 * @date   2013-05-07 13:48:20 d3g096
 * 
 * @brief  A simple MPI unit test to serve as an example
 * 
 * 
 */

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE ( test ) 

BOOST_AUTO_TEST_CASE( test )
{
  boost::mpi::communicator world;
  int i(1), sumi(0);
  boost::mpi::all_reduce(world, i, sumi, std::plus<int>());
  BOOST_CHECK_EQUAL(sumi, static_cast<int>(world.size()));
}

BOOST_AUTO_TEST_SUITE_END( )

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
  boost::mpi::environment env(argc, argv);
  return ::boost::unit_test::unit_test_main( &init_function, argc, argv );
}
  
