/**
 * @file   vector_construction_test.cpp
 * @author William A. Perkins
 * @date   2013-06-04 13:04:23 d3g096
 * 
 * @brief  Construction/clone unit testing for gridpack::math::Vector
 * 
 * 
 */

#include <iostream>
#include "PTI23_parser.hpp"

#define BOOST_TEST_NO_MAIN
//#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>


//____________________________________________________________________________//

// most frequently you implement test cases as a free functions with automatic registration
BOOST_AUTO_TEST_CASE( test1 )
{
    // reports 'error in "test1": test 2 == 1 failed'
    BOOST_CHECK( 2 == 1 );
}

//____________________________________________________________________________//

// each test file may contain any number of test cases; each test case has to have unique name
BOOST_AUTO_TEST_CASE( test2 )
{
    int i = 0;

    // reports 'error in "test2": check i == 2 failed [0 != 2]'
    BOOST_CHECK_EQUAL( i, 2 );

    BOOST_CHECK_EQUAL( i, 0 );
}

//____________________________________________________________________________//

// EOF
/*
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
*/
