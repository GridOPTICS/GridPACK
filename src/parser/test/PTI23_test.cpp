/**
 * @file   vector_construction_test.cpp
 * @author William A. Perkins
 * @date   2013-06-17 12:09:38 d3g096
 * 
 * @brief  Construction/clone unit testing for gridpack::math::Vector
 * 
 * 
 */

#include <iostream>
#include "parser.hpp"
#include "PTI23_parser.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Parser)

//____________________________________________________________________________//

// most frequently you implement test cases as a free functions with automatic registration
BOOST_AUTO_TEST_CASE(openFailure)
{
    bool                    opened         = true;
    try {
        Parser<PTI22_parser> parser("");
    } catch (ios_base::failure & e) {
        opened     = false;
    }

    BOOST_CHECK_EQUAL(opened, false);

}

BOOST_AUTO_TEST_CASE(openSuccess)
{
    bool                    opened        = false;
    try {
        Parser<PTI22_parser> parser("PTI23_seqtest.raw");
    } catch (ios_base::failure & e) {
        opened     = true;
    }

    BOOST_CHECK_EQUAL(opened, true);

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
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
}


