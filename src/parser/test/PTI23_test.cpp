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
#include <string>
#include <gridpack/parser/Parser.hpp>
#include <gridpack/parser/PTI23_parser.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Parser)

//____________________________________________________________________________//


BOOST_AUTO_TEST_CASE(openFailure)
{
    bool                    opened         = true;
    try { 
	std::string         fileName        = "";
        gridpack::parser::Parser<gridpack::parser::PTI23_parser> parser;
	parser.getCaseData(fileName);
    } catch (gridpack::Exception & e) {
        opened     = false;
    }

    BOOST_CHECK_EQUAL(opened, false);
}

BOOST_AUTO_TEST_CASE(openSuccess)
{
    bool                    opened        = false;
    try {
        std::string          fileName      = "../PTI_seqtest.raw";
        gridpack::parser::Parser<gridpack::parser::PTI23_parser> parser;
        parser.getCaseData(fileName);
    } catch (gridpack::Exception & e) {
        opened     = false;
    }

    BOOST_CHECK_EQUAL(opened, true);
}

BOOST_AUTO_TEST_CASE(readValidFile)
{
    bool                    opened        = false;
    try {
        std::string          fileName      = "../PTI_seqtest.raw";
        gridpack::parser::Parser<gridpack::parser::PTI23_parser> parser;
        parser.getCaseData(fileName);
	
    } catch (gridpack::Exception & e) {
        opened     = false;
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
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
}


