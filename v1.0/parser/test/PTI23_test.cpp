/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   PTI23_test.cpp
 * @author Kevin Glass
 * @date   2013-06-17 12:09:38 d3k427
 * 
 * @brief  Test PTI23_parser capability. Currently not implemented.
 * 
 * 
 */

#include <iostream>
#include <string>
//#include <gridpack/parser/Parser.hpp>
#include <gridpack/parser/PTI23_parser.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#define EPSILON     0.0000001
#define TOLERANCE(x, y, eps) (x - EPSILON > y && x + EPSILON < y)

//BOOST_AUTO_TEST_SUITE(Parser)

//____________________________________________________________________________//

//BOOST_AUTO_TEST_SUITE_END()

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


