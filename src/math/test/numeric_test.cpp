// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   numeric_test.cpp
 * @author William A. Perkins
 * @date   2015-01-20 10:55:18 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created January 20, 2015 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------

#include <iostream>

#include "numeric_type_check.hpp"
#include "value_transfer.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Numeric)

BOOST_AUTO_TEST_CASE(TypeCheck)
{
  BOOST_CHECK(gridpack::math::TypeCheck<gridpack::RealType>::check());
  BOOST_CHECK(gridpack::math::TypeCheck<gridpack::ComplexType>::check());
  BOOST_CHECK(!gridpack::math::TypeCheck<float>::check());
  BOOST_CHECK(!gridpack::math::TypeCheck<int>::check());
  BOOST_CHECK(gridpack::math::TypeCheck<double>::check());
  BOOST_CHECK(gridpack::math::TypeCheck< std::complex<double> >::check()); 
  BOOST_CHECK(!gridpack::math::TypeCheck< std::complex<float> >::check()); 
}
  
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Transfer)

BOOST_AUTO_TEST_CASE(Real_Complex)
{
  gridpack::RealType r[2] = { 1.0, 1.0 };
  gridpack::ComplexType *c;
  gridpack::math::ValueTransfer<gridpack::RealType, gridpack::ComplexType> trans(2, &r[0]);
  c = trans.to();
  std::cout << *c << std::endl;
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
  return result;
}
