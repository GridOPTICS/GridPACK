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
 * @date   2015-01-27 09:48:00 d3g096
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

static const double delta(0.0001);

BOOST_AUTO_TEST_SUITE(Numeric)

BOOST_AUTO_TEST_CASE(TypeChecking)
{
  using namespace gridpack::math;
  BOOST_CHECK(!TypeCheck<gridpack::RealType>::isSame<gridpack::ComplexType>::value);
  BOOST_CHECK(TypeCheck<gridpack::RealType>::isSame<gridpack::RealType>::value);
  BOOST_CHECK(TypeCheck<gridpack::ComplexType>::isSame<gridpack::ComplexType>::value);
  BOOST_CHECK(TypeCheck<gridpack::RealType>::check);
  BOOST_CHECK(TypeCheck<gridpack::ComplexType>::check);
  BOOST_CHECK(!TypeCheck<float>::check);
  BOOST_CHECK(!TypeCheck<int>::check);
  BOOST_CHECK(TypeCheck<double>::check);
  BOOST_CHECK(TypeCheck< std::complex<double> >::check); 
  BOOST_CHECK(!TypeCheck< std::complex<float> >::check); 
}

BOOST_AUTO_TEST_CASE(MPL) 
{
  using namespace gridpack;
  typedef 
    boost::mpl::bool_<
      boost::is_same<RealType, RealType>::value ||
      boost::is_same<RealType, ComplexType>::value
    >::type test1;
  BOOST_CHECK(test1::value);
}
  
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Transfer)

BOOST_AUTO_TEST_CASE(Real_Complex)
{
  gridpack::RealType r[2] = { 1.0, 1.0 };
  gridpack::ComplexType *c;
  gridpack::math::ValueTransfer<gridpack::RealType, gridpack::ComplexType> 
    trans(2, &r[0]);
  trans.go();
  c = trans.to();
  std::cout << *c << std::endl;
  BOOST_CHECK_CLOSE(std::real(*c), r[0], delta);
  BOOST_CHECK_CLOSE(std::imag(*c), r[1], delta);

  (*c) *= 2.0;
  gridpack::math::ValueTransfer<gridpack::ComplexType, gridpack::RealType> 
    rtrans(1, c, &r[0]); 
  rtrans.go();
  std::cout << "(" << r[0] << "," << r[1] << ")" << std::endl;


  gridpack::RealType *x;
  gridpack::math::ValueTransfer<gridpack::RealType, gridpack::RealType>
    real_trans(2, &r[0]);
  real_trans.go();
  x = real_trans.to();
  std::cout << "(" << x[0] << "," << x[1] << ")" << std::endl;

  gridpack::ComplexType cx;
  gridpack::math::ValueTransfer<gridpack::ComplexType, gridpack::ComplexType>
    complex_trans(1, c, &cx);
  complex_trans.go();
  std::cout << cx << std::endl;
}

BOOST_AUTO_TEST_CASE(StorageSize)
{
  using namespace gridpack;
  using namespace gridpack::math;

  BOOST_CHECK_EQUAL((storage_size<RealType, RealType>::value), 1);
  BOOST_CHECK_EQUAL((storage_size<ComplexType, ComplexType>::value), 1);
  BOOST_CHECK_EQUAL((storage_size<ComplexType, RealType>::value), 2);
  BOOST_CHECK_EQUAL((storage_size<RealType, ComplexType>::value), 1);
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
