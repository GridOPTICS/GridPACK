/**
 * @file   vector_construction_test.cpp
 * @author William A. Perkins
 * @date   2013-05-10 12:02:29 d3g096
 * 
 * @brief  Construction/clone unit testing for gridpack::math::Vector
 * 
 * 
 */

#include <iostream>
#include "gridpack/parallel/parallel.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/math/vector.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Vector)

BOOST_AUTO_TEST_CASE( construction )
{
  static const int local_size(5);
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v(world, local_size);

  BOOST_CHECK_EQUAL( v.local_size(), local_size);
  BOOST_CHECK_EQUAL( v.size(), world.size()*local_size);
}

BOOST_AUTO_TEST_CASE( fill_and_clone )
{
  static const int local_size(5);
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v(world, local_size);

  static const gridpack::math::complex_type value(1.0, 1.0);

  v.fill(value);
  v.ready();
  
  int lo, hi;
  gridpack::math::complex_type xlo, xhi;

  v.local_index_range(lo, hi);

  v.get_element(lo, xlo);
  v.get_element(hi-1, xhi);

  BOOST_CHECK_CLOSE(real(xlo), real(value), 0.001);
  BOOST_CHECK_CLOSE(abs(xlo), abs(value), 0.001);
  BOOST_CHECK_CLOSE(real(xhi), real(value), 0.001);
  BOOST_CHECK_CLOSE(abs(xhi), abs(value), 0.001);

  std::cerr << "About to clone vector" << std::endl;
  std::auto_ptr<gridpack::math::Vector> 
    vcopy(clone(v));

  vcopy->get_element(lo, xlo);
  vcopy->get_element(hi-1, xhi);

  BOOST_CHECK_CLOSE(real(xlo), real(value), 0.001);
  BOOST_CHECK_CLOSE(abs(xlo), abs(value), 0.001);
  BOOST_CHECK_CLOSE(real(xhi), real(value), 0.001);
  BOOST_CHECK_CLOSE(abs(xhi), abs(value), 0.001);
}

BOOST_AUTO_TEST_CASE( set_and_get )
{
  static const int local_size(5);
  gridpack::parallel::Communicator world;
  gridpack::math::Vector 
    v1(world, local_size),
    v2(world, local_size);

  int lo, hi;
  v1.local_index_range(lo, hi);
  
  std::vector<int> idx; 
  idx.reserve(local_size);
  std::vector<gridpack::math::complex_type> val; 
  val.reserve(local_size);

  for (int i = lo; i < hi; ++i) {
    gridpack::math::complex_type x;
    x = static_cast<double>(i);
    v1.set_element(i, x);
    idx.push_back(i);
    val.push_back(x);
  }
  v2.set_elements(idx.size(), &idx[0], &val[0]);

  v1.ready();
  v2.ready();

  val.clear();
  v1.get_elements(idx.size(), &idx[0], &val[0]);

  for (int i = lo; i < hi; ++i) {
    gridpack::math::complex_type x( static_cast<double>(i));
    gridpack::math::complex_type x1(val[i]);
    gridpack::math::complex_type x2;

    v2.get_element(i, x2);

    BOOST_CHECK_CLOSE(real(x), real(x1), 0.001);
    BOOST_CHECK_CLOSE(real(x), real(x2), 0.001);
  }
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
