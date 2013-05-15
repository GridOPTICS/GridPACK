/**
 * @file   vector_construction_test.cpp
 * @author William A. Perkins
 * @date   2013-05-15 11:33:40 d3g096
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

static const int local_size(5);
static const double delta(0.0001);

BOOST_AUTO_TEST_SUITE(Vector)

BOOST_AUTO_TEST_CASE( construction )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v(world, local_size);

  BOOST_CHECK_EQUAL( v.local_size(), local_size);
  BOOST_CHECK_EQUAL( v.size(), world.size()*local_size);
}

BOOST_AUTO_TEST_CASE( fill_and_clone )
{
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

  BOOST_CHECK_CLOSE(real(xlo), real(value), delta);
  BOOST_CHECK_CLOSE(abs(xlo), abs(value), delta);
  BOOST_CHECK_CLOSE(real(xhi), real(value), delta);
  BOOST_CHECK_CLOSE(abs(xhi), abs(value), delta);

  // std::cerr << "About to clone vector" << std::endl;
  std::auto_ptr<gridpack::math::Vector> 
    vcopy(gridpack::math::clone(v));

  vcopy->get_element(lo, xlo);
  vcopy->get_element(hi-1, xhi);

  BOOST_CHECK_CLOSE(real(xlo), real(value), delta);
  BOOST_CHECK_CLOSE(abs(xlo), abs(value), delta);
  BOOST_CHECK_CLOSE(real(xhi), real(value), delta);
  BOOST_CHECK_CLOSE(abs(xhi), abs(value), delta);
}

BOOST_AUTO_TEST_CASE( set_and_get_element )
{
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

    BOOST_CHECK_CLOSE(real(x), real(x1), delta);
    BOOST_CHECK_CLOSE(real(x), real(x2), delta);
  }
}

BOOST_AUTO_TEST_CASE( add_and_get_element )
{
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
    v1.add_element(i, x);
    idx.push_back(i);
    val.push_back(x);
  }
  v2.set_elements(idx.size(), &idx[0], &val[0]);
  v2.add_elements(idx.size(), &idx[0], &val[0]);

  v1.ready();
  v2.ready();

  val.clear();
  v1.get_elements(idx.size(), &idx[0], &val[0]);

  for (int i = lo; i < hi; ++i) {
    gridpack::math::complex_type x( static_cast<double>(i));
    gridpack::math::complex_type x1(val[i]);
    gridpack::math::complex_type x2;

    v2.get_element(i, x2);

    BOOST_CHECK_CLOSE(2.0*real(x), real(x1), delta);
    BOOST_CHECK_CLOSE(2.0*real(x), real(x2), delta);
  }
}

BOOST_AUTO_TEST_CASE( add )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector 
    v1(world, local_size),
    v2(world, local_size),
    *vsum1, *vsum2;
  
  v1.fill(1.0);
  v2.fill(2.0);
  
  v1.ready();
  v2.ready();

  vsum1 = gridpack::math::clone(v1);
  vsum1->add(v2);

  vsum2 = gridpack::math::add(v1, v2);

  int lo, hi;
  v1.local_index_range(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::math::complex_type x1, x2, x3;

    v1.get_element(i, x1);
    v2.get_element(i, x2);
    x3 = x1 + x2;

    vsum1->get_element(i, x1);
    vsum2->get_element(i, x2);

    BOOST_CHECK_CLOSE(real(x1), real(x2), delta);
    BOOST_CHECK_CLOSE(real(x1), real(x3), delta);
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
