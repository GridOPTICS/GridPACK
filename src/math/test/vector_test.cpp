/**
 * @file   vector_construction_test.cpp
 * @author William A. Perkins
 * @date   2013-08-13 12:22:31 d3g096
 * 
 * @brief  Construction/clone unit testing for gridpack::math::Vector
 * 
 * 
 */

#include <iostream>
#include <iterator>
#include "gridpack/utilities/exception.hpp"
#include "gridpack/parallel/parallel.hpp"
#include "math.hpp"
#include "vector.hpp"

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

  static const gridpack::ComplexType value(1.0, 1.0);

  v.fill(value);
  v.ready();
  
  int lo, hi;
  gridpack::ComplexType xlo, xhi;

  v.local_index_range(lo, hi);

  v.get_element(lo, xlo);
  v.get_element(hi-1, xhi);

  BOOST_CHECK_CLOSE(real(xlo), real(value), delta);
  BOOST_CHECK_CLOSE(abs(xlo), abs(value), delta);
  BOOST_CHECK_CLOSE(real(xhi), real(value), delta);
  BOOST_CHECK_CLOSE(abs(xhi), abs(value), delta);

  // std::cerr << "About to clone vector" << std::endl;
  std::auto_ptr<gridpack::math::Vector> vcopy(v.clone());

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
  std::vector<gridpack::ComplexType> val; 
  val.reserve(local_size);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x;
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
    gridpack::ComplexType x( static_cast<double>(i));
    gridpack::ComplexType x1;
    gridpack::ComplexType x2;

    v1.get_element(i, x1);
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
  std::vector<gridpack::ComplexType> val; 
  val.reserve(local_size);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x;
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
    gridpack::ComplexType x( static_cast<double>(i));
    gridpack::ComplexType x1;
    gridpack::ComplexType x2;

    v1.get_element(i, x1);
    v2.get_element(i, x2);

    BOOST_CHECK_CLOSE(2.0*real(x), real(x1), delta);
    BOOST_CHECK_CLOSE(2.0*real(x), real(x2), delta);
  }
}

BOOST_AUTO_TEST_CASE( get_all )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v(world, local_size);

  int lo, hi;
  v.local_index_range(lo, hi);
  
  std::vector<gridpack::ComplexType> x(hi-lo, static_cast<double>(world.rank()));
  v.set_element_range(lo, hi, &x[0]);

  std::vector<gridpack::ComplexType> all(world.size()*local_size);
  v.get_all_elements(&all[0]);

  for (int p = 0; p < world.size(); ++p) {
    if (p == world.rank()) {
      std::cout << p << ": ";
      std::copy(all.begin(), all.end(),
                std::ostream_iterator<gridpack::ComplexType>(std::cout, ", "));
      std::cout << std::endl;
    }
    world.barrier();
  }

  for (int p = 0; p < world.size(); ++p) {
    for (int i = 0; i < local_size; ++i) {
      int index(p*local_size + i);
      BOOST_CHECK_EQUAL(p, static_cast<int>(real(all[index])));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(VectorOperations)

BOOST_AUTO_TEST_CASE( add )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector 
    v1(world, local_size),
    v2(world, local_size);
  std::auto_ptr<gridpack::math::Vector> vsum1, vsum2;
  
  v1.fill(1.0);
  v2.fill(2.0);
  
  v1.ready();
  v2.ready();

  vsum1.reset(v1.clone());
  vsum1->add(v2);

  vsum2.reset(gridpack::math::add(v1, v2));

  int lo, hi;
  v1.local_index_range(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x1, x2, x3;

    v1.get_element(i, x1);
    v2.get_element(i, x2);
    x3 = x1 + x2;

    vsum1->get_element(i, x1);
    vsum2->get_element(i, x2);

    BOOST_CHECK_CLOSE(real(x1), real(x2), delta);
    BOOST_CHECK_CLOSE(real(x1), real(x3), delta);
  }

  // try using a vector that is the wrong size
  gridpack::math::Vector v3(world, local_size-1);
  BOOST_CHECK_THROW(v1.add(v3), gridpack::Exception);
}

BOOST_AUTO_TEST_CASE( add_or_scale_scalar )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v1(world, local_size), v2(world, local_size);
  gridpack::ComplexType x(2.0), offset(2.0), scale(3.0);

  v1.fill(x);
  v1.add(offset);

  v2.fill(x);
  v2.scale(scale);

  int lo, hi;
  v1.local_index_range(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType y;
    v1.get_element(i, y);
    BOOST_CHECK_CLOSE(real(x+offset), real(y), delta);
    BOOST_CHECK_CLOSE( abs(x+offset), abs(y), delta);
    v2.get_element(i, y);
    BOOST_CHECK_CLOSE(real(x*scale), real(y), delta);
    BOOST_CHECK_CLOSE( abs(x*scale), abs(y), delta);
  }
}


BOOST_AUTO_TEST_CASE( reciprocal )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v1(world, local_size);
  gridpack::ComplexType x(2.0, 5.0);
  gridpack::ComplexType rx(1.0/x);

  v1.fill(x);
  v1.reciprocal();

  int lo, hi;
  v1.local_index_range(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType y;
    v1.get_element(i, y);
    BOOST_CHECK_CLOSE(real(rx), real(y), delta);
    BOOST_CHECK_CLOSE( abs(rx), abs(y), delta);
  }
}

BOOST_AUTO_TEST_CASE( print )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v1(world, local_size);
  gridpack::ComplexType x(2.0, 1.0);
  v1.fill(x);
  v1.print();

  std::string out;
  if (world.size() > 1) {
    out = "vector_parallel.out";
  } else {
    out = "vector_serial.out";
  }
  v1.print(out.c_str());

  if (world.size() > 1) {
    out = "vector_parallel.mat";
  } else {
    out = "vector_serial.mat";
  }
  v1.save(out.c_str());
  
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
