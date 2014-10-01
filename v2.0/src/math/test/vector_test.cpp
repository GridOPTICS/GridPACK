// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   vector_construction_test.cpp
 * @author William A. Perkins
 * @date   2014-03-05 12:11:53 d3g096
 * 
 * @brief  Unit tests for gridpack::math::Vector
 * 
 * @test
 */
// -------------------------------------------------------------

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

BOOST_AUTO_TEST_SUITE(VectorTest)

BOOST_AUTO_TEST_CASE( construction )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v(world, local_size);

  BOOST_CHECK_EQUAL( v.localSize(), local_size);
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

  v.localIndexRange(lo, hi);

  v.getElement(lo, xlo);
  v.getElement(hi-1, xhi);

  BOOST_CHECK_CLOSE(real(xlo), real(value), delta);
  BOOST_CHECK_CLOSE(abs(xlo), abs(value), delta);
  BOOST_CHECK_CLOSE(real(xhi), real(value), delta);
  BOOST_CHECK_CLOSE(abs(xhi), abs(value), delta);

  // std::cerr << "About to clone vector" << std::endl;
  std::auto_ptr<gridpack::math::Vector> vcopy(v.clone());

  vcopy->getElement(lo, xlo);
  vcopy->getElement(hi-1, xhi);

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
  v1.localIndexRange(lo, hi);
  
  std::vector<int> idx; 
  idx.reserve(local_size);
  std::vector<gridpack::ComplexType> val; 
  val.reserve(local_size);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x;
    x = static_cast<double>(i);
    v1.setElement(i, x);
    idx.push_back(i);
    val.push_back(x);
  }
  v2.setElements(idx.size(), &idx[0], &val[0]);

  v1.ready();
  v2.ready();

  val.clear();
  v1.getElements(idx.size(), &idx[0], &val[0]);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x( static_cast<double>(i));
    gridpack::ComplexType x1;
    gridpack::ComplexType x2;

    v1.getElement(i, x1);
    v2.getElement(i, x2);

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
  v1.localIndexRange(lo, hi);
  
  std::vector<int> idx; 
  idx.reserve(local_size);
  std::vector<gridpack::ComplexType> val; 
  val.reserve(local_size);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x;
    x = static_cast<double>(i);
    v1.setElement(i, x);
    v1.addElement(i, x);
    idx.push_back(i);
    val.push_back(x);
  }
  v2.setElements(idx.size(), &idx[0], &val[0]);
  v2.addElements(idx.size(), &idx[0], &val[0]);

  v1.ready();
  v2.ready();

  val.clear();
  v1.getElements(idx.size(), &idx[0], &val[0]);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x( static_cast<double>(i));
    gridpack::ComplexType x1;
    gridpack::ComplexType x2;

    v1.getElement(i, x1);
    v2.getElement(i, x2);

    BOOST_CHECK_CLOSE(2.0*real(x), real(x1), delta);
    BOOST_CHECK_CLOSE(2.0*real(x), real(x2), delta);
  }
}

BOOST_AUTO_TEST_CASE( get_all )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v(world, local_size);

  int lo, hi;
  v.localIndexRange(lo, hi);
  
  std::vector<gridpack::ComplexType> x(hi-lo, static_cast<double>(world.rank()));
  v.setElementRange(lo, hi, &x[0]);

  std::vector<gridpack::ComplexType> all(world.size()*local_size);
  v.getAllElements(&all[0]);

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

BOOST_AUTO_TEST_CASE( complex_operations )
{
  gridpack::parallel::Communicator world;
  gridpack::ComplexType x0(static_cast<double>(world.rank()+1),
                           static_cast<double>(world.rank()+1));

  gridpack::math::Vector v(world, local_size);

  int lo, hi;
  v.localIndexRange(lo, hi);
  
  std::vector<gridpack::ComplexType> x(hi-lo, x0);
  v.setElementRange(lo, hi, &x[0]);
  v.ready();
  
  std::auto_ptr<gridpack::math::Vector> 
    vreal(gridpack::math::real(v)),
    vimag(gridpack::math::imaginary(v)),
    vconj(gridpack::math::conjugate(v)),
    vabs(gridpack::math::abs(v));

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType z, zreal, zimag, zconj, zabs;
    v.getElement(i, z);
    vreal->getElement(i, zreal);
    vimag->getElement(i, zimag);
    vconj->getElement(i, zconj);
    vabs->getElement(i, zabs);

    BOOST_CHECK_CLOSE(real(z), abs(zreal), delta);
    BOOST_CHECK_CLOSE(imag(z), abs(zimag), delta);
    BOOST_CHECK_CLOSE(abs(z), abs(zconj), delta);
    BOOST_CHECK_CLOSE(imag(z), -imag(zconj), delta);
    BOOST_CHECK_CLOSE(abs(z), real(zabs), delta);
  }
    
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(VectorOperationsTest)

BOOST_AUTO_TEST_CASE( add )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector 
    v1(world, local_size),
    v2(world, local_size),
    v3(world, local_size);
  std::auto_ptr<gridpack::math::Vector> vsum1, vsum2, vsum3;
  
  v1.fill(1.0);
  v2.fill(2.0);
  v3.fill(2.0);
  
  v1.ready();
  v2.ready();
  v3.ready();

  vsum1.reset(v1.clone());
  vsum1->add(v2);

  vsum2.reset(gridpack::math::add(v1, v2));

  vsum3.reset(v1.clone());
  vsum3->add(v3, 2.0);

  int lo, hi;
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x1, x2, x3;

    v1.getElement(i, x1);
    v2.getElement(i, x2);
    x3 = x1 + x2;

    vsum1->getElement(i, x1);
    vsum2->getElement(i, x2);

    BOOST_CHECK_CLOSE(real(x1), real(x2), delta);
    BOOST_CHECK_CLOSE(real(x1), real(x3), delta);

    v1.getElement(i, x1);
    v3.getElement(i, x2);
    x2 = x1 + 2.0*x2;

    vsum3->getElement(i, x3);
    BOOST_CHECK_CLOSE(real(x2), real(x3), delta);
  }

  // try using a vector that is the wrong size
  gridpack::math::Vector v4(world, local_size-1);
  BOOST_CHECK_THROW(v1.add(v4), gridpack::Exception);
}

BOOST_AUTO_TEST_CASE( add_or_scale_scalar )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v1(world, local_size), v2(world, local_size);
  gridpack::ComplexType x(2.0), offset(2.0), scale(3.0);

  v1.fill(x);
  v1.add(offset);
  v1.print();

  v2.fill(x);
  v2.scale(scale);
  v2.print();

  int lo, hi;
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType y;
    v1.getElement(i, y);
    BOOST_CHECK_CLOSE(real(x+offset), real(y), delta);
    BOOST_CHECK_CLOSE( abs(x+offset), abs(y), delta);
    v2.getElement(i, y);
    BOOST_CHECK_CLOSE(real(x*scale), real(y), delta);
    BOOST_CHECK_CLOSE( abs(x*scale), abs(y), delta);
  }
}

BOOST_AUTO_TEST_CASE( subtract_scalar )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v1(world, local_size);
  gridpack::ComplexType x(1.0);

  v1.fill(x);
  v1.add(-x);
  v1.print();

  int lo, hi;
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType y;
    v1.getElement(i, y);
    BOOST_CHECK_CLOSE(real(y), 0.0, delta);
    BOOST_CHECK_CLOSE( abs(y), 0.0, delta);
  }
}

BOOST_AUTO_TEST_CASE( elementwise )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v1(world, local_size);
  gridpack::math::Vector v2(world, local_size);
  gridpack::math::Vector v3(world, local_size);
  gridpack::ComplexType x(2.0), y(0.5);

  v1.fill(x);
  v2.fill(x);
  v3.fill(y);

  v1.elementMultiply(v3);
  v2.elementDivide(v3);

  int lo, hi;
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType z1, z2;
    v1.getElement(i, z1);
    v2.getElement(i, z2);
    BOOST_CHECK_CLOSE(real(z1), real(x*y), delta);
    BOOST_CHECK_CLOSE( abs(z1), abs(x*y), delta);
    BOOST_CHECK_CLOSE(real(z2), real(x/y), delta);
    BOOST_CHECK_CLOSE( abs(z2), abs(x/y), delta);
  }
}

BOOST_AUTO_TEST_CASE( elementwise2 )
{
  gridpack::parallel::Communicator world;
  static int three(3);
  gridpack::math::Vector v1(world, three);
  gridpack::math::Vector v2(world, three);
  gridpack::math::Vector v3(world, three);
  gridpack::math::Vector v4(world, three);

  gridpack::ComplexType x1[three];
  x1[0] = gridpack::ComplexType(0.348262, 3.4343 );
  x1[1] = gridpack::ComplexType(1.50794,  2.76069);
  x1[2] = gridpack::ComplexType(1.04059,  4.50791);

  gridpack::ComplexType x2[three];
  x2[0] = gridpack::ComplexType(1.08099,  0.0391692);
  x2[1] = gridpack::ComplexType(1.04585,  0.378384);
  x2[2] = gridpack::ComplexType(1.08145,  0.249638);

  int lo, hi;
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    v1.setElement(i, x1[i]);
    v2.setElement(i, x2[i]);
    v3.setElement(i, x1[i]*x2[i]);
  }
  v1.ready();
  v2.ready();
  v3.ready();

  v3.print();
  
  v4.equate(v1);
  v4.elementMultiply(v2);

  v4.print();

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x, y;
    v3.getElement(i, x);
    v4.getElement(i, y);
    BOOST_CHECK_CLOSE(real(y), real(x), delta);
    BOOST_CHECK_CLOSE( abs(y), abs(x), delta);
  }

  // it should be commutative

  v4.equate(v2);
  v4.elementMultiply(v1);
  v4.print();

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x, y;
    v3.getElement(i, x);
    v4.getElement(i, y);
    BOOST_CHECK_CLOSE(real(y), real(x), delta);
    BOOST_CHECK_CLOSE( abs(y), abs(x), delta);
  }
}

BOOST_AUTO_TEST_CASE( exp )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v1(world, local_size);
  gridpack::ComplexType x(7.3);

  v1.fill(x);
  v1.exp();

  int lo, hi;
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType y;
    v1.getElement(i, y);
    BOOST_CHECK_CLOSE(real(log(y)), real(x), delta);
    BOOST_CHECK_CLOSE( abs(log(y)), abs(x), delta);
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
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType y;
    v1.getElement(i, y);
    BOOST_CHECK_CLOSE(real(rx), real(y), delta);
    BOOST_CHECK_CLOSE( abs(rx), abs(y), delta);
  }
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( VectorIOTest )

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

BOOST_AUTO_TEST_CASE( load_save )
{
  gridpack::parallel::Communicator world;
  gridpack::math::Vector v1(world, local_size);
  gridpack::ComplexType x(2.0, 1.0);
  v1.fill(x);
  v1.print();

  std::string out;
  if (world.size() > 1) {
    out = "vector_binary_parallel.out";
  } else {
    out = "vector_binary_serial.out";
  }

  v1.saveBinary(out.c_str());

  gridpack::math::Vector v2(world, local_size);
  v2.loadBinary(out.c_str());

  int lo, hi;
  v2.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x, y;
    v1.getElement(i, x);
    v2.getElement(i, y);
    BOOST_CHECK_CLOSE(real(y), real(x), delta);
    BOOST_CHECK_CLOSE( abs(y), abs(x), delta);
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
  return result;
}
