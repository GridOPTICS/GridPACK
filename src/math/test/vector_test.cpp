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
 * @date   2019-11-20 11:51:55 d3g096
 * 
 * @brief  Unit tests for gridpack::math::Vector
 * 
 * @test
 */
// -------------------------------------------------------------

#include <iostream>
#include <iterator>
#include <boost/scoped_ptr.hpp>
#include "vector.hpp"

#include "test_main.cpp"

static const int local_size(5);
static const double delta(0.0001);

#ifdef TEST_REAL

typedef gridpack::RealType TestType;

#define TEST_VALUE_CLOSE(x, y, delta) \
  BOOST_CHECK_CLOSE((y), (x), delta);

#define TEST_VALUE(r, i) (r)

#else  // TEST_REAL

typedef gridpack::ComplexType TestType;
#define TEST_VALUE_CLOSE(x, y, delta) \
  BOOST_CHECK_CLOSE(real(y), real(x), delta); \
  BOOST_CHECK_CLOSE( abs(y), abs(x), delta);

#define TEST_VALUE(r, i) TestType(r,i)

#endif  // TEST_REAL


BOOST_AUTO_TEST_SUITE(VectorTest)

BOOST_AUTO_TEST_CASE( construction )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> v(world, local_size);

  BOOST_CHECK_EQUAL( v.localSize(), local_size);
  BOOST_CHECK_EQUAL( v.size(), world.size()*local_size);
}

BOOST_AUTO_TEST_CASE( fill_and_clone )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> v(world, local_size);

  static const TestType value(TEST_VALUE(1.0, 1.0));

  v.fill(value);
  v.ready();
  
  int lo, hi;
  TestType xlo, xhi;

  v.localIndexRange(lo, hi);

  v.getElement(lo, xlo);
  v.getElement(hi-1, xhi);

  TEST_VALUE_CLOSE(xlo, value, delta);
  TEST_VALUE_CLOSE(xhi, value, delta);

  // std::cerr << "About to clone vector" << std::endl;
  boost::scoped_ptr< gridpack::math::VectorT<TestType> > vcopy(v.clone());

  vcopy->getElement(lo, xlo);
  vcopy->getElement(hi-1, xhi);

  TEST_VALUE_CLOSE(xlo, value, delta);
  TEST_VALUE_CLOSE(xhi, value, delta);
}

BOOST_AUTO_TEST_CASE( set_and_get_element )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> 
    v1(world, local_size),
    v2(world, local_size);

  int lo, hi;
  v1.localIndexRange(lo, hi);
  
  std::vector<int> idx; 
  idx.reserve(local_size);
  std::vector<TestType> val; 
  val.reserve(local_size);

  for (int i = lo; i < hi; ++i) {
    TestType x;
    x = static_cast<double>(i);
    v1.setElement(i, x);
    idx.push_back(i);
    val.push_back(x);
  }
  v2.setElements(idx.size(), &idx[0], &val[0]);

  v1.ready();
  v2.ready();

  val.clear();
  val.reserve(local_size);
  v1.getElements(idx.size(), &idx[0], &val[0]);

  for (int i = lo; i < hi; ++i) {
    TestType x( static_cast<double>(i));
    TestType x1;
    TestType x2;

    v1.getElement(i, x1);
    v2.getElement(i, x2);

    TEST_VALUE_CLOSE(x, x1, delta);
    TEST_VALUE_CLOSE(x, x2, delta);
  }
}

BOOST_AUTO_TEST_CASE( add_and_get_element )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> 
    v1(world, local_size),
    v2(world, local_size);

  int lo, hi;
  v1.localIndexRange(lo, hi);
  
  std::vector<int> idx; 
  idx.reserve(local_size);
  std::vector<TestType> val; 
  val.reserve(local_size);

  for (int i = lo; i < hi; ++i) {
    TestType x;
    x = static_cast<double>(i);
    v1.setElement(i, x);
    idx.push_back(i);
    val.push_back(x);
  }
  v2.setElements(idx.size(), &idx[0], &val[0]);

  v1.ready();
  v2.ready();

  for (int i = lo; i < hi; ++i) {
    TestType x;
    x = static_cast<double>(i);
    v1.addElement(i, x);
  }
  v2.addElements(idx.size(), &idx[0], &val[0]);

  v1.ready();
  v2.ready();

  val.clear();
  v1.getElements(idx.size(), &idx[0], &val[0]);

  for (int i = lo; i < hi; ++i) {
    TestType x( static_cast<double>(i));
    TestType x1;
    TestType x2;

    v1.getElement(i, x1);
    v2.getElement(i, x2);

    TEST_VALUE_CLOSE(2.0*x, x1, delta);
    TEST_VALUE_CLOSE(2.0*x, x2, delta);
  }
}

BOOST_AUTO_TEST_CASE( set_add_get_range )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> 
    v1(world, local_size),
    v2(world, local_size);

  int lo, hi;
  v1.localIndexRange(lo, hi);
  
  std::vector<TestType> in, out; 
  in.reserve(local_size);

  for (int i = lo; i < hi; ++i) {
    TestType x;
    x = static_cast<double>(i);
    in.push_back(x);
  }
  v1.setElementRange(lo, hi, &in[0]);
  v1.ready();
  v1.addElementRange(lo, hi, &in[0]);
  v1.ready();

  out.reserve(local_size);
  v1.getElementRange(lo, hi, &out[0]);

  for (int i = lo; i < hi; ++i) {
    TestType x1;

    v1.getElement(i, x1);

    TEST_VALUE_CLOSE(2.0*in[i-lo], x1, delta);
    TEST_VALUE_CLOSE(out[i-lo], x1, delta);
  }
}

BOOST_AUTO_TEST_CASE( get_local )
{
  gridpack::parallel::Communicator world;
  int me(world.rank());
  double p(me + 1);

  if (me == 0) BOOST_TEST_MESSAGE( "Vector get local test:" );

  gridpack::math::VectorT<TestType> v(world, local_size);

  v.fill((double)0.0);
  if (me == 0) BOOST_TEST_MESSAGE( "Before vector:" );
  v.print();
  int lsize(v.localSize());
  TestType *l = v.getLocalElements();

  for (int i = 0; i < lsize; ++i) {
    l[i] = p;
  }
  v.releaseLocalElements(l);
  if (me == 0) BOOST_TEST_MESSAGE( "After vector:" );
  v.print();

  int lo, hi;
  v.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    TestType x, y(p);
    v.getElement(i, x);
    TEST_VALUE_CLOSE(x, y, delta);
  }
}

BOOST_AUTO_TEST_CASE( get_all )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> v(world, local_size);

  int lo, hi;
  v.localIndexRange(lo, hi);
  
  std::vector<TestType> x(hi-lo, static_cast<double>(world.rank()));
  v.setElementRange(lo, hi, &x[0]);

  std::vector<TestType> all(world.size()*local_size);
  v.getAllElements(&all[0]);

  for (int p = 0; p < world.size(); ++p) {
    if (p == world.rank()) {
      std::cout << p << ": ";
      std::copy(all.begin(), all.end(),
                std::ostream_iterator<TestType>(std::cout, ", "));
      std::cout << std::endl;
    }
    world.barrier();
  }

  for (int p = 0; p < world.size(); ++p) {
    for (int i = 0; i < local_size; ++i) {
      int index(p*local_size + i);
      gridpack::ComplexType x = all[index];
      BOOST_CHECK_EQUAL(p, static_cast<int>(real(x)));
    }
  }
}

BOOST_AUTO_TEST_CASE( local_clone )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> v1(world, local_size);

  int lo, hi;
  v1.localIndexRange(lo, hi);
  
  for (int i = lo; i < hi; ++i) {
    TestType x(TEST_VALUE(static_cast<double>(i), static_cast<double>(i)));
    v1.setElement(i, x);
  }
  v1.ready();

  boost::scoped_ptr< gridpack::math::VectorT<TestType> > 
    v2(v1.localClone());

  BOOST_CHECK_EQUAL(v2->processor_size(), 1);
  
  for (int i = lo; i < hi; ++i) {
    TestType x, y;
    v1.getElement(i, x);
    v2->getElement(i, y);

    TEST_VALUE_CLOSE(x, y, delta);
  }
}

BOOST_AUTO_TEST_CASE( complex_operations )
{
  gridpack::parallel::Communicator world;
  TestType x0(TEST_VALUE(static_cast<double>(world.rank()+1),
                         static_cast<double>(world.rank()+1)));

  gridpack::math::VectorT<TestType> v(world, local_size);

  int lo, hi;
  v.localIndexRange(lo, hi);
  
  std::vector<TestType> x(hi-lo, x0);
  v.setElementRange(lo, hi, &x[0]);
  v.ready();
  
  boost::scoped_ptr< gridpack::math::VectorT<TestType> > 
    vreal(gridpack::math::real(v)),
    vimag(gridpack::math::imaginary(v)),
    vconj(gridpack::math::conjugate(v)),
    vabs(gridpack::math::abs(v));

  for (int i = lo; i < hi; ++i) {
    TestType y;
    gridpack::ComplexType  z, zreal, zimag, zconj, zabs;
    v.getElement(i, y); z = y;
    vreal->getElement(i, y); zreal = y;
    vimag->getElement(i, y); zimag = y;
    vconj->getElement(i, y); zconj = y;
    vabs->getElement(i, y); zabs = y;

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
  gridpack::math::VectorT<TestType> 
    v1(world, local_size),
    v2(world, local_size),
    v3(world, local_size);
  boost::scoped_ptr< gridpack::math::VectorT<TestType> > vsum1, vsum2, vsum3;
  
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
    TestType x1, x2, x3;

    v1.getElement(i, x1);
    v2.getElement(i, x2);
    x3 = x1 + x2;

    vsum1->getElement(i, x1);
    vsum2->getElement(i, x2);

    TEST_VALUE_CLOSE(x1, x2, delta);

    v1.getElement(i, x1);
    v3.getElement(i, x2);
    x2 = x1 + 2.0*x2;

    vsum3->getElement(i, x3);
    TEST_VALUE_CLOSE(x2, x3, delta);
  }

  // try using a vector that is the wrong size
  gridpack::math::VectorT<TestType> v4(world, local_size-1);
  BOOST_CHECK_THROW(v1.add(v4), gridpack::Exception);
}

BOOST_AUTO_TEST_CASE( add_or_scale_scalar )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> v1(world, local_size), v2(world, local_size);
  TestType x(2.0), offset(2.0), scale(3.0);

  v1.fill(x);
  v1.add(offset);
  v1.print();

  v2.fill(x);
  v2.scale(scale);
  v2.print();

  int lo, hi;
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    TestType y;
    v1.getElement(i, y);
    TEST_VALUE_CLOSE(x+offset, y, delta);
    v2.getElement(i, y);
    TEST_VALUE_CLOSE(x*scale, y, delta);
  }
}

BOOST_AUTO_TEST_CASE( subtract_scalar )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> v1(world, local_size);
  TestType x(1.0);

  v1.fill(x);
  v1.add(-x);
  v1.print();

  int lo, hi;
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    TestType y, zero(0.0);
    v1.getElement(i, y);
    TEST_VALUE_CLOSE(y, zero, delta);
  }
}

BOOST_AUTO_TEST_CASE( elementwise )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> v1(world, local_size);
  gridpack::math::VectorT<TestType> v2(world, local_size);
  gridpack::math::VectorT<TestType> v3(world, local_size);
  TestType x(2.0), y(0.5);

  v1.fill(x);
  v2.fill(x);
  v3.fill(y);

  v1.elementMultiply(v3);
  v2.elementDivide(v3);

  int lo, hi;
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    TestType z1, z2;
    v1.getElement(i, z1);
    v2.getElement(i, z2);
    TEST_VALUE_CLOSE(z1, x*y, delta);
    TEST_VALUE_CLOSE(z2, x/y, delta);
  }
}

BOOST_AUTO_TEST_CASE( elementwise2 )
{
  gridpack::parallel::Communicator world;
  static int three(3);
  gridpack::math::VectorT<TestType> v1(world, three);
  gridpack::math::VectorT<TestType> v2(world, three);
  gridpack::math::VectorT<TestType> v3(world, three);
  gridpack::math::VectorT<TestType> v4(world, three);

  std::vector<TestType> x1(three);
  x1[0] = TEST_VALUE(0.348262, 3.4343 );
  x1[1] = TEST_VALUE(1.50794,  2.76069);
  x1[2] = TEST_VALUE(1.04059,  4.50791);

  std::vector<TestType> x2(three);
  x2[0] = TEST_VALUE(1.08099,  0.0391692);
  x2[1] = TEST_VALUE(1.04585,  0.378384);
  x2[2] = TEST_VALUE(1.08145,  0.249638);

  int lo, hi;
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    v1.setElement(i, x1[i-lo]);
    v2.setElement(i, x2[i-lo]);
    v3.setElement(i, x1[i-lo]*x2[i-lo]);
  }
  v1.ready();
  v2.ready();
  v3.ready();

  v3.print();
  
  v4.equate(v1);
  v4.elementMultiply(v2);

  v4.print();

  for (int i = lo; i < hi; ++i) {
    TestType x, y;
    v3.getElement(i, x);
    v4.getElement(i, y);
    TEST_VALUE_CLOSE(y, x, delta);
  }

  // it should be commutative

  v4.equate(v2);
  v4.elementMultiply(v1);
  v4.print();

  for (int i = lo; i < hi; ++i) {
    TestType x, y;
    v3.getElement(i, x);
    v4.getElement(i, y);
    TEST_VALUE_CLOSE(y, x, delta);
  }
}

BOOST_AUTO_TEST_CASE( exp )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> v1(world, local_size);
  TestType x(7.3);

  v1.fill(x);
  v1.exp();

  int lo, hi;
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    TestType y;
    v1.getElement(i, y);
    TEST_VALUE_CLOSE(log(y), x, delta);
  }
}  


BOOST_AUTO_TEST_CASE( reciprocal )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> v1(world, local_size);
  TestType x(TEST_VALUE(2.0, 5.0));
  TestType rx(1.0/x);

  v1.fill(x);
  v1.reciprocal();

  int lo, hi;
  v1.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    TestType y;
    v1.getElement(i, y);
    TEST_VALUE_CLOSE(rx, y, delta);
  }
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( VectorIOTest )

BOOST_AUTO_TEST_CASE( print )
{
  gridpack::parallel::Communicator world;
  gridpack::math::VectorT<TestType> v1(world, local_size);
  TestType x(TEST_VALUE(2.0, 1.0));
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
  gridpack::math::VectorT<TestType> v1(world, local_size);
  TestType x(TEST_VALUE(2.0, 1.0));
  v1.fill(x);
  v1.print();

  std::string out;
  if (world.size() > 1) {
    out = "vector_binary_parallel.out";
  } else {
    out = "vector_binary_serial.out";
  }

  v1.saveBinary(out.c_str());

  gridpack::math::VectorT<TestType> v2(world, local_size);
  v2.loadBinary(out.c_str());

  int lo, hi;
  v2.localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    TestType x, y;
    v1.getElement(i, x);
    v2.getElement(i, y);
    TEST_VALUE_CLOSE(y, x, delta);
  }
}
 
BOOST_AUTO_TEST_SUITE_END()


