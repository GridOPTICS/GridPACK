// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   matrix_test.cpp
 * @author William A. Perkins
 * @date   2018-12-10 13:42:57 d3g096
 * 
 * @brief  Unit tests for Matrix
 * 
 * @test
 */
// -------------------------------------------------------------

#include <iostream>
#include <iterator>
#include <boost/assert.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "gridpack/parallel/random.hpp"
#include "math.hpp"
#include "matrix.hpp"
#include "gridpack/utilities/exception.hpp"

#include "test_main.cpp"

static const int local_size(15);
static const double delta(0.0001);

#ifdef TEST_REAL

typedef gridpack::RealType TestType;

#define TEST_VALUE_CLOSE(x, y, delta) \
  BOOST_CHECK_CLOSE((y), (x), delta);

#define TEST_VALUE(r, i) (r)

#else 

typedef gridpack::ComplexType TestType;
#define TEST_VALUE_CLOSE(x, y, delta) \
  BOOST_CHECK_CLOSE(std::real(y), std::real(x), delta); \
  BOOST_CHECK_CLOSE(std::imag(y), std::imag(x), delta); 

#define TEST_VALUE(r, i) TestType(r,i)

#endif

typedef gridpack::math::MatrixT<TestType> TestMatrixType;
typedef gridpack::math::VectorT<TestType> TestVectorType;

static const gridpack::math::MatrixStorageType the_storage_type =
#ifdef TEST_DENSE
  gridpack::math::Dense;
#else
  gridpack::math::Sparse;
#endif

static const gridpack::math::MatrixStorageType the_other_storage_type =
#ifdef TEST_DENSE
  gridpack::math::Sparse;
#else
  gridpack::math::Dense;
#endif

static const std::string print_prefix = 
#ifdef TEST_DENSE
  "dense_";
#else
  "sparse_";
#endif

// -------------------------------------------------------------
// make_test_matrix
// -------------------------------------------------------------
static TestMatrixType *
make_test_matrix(const gridpack::parallel::Communicator& comm,
                 int& global_size)
{
  boost::mpi::all_reduce(comm, local_size, global_size, std::plus<int>());

  TestMatrixType *A =
    new TestMatrixType(comm, local_size, local_size, the_storage_type);
  return A;
}

// -------------------------------------------------------------
// random_test_value
// -------------------------------------------------------------
static boost::random::mt19937 rng;
static TestType
random_test_value()
{
  boost::random::uniform_real_distribution<> rdist(-1, 1);
  double r(rdist(rng));
  double i(rdist(rng));
  return TEST_VALUE(r, i);
}
 

BOOST_AUTO_TEST_SUITE(MatrixTransposeTest)

BOOST_AUTO_TEST_CASE( TransposeRandom )
{
  BOOST_TEST_MESSAGE("Transpose");
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_test_matrix(world, global_size));

  int bandwidth(9);
  int halfbw(std::max((bandwidth - 1)/2, 0));

  int lo, hi;
  A->localRowRange(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    int jmin(std::max(i-halfbw, 0)), jmax(std::min(i+halfbw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      TestType x(random_test_value());
      A->setElement(i, j, x);
    }
  }
  A->ready();
  A->print();

  boost::scoped_ptr<TestMatrixType> 
    Aloc(A->localClone());

  boost::scoped_ptr<TestMatrixType> B(gridpack::math::transpose(*A));

  B->print();
  B->localRowRange(lo, hi);
  for (int i = lo; i < hi; ++i) {
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      TestType x, y;
      Aloc->getElement(j, i, x);
      B->getElement(i, j, y);
      TEST_VALUE_CLOSE(x, y, delta);
    }
  }
}




BOOST_AUTO_TEST_SUITE_END()



