// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   matrix_transpose_test.cpp
 * @author William A. Perkins
 * @date   2018-12-18 09:31:40 d3g096
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

static const int local_rows(5), local_cols(6);
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
                 int& global_rows, int& global_cols)
{
  int myrows(local_rows);
  int mycols(local_cols);
  if (comm.rank() == 0) {       // create an imbalance)
    --myrows;
    --mycols;
  }
  boost::mpi::all_reduce(comm, myrows, global_rows, std::plus<int>());
  boost::mpi::all_reduce(comm, mycols, global_cols, std::plus<int>());

  TestMatrixType *A =
    new TestMatrixType(comm, myrows, mycols, the_storage_type);
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
  int global_rows, global_cols;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_test_matrix(world, global_rows, global_cols));

  boost::random::uniform_real_distribution<> jdist(0, global_cols - 1);
  
  int lo, hi;
  A->localRowRange(lo, hi);

#ifdef TEST_DENSE
  // fill the entire matrix
  for (int i = lo; i < hi; ++i) {
    for (int j = 0; j < global_cols; ++j) {
      TestType x(random_test_value());
      A->setElement(i, j, x);
    }
  }
#else
  // fill approximately 25% of the matrix
  int nentry(global_cols/4 + 1); 
  for (int i = lo; i < hi; ++i) {
    for (int e = 0; e < nentry; ++e) {
      TestType x(random_test_value());
      int j(jdist(rng));
      A->setElement(i, j, x);
    }
  }
#endif
  A->ready();
  A->print();
  A->save("A.m");

  boost::scoped_ptr<TestMatrixType> 
    Aloc(A->localClone());

  boost::scoped_ptr<TestMatrixType> B(gridpack::math::transpose(*A));

  B->print();
  B->save("B.m");

  // check every element of the transpose
  
  B->localRowRange(lo, hi);
  for (int i = lo; i < hi; ++i) {
    for (int j = 0; j < global_rows; ++j) {
      TestType x, y;
      Aloc->getElement(j, i, x);
      B->getElement(i, j, y);
      TEST_VALUE_CLOSE(x, y, delta);
    }
  }

  // see if the transpose is reversible (transpose of transpose needs
  // to be distributed the same as A
  
  boost::scoped_ptr<TestMatrixType> C(A->clone());
  gridpack::math::transpose(*B, *C);
  C->scale(-1.0);
  boost::scoped_ptr<TestMatrixType> E(gridpack::math::add(*A, *C));

  BOOST_CHECK_CLOSE(E->norm2(), 0.0, delta);
}

BOOST_AUTO_TEST_SUITE_END()



