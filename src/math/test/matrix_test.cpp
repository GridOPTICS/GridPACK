/**
 * @file   matrix_test.cpp
 * @author William A. Perkins
 * @date   2013-06-04 14:30:53 d3g096
 * 
 * @brief  Unit tests for Matrix
 * 
 * 
 */

#include <iostream>
#include <boost/assert.hpp>
#include <boost/mpi/collectives.hpp>
#include "gridpack/parallel/parallel.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/math/matrix.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

static const int local_size(5);
static const double delta(0.0001);

static const gridpack::math::Matrix::StorageType the_storage_type =
#ifdef TEST_DENSE
  gridpack::math::Matrix::Dense;
#else
  gridpack::math::Matrix::Sparse;
#endif


BOOST_AUTO_TEST_SUITE(Matrix)

BOOST_AUTO_TEST_CASE( construction )
{
  gridpack::parallel::Communicator world;
  int global_size;
  boost::mpi::all_reduce(world, local_size, global_size, std::plus<int>());

  gridpack::math::Matrix A(world, local_size, global_size, the_storage_type);

  int lo, hi;
  A.local_row_range(lo, hi);

  BOOST_CHECK_GT(hi, lo);

  BOOST_CHECK_EQUAL(A.local_rows(), local_size);
  BOOST_CHECK_EQUAL(A.rows(), global_size);
  BOOST_CHECK_EQUAL(A.cols(), global_size);
  
}

BOOST_AUTO_TEST_CASE( set_and_get )
{
  gridpack::parallel::Communicator world;
  int global_size;
  boost::mpi::all_reduce(world, local_size, global_size, std::plus<int>());

  gridpack::math::Matrix A(world, local_size, global_size, the_storage_type);

  int lo, hi;
  A.local_row_range(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::math::complex_type x(static_cast<double>(i));
    A.set_element(i, i, x);
  }
  A.ready();

  for (int i = lo; i < hi; ++i) {
    gridpack::math::complex_type x(static_cast<double>(i));
    gridpack::math::complex_type y;
    A.get_element(i, i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
  }
}

BOOST_AUTO_TEST_CASE( multiple_set_and_get )
{
  const int nentries(3);
  gridpack::parallel::Communicator world;
  int global_size;
  boost::mpi::all_reduce(world, local_size, global_size, std::plus<int>());

  gridpack::math::Matrix A(world, local_size, global_size, the_storage_type);

  int lo, hi;
  A.local_row_range(lo, hi);

  for (int i = lo; i < hi; ++i) {
    int n(nentries);
    std::vector<gridpack::math::complex_type> x(n, static_cast<double>(i));
    std::vector<int> iidx(n, i);
    std::vector<int> jidx(n);
    jidx[0] = i-1; jidx[1] = i; jidx[2] = i+1;
    
    int startidx(0);
    if (i <= 0) {
      startidx = 1;
      n = 2;
    } else if (i >= global_size - 1) {
      n = 2;
    }
    A.set_elements(n, &iidx[startidx], &jidx[startidx], &x[startidx]);
  }
  A.ready();

  for (int i = lo; i < hi; ++i) {
    int n(nentries);
    std::vector<int> iidx(n, i);
    std::vector<int> jidx(n);
    jidx[0] = i-1; jidx[1] = i; jidx[2] = i+1;

    gridpack::math::complex_type x(static_cast<double>(i));

    // fill w/ bogus values 
    std::vector<gridpack::math::complex_type> 
      y(n, gridpack::math::complex_type(-1.0, 1.0));
    
    int startidx(0);
    if (i <= 0) {
      startidx = 1;
      n = 2;
    } else if (i >= global_size - 1) {
      n = 2;
    }
    A.get_elements(n, &iidx[startidx], &jidx[startidx], &y[startidx]);

    for (int j = startidx; j < n - startidx; ++j) {
      BOOST_CHECK_CLOSE(real(x), real(y[j]), delta);
      BOOST_CHECK_CLOSE(abs(x), abs(y[j]), delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( accumulate )
{
  gridpack::parallel::Communicator world;
  int global_size;
  boost::mpi::all_reduce(world, local_size, global_size, std::plus<int>());

  gridpack::math::Matrix A(world, local_size, global_size, the_storage_type);

  int lo, hi;
  A.local_row_range(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::math::complex_type x(static_cast<double>(i));
    A.add_element(i, i, x);
    A.add_element(i, i, x);
  }
  A.ready();

  for (int i = lo; i < hi; ++i) {
    gridpack::math::complex_type x(static_cast<double>(2*i));
    gridpack::math::complex_type y;
    A.get_element(i, i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
  }
}

BOOST_AUTO_TEST_CASE( clone )
{
  gridpack::parallel::Communicator world;
  int global_size;
  boost::mpi::all_reduce(world, local_size, global_size, std::plus<int>());

  gridpack::math::Matrix A(world, local_size, global_size, the_storage_type);

  int lo, hi;
  A.local_row_range(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    gridpack::math::complex_type x(static_cast<double>(i));
    int jmin(std::max(i, 0)), jmax(std::min(i-1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A.set_element(i, j, x);
    }
  }
  A.ready();

  std::auto_ptr<gridpack::math::Matrix> Aclone(A.clone());

  BOOST_CHECK_EQUAL(A.rows(), Aclone->rows());
  BOOST_CHECK_EQUAL(A.local_rows(), Aclone->local_rows());
  BOOST_CHECK_EQUAL(A.cols(), Aclone->cols());

  for (int i = lo; i < hi; ++i) {
    gridpack::math::complex_type x;
    gridpack::math::complex_type y;
    int jmin(std::max(i, 0)), jmax(std::min(i-1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A.get_element(i, j, x);
      Aclone->get_element(i, j, y);
      BOOST_CHECK_CLOSE(real(x), real(y), delta);
      BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( add )
{
  gridpack::parallel::Communicator world;
  int global_size;
  boost::mpi::all_reduce(world, local_size, global_size, std::plus<int>());

  gridpack::math::Matrix A(world, local_size, global_size, the_storage_type);

  int lo, hi;
  A.local_row_range(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    gridpack::math::complex_type x(static_cast<double>(i));
    int jmin(std::max(i, 0)), jmax(std::min(i-1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A.set_element(i, j, x);
    }
  }
  A.ready();

  std::auto_ptr<gridpack::math::Matrix> B(A.clone());
  std::auto_ptr<gridpack::math::Matrix> C(gridpack::math::add(A, *B));

  B->add(A);

  for (int i = lo; i < hi; ++i) {
    gridpack::math::complex_type x;
    gridpack::math::complex_type y;
    int jmin(std::max(i, 0)), jmax(std::min(i-1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A.get_element(i, j, x);
      C->get_element(i, j, y);
      BOOST_CHECK_CLOSE(2.0*real(x), real(y), delta);
      BOOST_CHECK_CLOSE(2.0*abs(x), abs(y), delta);
      B->get_element(i, j, y);
      BOOST_CHECK_CLOSE(2.0*real(x), real(y), delta);
      BOOST_CHECK_CLOSE(2.0*abs(x), abs(y), delta);
    }
  }

}

BOOST_AUTO_TEST_SUITE_END()


// -------------------------------------------------------------
// init_function
// -------------------------------------------------------------
bool init_function()
{
  switch (the_storage_type) {
  case gridpack::math::Matrix::Dense:
    BOOST_TEST_MESSAGE("Testing Dense Matrices ...");
    break;
  case gridpack::math::Matrix::Sparse:
    BOOST_TEST_MESSAGE("Testing Sparse Matrices ...");
    break;
  default:
    BOOST_ASSERT(false);
  }

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
