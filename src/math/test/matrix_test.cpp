/**
 * @file   matrix_test.cpp
 * @author William A. Perkins
 * @date   2013-06-11 12:10:05 d3g096
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
#include "gridpack/utilities/exception.hpp"

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

// -------------------------------------------------------------
// make_test_matrix
// -------------------------------------------------------------
static gridpack::math::Matrix *
make_test_matrix(int& global_size)
{
  gridpack::parallel::Communicator world;
  boost::mpi::all_reduce(world, local_size, global_size, std::plus<int>());

  gridpack::math::Matrix *A =
    new gridpack::math::Matrix(world, local_size, global_size, the_storage_type);
  return A;
}

// -------------------------------------------------------------
// make_and_fill_test_matrix
// -------------------------------------------------------------
static gridpack::math::Matrix *
make_and_fill_test_matrix(const int& bandwidth, int& global_size)
{
  gridpack::math::Matrix *A = make_test_matrix(global_size);

  int halfbw(std::max((bandwidth - 1)/2, 0));

  int lo, hi;
  A->local_row_range(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(i));
    int jmin(std::max(i-halfbw, 0)), jmax(std::min(i+halfbw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->set_element(i, j, x);
    }
  }
  A->ready();
  return A;
}

BOOST_AUTO_TEST_SUITE(Matrix)

BOOST_AUTO_TEST_CASE( construction )
{
  int global_size;
  std::auto_ptr<gridpack::math::Matrix> A(make_test_matrix(global_size));

  int lo, hi;
  A->local_row_range(lo, hi);

  BOOST_CHECK_GT(hi, lo);

  BOOST_CHECK_EQUAL(A->local_rows(), local_size);
  BOOST_CHECK_EQUAL(A->rows(), global_size);
  BOOST_CHECK_EQUAL(A->cols(), global_size);
  
}

BOOST_AUTO_TEST_CASE( set_and_get )
{
  int global_size;
  std::auto_ptr<gridpack::math::Matrix> A(make_test_matrix(global_size));

  int lo, hi;
  A->local_row_range(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(i));
    A->set_element(i, i, x);
  }
  A->ready();

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(i));
    gridpack::ComplexType y;
    A->get_element(i, i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
  }
}

BOOST_AUTO_TEST_CASE( bad_get )
{
  int global_size;
  std::auto_ptr<gridpack::math::Matrix> A(make_test_matrix(global_size));

  int lo, hi;
  A->local_row_range(lo, hi);

  gridpack::ComplexType x(1.0);

  for (int i = lo; i < hi; ++i) {
    x = static_cast<gridpack::ComplexType>(i);
    A->set_element(i, i, x);
  }
  A->ready();

  BOOST_CHECK_THROW( A->get_element(0, global_size, x), gridpack::Exception );
  BOOST_CHECK_THROW( A->get_element(global_size, 0, x), gridpack::Exception );
  
}

BOOST_AUTO_TEST_CASE( bad_set )
{
  int global_size;
  std::auto_ptr<gridpack::math::Matrix> A(make_test_matrix(global_size));

  gridpack::ComplexType x(1.0);

  // this does not work for PETSc, apparently indexes are not checked for set?
  // BOOST_CHECK_THROW( A->set_element(0, global_size, x), gridpack::Exception );
  // BOOST_CHECK_THROW( A->set_element(global_size, 0, x), gridpack::Exception );

  // A->set_element(0, global_size, x);
  // A->set_element(global_size, 0, x);
  
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
    std::vector<gridpack::ComplexType> x(n, static_cast<double>(i));
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

    gridpack::ComplexType x(static_cast<double>(i));

    // fill w/ bogus values 
    std::vector<gridpack::ComplexType> 
      y(n, gridpack::ComplexType(-1.0, 1.0));
    
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
    gridpack::ComplexType x(static_cast<double>(i));
    A.add_element(i, i, x);
    A.add_element(i, i, x);
  }
  A.ready();

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(2*i));
    gridpack::ComplexType y;
    A.get_element(i, i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
  }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(MatrixOperations)

BOOST_AUTO_TEST_CASE( clone )
{
  int global_size;
  std::auto_ptr<gridpack::math::Matrix> A(make_test_matrix(global_size));

  const int bw(1);
  int lo, hi;
  A->local_row_range(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(i));
    int jmin(std::max(i-bw, 0)), jmax(std::min(i+bw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->set_element(i, j, x);
    }
  }
  A->ready();

  std::auto_ptr<gridpack::math::Matrix> Aclone(A->clone());

  BOOST_CHECK_EQUAL(A->rows(), Aclone->rows());
  BOOST_CHECK_EQUAL(A->local_rows(), Aclone->local_rows());
  BOOST_CHECK_EQUAL(A->cols(), Aclone->cols());

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x;
    gridpack::ComplexType y;
    int jmin(std::max(i-bw, 0)), jmax(std::min(i+bw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->get_element(i, j, x);
      Aclone->get_element(i, j, y);
      BOOST_CHECK_CLOSE(real(x), real(y), delta);
      BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( add )
{
  
  int global_size;
  std::auto_ptr<gridpack::math::Matrix> A(make_test_matrix(global_size));

  int lo, hi;
  A->local_row_range(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(i));
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->set_element(i, j, x);
    }
  }
  A->ready();

  std::auto_ptr<gridpack::math::Matrix> B(A->clone());
  std::auto_ptr<gridpack::math::Matrix> C(gridpack::math::add(*A, *B));

  B->add(*A);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x;
    gridpack::ComplexType y;
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->get_element(i, j, x);
      C->get_element(i, j, y);
      BOOST_CHECK_CLOSE(2.0*real(x), real(y), delta);
      BOOST_CHECK_CLOSE(2.0*abs(x), abs(y), delta);
      B->get_element(i, j, y);
      BOOST_CHECK_CLOSE(2.0*real(x), real(y), delta);
      BOOST_CHECK_CLOSE(2.0*abs(x), abs(y), delta);
    }
  }

}

BOOST_AUTO_TEST_CASE( identity )
{
  int global_size;
  std::auto_ptr<gridpack::math::Matrix> A(make_test_matrix(global_size));

  int lo, hi;
  A->local_row_range(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(1.0);
    A->set_element(i, i, x);
  }
  A->ready();

  A->identity();

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(1.0);
    gridpack::ComplexType y;
    A->get_element(i, i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
  }
}

BOOST_AUTO_TEST_CASE( scale )
{
  int global_size;
  std::auto_ptr<gridpack::math::Matrix> 
    A(make_and_fill_test_matrix(3, global_size));

  gridpack::ComplexType z(2.0);
  A->scale(z);
  
  int lo, hi;
  A->local_row_range(lo, hi);

  for (int i = lo; i < hi; ++i) {
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      gridpack::ComplexType 
        x(static_cast<gridpack::ComplexType>(i)*z), y;
      A->get_element(i, j, y);
      BOOST_CHECK_CLOSE(real(x), real(y), delta);
      BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( transpose )
{
  int global_size;
  std::auto_ptr<gridpack::math::Matrix> 
    A(make_and_fill_test_matrix(3, global_size));

  std::auto_ptr<gridpack::math::Matrix> B(gridpack::math::transpose(*A));

  int lo, hi;
  A->local_row_range(lo, hi);

  for (int i = lo; i < hi; ++i) {
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      gridpack::ComplexType 
        x(static_cast<gridpack::ComplexType>(j)), y;
      B->get_element(i, j, y);
      BOOST_CHECK_CLOSE(real(x), real(y), delta);
      BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( column_diagonal )
{
  int global_size;
  std::auto_ptr<gridpack::math::Matrix> 
    A(make_and_fill_test_matrix(3, global_size));
  int icolumn(global_size/2);

  std::auto_ptr<gridpack::math::Vector>  
    cvector(gridpack::math::column(*A, icolumn)),
    dvector(gridpack::math::diagonal(*A));

  int lo, hi;
  cvector->local_index_range(lo, hi);

  for (int i = -1; i <= 1; ++i) {
    int idx(icolumn+i);
    if (lo <= idx && idx < hi) {
      gridpack::ComplexType 
        x(static_cast<gridpack::ComplexType>(idx));
      gridpack::ComplexType y;
      cvector->get_element(idx, y);
      BOOST_CHECK_CLOSE(real(x), real(y), delta);
      if (idx == icolumn) {
        dvector->get_element(idx, y);
        BOOST_CHECK_CLOSE(real(x), real(y), delta);
      }
    }
    (cvector->communicator()).barrier();
  }

}

BOOST_AUTO_TEST_CASE( matrix_vector_multiply )
{
  static const int bandwidth(3);
  static const gridpack::ComplexType scale(2.0);
  int global_size;
  std::auto_ptr<gridpack::math::Matrix> 
    A(make_and_fill_test_matrix(bandwidth, global_size));

  std::auto_ptr<gridpack::math::Vector>  
    xvector(new gridpack::math::Vector(A->communicator(), A->local_rows())),
    yvector;

  xvector->fill(scale);
  yvector.reset(multiply(*A, *xvector));
            
  int lo, hi;
  xvector->local_index_range(lo, hi);

  for (int i = lo; i < hi; ++i) {
    int bw(bandwidth);
    if (i == 0 || i == global_size - 1) bw--;
    gridpack::ComplexType 
      x(static_cast<gridpack::ComplexType>(i*bw)*scale), y;
    yvector->get_element(i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
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
