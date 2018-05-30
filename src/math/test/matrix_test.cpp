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
 * @date   2016-12-16 09:35:46 d3g096
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
#include "gridpack/parallel/random.hpp"
#include "math.hpp"
#include "matrix.hpp"
#include "gridpack/utilities/exception.hpp"

#include "test_main.cpp"

static const int local_size(5);
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
// make_and_fill_test_matrix
// -------------------------------------------------------------
static TestMatrixType *
make_and_fill_test_matrix(const gridpack::parallel::Communicator& comm,
                          const int& bandwidth, int& global_size)
{
  TestMatrixType *A = make_test_matrix(comm, global_size);

  int halfbw(std::max((bandwidth - 1)/2, 0));

  int lo, hi;
  A->localRowRange(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    TestType x(static_cast<double>(i));
    int jmin(std::max(i-halfbw, 0)), jmax(std::min(i+halfbw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->setElement(i, j, x);
    }
  }
  A->ready();
  return A;
}

BOOST_AUTO_TEST_SUITE(MatrixTest)

BOOST_AUTO_TEST_CASE( construction )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr< TestMatrixType > 
    A(make_test_matrix(world, global_size));

  int lo, hi;
  A->localRowRange(lo, hi);

  BOOST_CHECK_GT(hi, lo);

  BOOST_CHECK_EQUAL(A->localRows(), local_size);
  BOOST_CHECK_EQUAL(A->rows(), global_size);
  BOOST_CHECK_EQUAL(A->cols(), global_size);
  gridpack::math::MatrixStorageType tmptype(A->storageType());
  BOOST_CHECK_EQUAL(tmptype, the_storage_type);

  gridpack::parallel::Communicator self(world.divide(1));
  boost::scoped_ptr< TestMatrixType > 
    B(make_test_matrix(self, global_size));

  BOOST_CHECK_EQUAL(B->localRows(), local_size);
  BOOST_CHECK_EQUAL(B->rows(), local_size);
  BOOST_CHECK_EQUAL(B->cols(), local_size);

}

#ifdef TEST_DENSE

static void 
fillAndCheckDense(TestMatrixType *A)
{
  int cols(A->cols());
  int lo, hi;
  A->localRowRange(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    for (int j = 0; j < cols; ++j) {
      TestType x(static_cast<double>(i*(j+1)));
      A->setElement(i, j, x);
    }
  }
  A->ready();

  for (int i = lo; i < hi; ++i) {
    for (int j = 0; j < cols; ++j) {
      TestType x(static_cast<double>(i*(j+1))), y;
      A->getElement(i, j, y);
      TEST_VALUE_CLOSE(x, y, delta);
    }
  }
}

BOOST_AUTO_TEST_CASE(denseConstruction)
{
  gridpack::parallel::Communicator world;
  int lrows(5), rows(lrows*world.size()), cols(3);
  boost::scoped_ptr<TestMatrixType> 
    A(TestMatrixType::createDense(world, rows, cols, 0, 0)),
    B(TestMatrixType::createDense(world, 0, cols, lrows, 0));
  fillAndCheckDense(A.get());
  fillAndCheckDense(B.get());
}

#endif

BOOST_AUTO_TEST_CASE( storage )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr< TestMatrixType > 
    A(make_and_fill_test_matrix(world, 3, global_size));

  gridpack::math::MatrixStorageType tmptype(A->storageType());
  BOOST_CHECK_EQUAL(tmptype, the_storage_type);

  tmptype = the_other_storage_type;
  
  boost::scoped_ptr< TestMatrixType > 
    B(gridpack::math::storageType(*A, tmptype));
  BOOST_CHECK_EQUAL(B->storageType(), tmptype);

  double Anorm(A->norm2());
  double Bnorm(B->norm2());

  // norms should be identical
  BOOST_CHECK_CLOSE(Anorm, Bnorm, delta);

  B.reset();
  A.reset();
}

BOOST_AUTO_TEST_CASE( set_and_get )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr< TestMatrixType > 
    A(make_test_matrix(world, global_size));

  int lo, hi;
  A->localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    TestType x(static_cast<double>(i));
    A->setElement(i, i, x);
  }
  A->ready();

  for (int i = lo; i < hi; ++i) {
    TestType x(static_cast<double>(i));
    TestType y;
    A->getElement(i, i, y);
    TEST_VALUE_CLOSE(x, y, delta);
  }
}

BOOST_AUTO_TEST_CASE( bad_get )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr< TestMatrixType > 
    A(make_test_matrix(world, global_size));

  int lo, hi;
  A->localRowRange(lo, hi);

  TestType x(1.0);

  for (int i = lo; i < hi; ++i) {
    x = static_cast<TestType>(i);
    A->setElement(i, i, x);
  }
  A->ready();

  BOOST_CHECK_THROW( A->getElement(0, global_size, x), gridpack::Exception );
  BOOST_CHECK_THROW( A->getElement(global_size, 0, x), gridpack::Exception );
  
}

BOOST_AUTO_TEST_CASE( bad_set )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr< TestMatrixType > 
    A(make_test_matrix(world, global_size));

  // this does not work for PETSc, apparently indexes are not checked for set?
  // BOOST_CHECK_THROW( A->set_element(0, global_size, x), gridpack::Exception );
  // BOOST_CHECK_THROW( A->set_element(global_size, 0, x), gridpack::Exception );

  // A->setElement(0, global_size, x);
  // A->setElement(global_size, 0, x);
  
}

BOOST_AUTO_TEST_CASE( multiple_set_and_get )
{
  const int nentries(3);
  gridpack::parallel::Communicator world;
  int global_size;
  boost::mpi::all_reduce(world, local_size, global_size, std::plus<int>());

  TestMatrixType 
    A(world, local_size, global_size, the_storage_type);

  int lo, hi;
  A.localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    int n(nentries);
    std::vector<TestType> x(n, static_cast<double>(i));
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
    A.setElements(n, &iidx[startidx], &jidx[startidx], &x[startidx]);
  }
  A.ready();

  for (int i = lo; i < hi; ++i) {
    int n(nentries);
    std::vector<int> iidx(n, i);
    std::vector<int> jidx(n);
    jidx[0] = i-1; jidx[1] = i; jidx[2] = i+1;

    TestType x(static_cast<double>(i));

    // fill w/ bogus values 
    std::vector<TestType> 
      y(n, TEST_VALUE(-1.0, 1.0));
    
    int startidx(0);
    if (i <= 0) {
      startidx = 1;
      n = 2;
    } else if (i >= global_size - 1) {
      n = 2;
    }
    A.getElements(n, &iidx[startidx], &jidx[startidx], &y[startidx]);

    for (int j = startidx; j < n - startidx; ++j) {
      TEST_VALUE_CLOSE(x, y[j], delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( accumulate )
{
  gridpack::parallel::Communicator world;
  int global_size;
  boost::mpi::all_reduce(world, local_size, global_size, std::plus<int>());

  TestMatrixType 
    A(world, local_size, global_size, the_storage_type);

  int lo, hi;
  A.localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    TestType x(static_cast<double>(i));
    A.addElement(i, i, x);
    A.addElement(i, i, x);
  }
  A.ready();

  for (int i = lo; i < hi; ++i) {
    TestType x(static_cast<double>(2*i));
    TestType y;
    A.getElement(i, i, y);
    TEST_VALUE_CLOSE(x, y, delta);
  }
}

BOOST_AUTO_TEST_CASE( local_clone )
{
  gridpack::parallel::Communicator world;
  int global_size;
  boost::scoped_ptr< TestMatrixType > 
    A(make_and_fill_test_matrix(world, 3, global_size));
  boost::scoped_ptr< TestMatrixType > 
    B(A->localClone());

  BOOST_CHECK_EQUAL(B->processor_size(), 1);
  BOOST_CHECK_EQUAL(A->rows(), B->rows());
  BOOST_CHECK_EQUAL(A->cols(), B->cols());
  
  int lo, hi;
  A->localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    TestType x;
    TestType y;
    A->getElement(i, i, x);
    B->getElement(i, i, y);
    TEST_VALUE_CLOSE(x, y, delta);
  }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(MatrixOperationsTest)

BOOST_AUTO_TEST_CASE( clone )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr< TestMatrixType > 
    A(make_test_matrix(world, global_size));

  const int bw(1);
  int lo, hi;
  A->localRowRange(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    TestType x(static_cast<double>(i));
    int jmin(std::max(i-bw, 0)), jmax(std::min(i+bw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->setElement(i, j, x);
    }
  }
  A->ready();

  boost::scoped_ptr< TestMatrixType > 
    Aclone(A->clone());

  BOOST_CHECK_EQUAL(A->rows(), Aclone->rows());
  BOOST_CHECK_EQUAL(A->localRows(), Aclone->localRows());
  BOOST_CHECK_EQUAL(A->cols(), Aclone->cols());

  for (int i = lo; i < hi; ++i) {
    TestType x;
    TestType y;
    int jmin(std::max(i-bw, 0)), jmax(std::min(i+bw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->getElement(i, j, x);
      Aclone->getElement(i, j, y);
      TEST_VALUE_CLOSE(x, y, delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( add )
{
  
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr< TestMatrixType > 
    A(make_test_matrix(world, global_size));

  int lo, hi;
  A->localRowRange(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    TestType x(static_cast<double>(i));
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->setElement(i, j, x);
    }
  }
  A->ready();

  boost::scoped_ptr< TestMatrixType > B(A->clone());
  boost::scoped_ptr< TestMatrixType > C(gridpack::math::add(*A, *B));

  B->add(*A);

  for (int i = lo; i < hi; ++i) {
    TestType x;
    TestType y;
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->getElement(i, j, x);
      C->getElement(i, j, y);
      TEST_VALUE_CLOSE(2.0*x, y, delta);
      B->getElement(i, j, y);
      TEST_VALUE_CLOSE(2.0*x, y, delta);
    }
  }

}

BOOST_AUTO_TEST_CASE( identity )
{
  gridpack::parallel::Communicator world;
  int global_size;
  boost::scoped_ptr< TestMatrixType > 
    A(make_test_matrix(world, global_size));

  int lo, hi;
  A->localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    TestType x(1.0);
    A->setElement(i, i, x);
  }
  A->ready();

  A->identity();

  boost::scoped_ptr< TestMatrixType > 
    B(gridpack::math::identity(*A));
  
  boost::scoped_ptr< TestMatrixType > 
    C(make_and_fill_test_matrix(world, 3, global_size));
  C->identity();

  for (int i = lo; i < hi; ++i) {
    TestType x(1.0);
    TestType y;
    A->getElement(i, i, y);
    TEST_VALUE_CLOSE(x, y, delta);
    B->getElement(i, i, y);
    TEST_VALUE_CLOSE(x, y, delta);
    C->getElement(i, i, y);
    TEST_VALUE_CLOSE(x, y, delta);
  }

}

BOOST_AUTO_TEST_CASE( scale )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_and_fill_test_matrix(world, 3, global_size));

  TestType z(TEST_VALUE(2.0, 2.0));
  // A->print();
  A->scale(z);
  // A->print();

  int lo, hi;
  A->localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      TestType 
        x(static_cast<TestType>(i)*z), y;
      A->getElement(i, j, y);
      TEST_VALUE_CLOSE(x, y, delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( Transpose )
{
  BOOST_TEST_MESSAGE("Transpose");
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_test_matrix(world, global_size));

  int bandwidth(3);
  int halfbw(std::max((bandwidth - 1)/2, 0));

  int lo, hi;
  A->localRowRange(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    TestType x(TEST_VALUE(i, i));
    int jmin(std::max(i-halfbw, 0)), jmax(std::min(i+halfbw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->setElement(i, j, x);
    }
  }
  A->ready();
  A->print();

  boost::scoped_ptr<TestMatrixType> B(gridpack::math::transpose(*A));

  B->print();
  B->localRowRange(lo, hi);
  for (int i = lo; i < hi; ++i) {
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      TestType 
        x(TEST_VALUE(j, j)), y;
      B->getElement(i, j, y);
      TEST_VALUE_CLOSE(x, y, delta);
    }
  }
}



BOOST_AUTO_TEST_CASE( GetRowLocal )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_and_fill_test_matrix(world, 3, global_size));

  int nrow(2);                  // rows per processor to extract
  int ncol(A->cols());
  TestMatrixType::IdxType lo, hi, j, i;
  std::vector<TestMatrixType::IdxType> ridx(nrow);
  A->localRowRange(lo, hi);

  // get the first row on each processor; each processor can get a
  // different row

  std::vector<TestType> vals(nrow*ncol);
  
  for (i = 0; i < nrow; ++i) {
    ridx[i] = lo + i + 1;
  }
  A->getRowBlock(nrow, &ridx[0], &vals[0]);

  // std::cout << world.rank() << ": ";
  // std::copy(vals.begin(), vals.end(), 
  //           std::ostream_iterator<TestType>(std::cout, ", "));
  // std::cout << std::endl;
  ncol = A->cols();
  for (i = 0; i < nrow; ++i) {
    for (j = 0; j < ncol; ++j) {
      TestType x;
      A->getElement(lo + i + 1, j, x);
      TEST_VALUE_CLOSE(vals[j + i*ncol], x, delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( GetRow )
{
  int global_size;
  gridpack::parallel::Communicator world;
  int me(world.rank()), nproc(world.size());
  boost::scoped_ptr<TestMatrixType> 
    A(make_and_fill_test_matrix(world, 3, global_size));

  int nrow(1);
  int ncol(A->cols());
  TestMatrixType::IdxType lo, hi, j, i;
  std::vector<TestMatrixType::IdxType> lridx(nproc, 0), ridx(nproc, 0);
  A->localRowRange(lo, hi);
  lridx[me] = lo;

  boost::mpi::all_reduce(world, &lridx[0], nproc, &ridx[0],
                         std::plus<TestMatrixType::IdxType>());

  // shift the row indices around so each process gets a row from the
  // next

  lo = ridx[0];
  for (i = 0; i < nproc-1; ++i) {
    ridx[i] = ridx[i+1];
  }
  ridx[nproc-1] = lo;

  std::vector<TestType> vals(nrow*ncol);
  A->getRowBlock(nrow, &ridx[me], &vals[0]);
  std::cout << world.rank() << ": ";
  std::copy(vals.begin(), vals.end(), 
            std::ostream_iterator<TestType>(std::cout, ", "));
  std::cout << std::endl;
}

  

BOOST_AUTO_TEST_CASE( ColumnOps )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_and_fill_test_matrix(world, 3, global_size));
  int icolumn(global_size/2);

  boost::scoped_ptr< TestVectorType >  
    cvector(gridpack::math::column(*A, icolumn));
  
  int lo, hi;
  cvector->localIndexRange(lo, hi);

  for (int i = -1; i <= 1; ++i) {
    int idx(icolumn+i);
    if (lo <= idx && idx < hi) {
      TestType 
        x(static_cast<TestType>(idx));
      TestType y;
      cvector->getElement(idx, y);
      TEST_VALUE_CLOSE(x, y, delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( ColumnDiagonalOps )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr< TestMatrixType > 
    A(make_and_fill_test_matrix(world, 3, global_size));
  int icolumn(global_size/2);

  boost::scoped_ptr< TestVectorType >  
    cvector(gridpack::math::column(*A, icolumn)),
    dvector(gridpack::math::diagonal(*A));

  int lo, hi;
  cvector->localIndexRange(lo, hi);

  for (int i = -1; i <= 1; ++i) {
    int idx(icolumn+i);
    if (lo <= idx && idx < hi) {
      TestType x(static_cast<TestType>(idx));
      TestType y;
      cvector->getElement(idx, y);
      TEST_VALUE_CLOSE(x, y, delta);
      if (idx == icolumn) {
        dvector->getElement(idx, y);
        TEST_VALUE_CLOSE(x, y, delta);
      }
    }
    (cvector->communicator()).barrier();
  }

  boost::scoped_ptr< TestMatrixType > 
    B(gridpack::math::diagonal(*dvector, the_storage_type));
  dvector->print();
  B->print();

  // norms of the diagonal matrix and original vector should be very
  // close
  double Bnorm(B->norm2());
  double vnorm(dvector->norm2());
  BOOST_CHECK_CLOSE(Bnorm, vnorm, delta);

  // make the diagonal matrix back into a vector and see that it has
  // not changed

  boost::scoped_ptr<TestVectorType>  
    bvector(gridpack::math::diagonal(*B));

  bvector->scale(-1.0);
  bvector->add(*dvector);
  vnorm = bvector->norm2();
  
  // norm should be really really small
  BOOST_CHECK(vnorm < delta*delta);
}

BOOST_AUTO_TEST_CASE( AddDiagonal )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_test_matrix(world, global_size));
  A->identity();

  A->print();

  boost::scoped_ptr<TestVectorType>  
    v(new TestVectorType(A->communicator(), A->localRows()));
  TestType z = TEST_VALUE(1.0, 1.0);
  v->fill(z);

  A->addDiagonalVector(*v);
  A->print();

  boost::scoped_ptr<TestVectorType> d(diagonal(*A));
  d->print();

  v->add(1.0);
  double vnorm(v->norm1());
  double dnorm(d->norm1());

  BOOST_CHECK_CLOSE(vnorm, dnorm, delta);

  vnorm = v->norm2();
  dnorm = A->norm2();

  BOOST_CHECK_CLOSE(vnorm, dnorm, delta);
}

BOOST_AUTO_TEST_CASE( AddDiagonal2 )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_test_matrix(world, global_size));
  A->identity();

  A->print();

  TestType z = TEST_VALUE(1.0, 1.0);

  A->addDiagonal(z);
  A->print();

  boost::scoped_ptr<TestVectorType> d(diagonal(*A));
  d->print();

  double vnorm = d->norm2();
  double dnorm = A->norm2();

  BOOST_CHECK_CLOSE(vnorm, dnorm, delta);
}
  
BOOST_AUTO_TEST_CASE( MatrixVectorMultiply )
{
  static const int bandwidth(3);
  static const TestType scale(2.0);
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_and_fill_test_matrix(world, bandwidth, global_size)),
    T(transpose(*A));

  // A->print();
  // T->print();

  boost::scoped_ptr<TestVectorType> 
    xvector(new TestVectorType(A->communicator(), A->localRows())),
    yvector, zvector;

  xvector->fill(scale);
  yvector.reset(multiply(*A, *xvector));
  zvector.reset(transposeMultiply(*T, *xvector));
            
  int lo, hi;
  xvector->localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    int bw(bandwidth);
    if (i == 0 || i == global_size - 1) bw--;
    TestType x(static_cast<TestType>(i*bw)*scale), y, z;
    yvector->getElement(i, y);
    TEST_VALUE_CLOSE(x, y, delta);

    zvector->getElement(i, z);
    TEST_VALUE_CLOSE(x, z, delta);
  }
}



BOOST_AUTO_TEST_CASE( MultiplyDiagonalTest )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_and_fill_test_matrix(world, 3, global_size));
  boost::scoped_ptr<TestVectorType>
    dscale(new TestVectorType(A->communicator(), A->localRows()));
  TestType z(2.0);
  dscale->fill(z);

  A->multiplyDiagonal(*dscale);

  std::cout << "From MultiplyDiagonalTest" << std::endl;
  // dscale->print();
  // A->print();

  int lo, hi;
  A->localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    TestType x(static_cast<TestType>(i)*z), y;
    A->getElement(i, i, y);
    TEST_VALUE_CLOSE(x, y, delta);
  }
}

BOOST_AUTO_TEST_CASE( MultiplyIdentity )
{
  static const int bandwidth(3);
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_and_fill_test_matrix(world, bandwidth, global_size)),
    B(new TestMatrixType(A->communicator(), A->localRows(), A->localCols(), 
                                 gridpack::math::Sparse));
  B->identity();
  boost::scoped_ptr<TestMatrixType> 
    C(gridpack::math::multiply(*A, *B));

  int lo, hi;
  A->localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      TestType x, y;
      A->getElement(i, j, x);
      C->getElement(i, j, y);
      TEST_VALUE_CLOSE(x, y, delta);
    }
  }
}

static void
testMatrixMultiply(TestMatrixType *A,
                   TestMatrixType *B)
{
  gridpack::parallel::Communicator comm(A->communicator());

  static const TestType avalues[] =
    { 1.0,  2.0,  6.0,
      3.0,  4.0,  5.0 };
  static const TestType bvalues[] =
    {  7.0, -1.0,
       0.0,  1.0,
      -3.0,  4.0 };
  // This should be the result
  // static const TestType cvalues[] =
  //   { -11.0, 25.0,
  //       6.0, 21.0 }; 

  std::vector<int> iidx(2*3) , jidx(2*3);
  int k(0);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      iidx[k] = i + comm.rank()*2;
      jidx[k] = j + comm.rank()*3;
      k++;
    }
  }
  A->setElements(2*3, &iidx[0], &jidx[0], &avalues[0]);
  A->ready();
  // A->print();

  k = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      iidx[k] = i + comm.rank()*3;
      jidx[k] = j + comm.rank()*2;
      k++;
    }
  }

  B->setElements(3*2, &iidx[0], &jidx[0], &bvalues[0]);
  B->ready();
  // B->print();

  try {
    boost::scoped_ptr<TestMatrixType> C;
    C.reset(gridpack::math::multiply(*A, *B));
    C->print();
    
    BOOST_CHECK_EQUAL(C->rows(), 2*comm.size());
    BOOST_CHECK_EQUAL(C->cols(), 2*comm.size());
    
    TestMatrixType::TheType x, y;
    C->getElement(2*comm.rank(), 2*comm.rank(), x);
    y = -11.0;
    TEST_VALUE_CLOSE(x, y, delta);
    TEST_VALUE_CLOSE(x, y, delta);
    
    C->getElement(2*comm.rank()+1, 2*comm.rank()+1, x);
    y = 21.0;
    TEST_VALUE_CLOSE(x, y, delta);
    TEST_VALUE_CLOSE(x, y, delta);
  } catch (const gridpack::Exception& e) {
    BOOST_ERROR("Matrix-matrix multiply failed with an exception");
  }
}

BOOST_AUTO_TEST_CASE ( MatrixMatrixMultiplySame )
{
  gridpack::parallel::Communicator world;

  boost::scoped_ptr<TestMatrixType> 
    A(new TestMatrixType(world, 2, 3, the_storage_type)),
    B(new TestMatrixType(world, 3, 2, the_storage_type));
  testMatrixMultiply(A.get(), B.get());
}

BOOST_AUTO_TEST_CASE ( MatrixMatrixMultiplyDifferent )
{
  gridpack::parallel::Communicator world;

  boost::scoped_ptr<TestMatrixType> 
    A(new TestMatrixType(world, 2, 3, the_storage_type)),
    B(new TestMatrixType(world, 3, 2, the_other_storage_type));
  testMatrixMultiply(A.get(), B.get());
}

BOOST_AUTO_TEST_CASE( NonSquareTranspose )
{
  gridpack::parallel::Communicator world;

  static const TestType avalues[] =
    { 1.0,  2.0,  3.0,
      4.0,  5.0,  6.0 };
  boost::scoped_ptr<TestMatrixType> 
    A(new TestMatrixType(world, 2, 3, the_storage_type));

  std::vector<int> iidx(2*3) , jidx(2*3);
  int k(0);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      iidx[k] = i + world.rank()*2;
      jidx[k] = j + world.rank()*3;
      k++;
    }
  }
  A->setElements(2*3, &iidx[0], &jidx[0], &avalues[0]);
  A->ready();
  A->print();

  boost::scoped_ptr<TestMatrixType> 
    B(gridpack::math::transpose(*A));
  B->print();

  // FIXME: check B contents

  boost::scoped_ptr<TestMatrixType> 
    C(new TestMatrixType(world, 3, 2, the_storage_type));
  transpose(*A, *C);
  C->print();

  // FIXME: check C contents

  C.reset(A->clone());
  BOOST_CHECK_THROW(transpose(*A, *C), gridpack::Exception);
}

BOOST_AUTO_TEST_CASE( AnotherNonSquareTranspose )
{

  // A test for transposing this PETSc matrix, which someone had
  // trouble with.
  // row 0: (0, 4.8903)  (1, -10.099)  (2, -0.990099)  (3, -8.16)  (4, -4) 
  // row 1: (0, 0.970297)  (1, -10.099)  (2, -0.990099) 
  // row 2: (0, 3.92)  (3, -8.16)  (4, -4) 
  // row 3: (1, 0)  (2, 1) 
  // row 4: (0, -1.0099)  (1, 18.4222)  (2, 5.1097)  (3, -8.3232)  (4, -4.08) 
  // row 5: (1, 8.3232)  (2, 4.08)  (3, -8.3232)  (4, -4.08) 
  // row 6: (3, 0)  (4, 1) 
  // row 7: (0, -4.08)  (1, -8.3232)  (2, -4.08)  (3, 16.4832)  (4, 8.24) 
  
  gridpack::parallel::Communicator world;

  const int isize(8), jsize(5);
  static const int iidx[] =
    {
      0, 0, 0, 0, 0,
      1, 1, 1,
      2, 2, 2,
      3, 3,
      4, 4, 4, 4, 4,
      5, 5, 5, 5,
      6, 6,
      7, 7, 7, 7, 7
    };
  static const int jidx[] =
    {
      0, 1, 2, 3, 4,
      0, 1, 2,
      0, 3, 4,
      1, 2, 
      0, 1, 2, 3, 4,
      1, 2, 3, 4, 
      3, 4,
      0, 1, 2, 3, 4
    };
      
  static const TestType avalues[] =
    { 
      4.8903, -10.099, -0.990099, -8.16, -4,
      0.970297, -10.099, -0.990099, 
      3.92, -8.16, -4,
      0.0, 1, 
      -1.0099, 18.4222, 5.1097, -8.3232, -4.08,
      8.3232, 4.08, -8.3232, -4.08,
      0, 1, 
      -4.08, -8.3232, -4.08, 16.4832, 8.24
    };
  static const int n(29);

  size_t me(world.rank());

  boost::scoped_ptr<TestMatrixType> A;
  switch (the_storage_type) {
  case gridpack::math::Dense:
    A.reset(new TestMatrixType(world, isize, jsize, the_storage_type));
    break;
  case gridpack::math::Sparse:
    A.reset(new TestMatrixType(world, isize, jsize, jsize*3));
    break;
  default:
    throw gridpack::Exception("Unknown MatrixT<TestType> storage type");
  }
  int lo, hi;
  A->localRowRange(lo, hi);

  BOOST_CHECK_EQUAL(lo, isize*me);
  
  for (int k = 0; k < n; ++k) {
    A->setElement(iidx[k] + isize*me, 
                  jidx[k] + jsize*me, 
                  avalues[k]);
  }
  A->ready();
  A->print();

  boost::scoped_ptr<TestMatrixType> 
    B(gridpack::math::transpose(*A));
  B->print();

  // FIXME: check B contents

  boost::scoped_ptr<TestMatrixType> C;
  switch (the_storage_type) {
  case gridpack::math::Dense:
    C.reset(new TestMatrixType(world, jsize, isize, the_storage_type));
    break;
  case gridpack::math::Sparse:
    C.reset(new TestMatrixType(world, jsize, isize, isize*2));
    break;
  default:
    throw gridpack::Exception("Unknown Matrix storage type");
  }

  transpose(*A, *C);
  C->print();

  // FIXME: check C contents

  C.reset(A->clone());
  BOOST_CHECK_THROW(transpose(*A, *C), gridpack::Exception);
}


BOOST_AUTO_TEST_CASE( ComplexOperations )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_test_matrix(world, global_size));

  int lo, hi;
  A->localRowRange(lo, hi);
 
  int halfbw(1);

  for (int i = lo; i < hi; ++i) {
    TestType 
      x(TEST_VALUE(static_cast<double>(i), static_cast<double>(i)));
    int jmin(std::max(i-halfbw, 0)), jmax(std::min(i+halfbw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->setElement(i, j, x);
    }
  }
  A->ready();

  boost::scoped_ptr<TestMatrixType> 
    Areal(real(*A)),
    Aimag(imaginary(*A)),
    Aconj(conjugate(*A));

  for (int i = lo; i < hi; ++i) {
    int jmin(std::max(i-halfbw, 0)), jmax(std::min(i+halfbw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      TestType y;
      gridpack::ComplexType a, areal, aimag, aconj; 
      A->getElement(i, j, y); a = y;
      Areal->getElement(i, j, y); areal = y;
      Aimag->getElement(i, j, y); aimag = y;
      Aconj->getElement(i, j, y); aconj = y;
      BOOST_CHECK_CLOSE(real(a), real(areal), delta);
      BOOST_CHECK_CLOSE(imag(a), real(aimag), delta);
      BOOST_CHECK_CLOSE(imag(a), -imag(aconj), delta);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(MatrixIOTest)

BOOST_AUTO_TEST_CASE( print)
{
  static const int bandwidth(3);
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_and_fill_test_matrix(world, bandwidth, global_size));

  A->print();

  std::string out(print_prefix);
  if (A->processor_size() > 1) {
    out += "matrix_parallel.out";
  } else {
    out += "matrix_serial.out";
  }
  A->print(out.c_str());

  out = print_prefix;
  if (A->processor_size() > 1) {
    out += "matrix_parallel.mat";
  } else {
    out += "matrix_serial.mat";
  }
  A->save(out.c_str());
}


BOOST_AUTO_TEST_CASE( load_save )
{
  static const int bandwidth(3);
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<TestMatrixType> 
    A(make_and_fill_test_matrix(world, bandwidth, global_size));

  A->print();

  std::string out(print_prefix);
  if (A->processor_size() > 1) {
    out += "binary_matrix_parallel.out";
  } else {
    out += "binary_matrix_serial.out";
  }
  A->saveBinary(out.c_str());

  boost::scoped_ptr<TestMatrixType> 
    B(new TestMatrixType(A->communicator(), 
                         A->localRows(), A->localCols(),
                         the_storage_type));
  B->loadBinary(out.c_str());

  B->scale(-1.0);
  B->add(*A);

  BOOST_CHECK_CLOSE(B->norm2(), 0.0, delta);

}

BOOST_AUTO_TEST_SUITE_END()



