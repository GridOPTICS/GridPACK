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
 * @date   2014-10-21 13:48:22 d3g096
 * 
 * @brief  Unit tests for Matrix
 * 
 * @test
 */
// -------------------------------------------------------------

#include <iostream>
#include <boost/assert.hpp>
#include <boost/mpi/collectives.hpp>
#include "gridpack/parallel/parallel.hpp"
#include "math.hpp"
#include "matrix.hpp"
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

static const std::string print_prefix = 
#ifdef TEST_DENSE
  "dense_";
#else
  "sparse_";
#endif

// -------------------------------------------------------------
// make_test_matrix
// -------------------------------------------------------------
static gridpack::math::Matrix *
make_test_matrix(const gridpack::parallel::Communicator& comm,
                 int& global_size)
{
  boost::mpi::all_reduce(comm, local_size, global_size, std::plus<int>());

  gridpack::math::Matrix *A =
    new gridpack::math::Matrix(comm, local_size, local_size, the_storage_type);
  return A;
}

// -------------------------------------------------------------
// make_and_fill_test_matrix
// -------------------------------------------------------------
static gridpack::math::Matrix *
make_and_fill_test_matrix(const gridpack::parallel::Communicator& comm,
                          const int& bandwidth, int& global_size)
{
  gridpack::math::Matrix *A = make_test_matrix(comm, global_size);

  int halfbw(std::max((bandwidth - 1)/2, 0));

  int lo, hi;
  A->localRowRange(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(i));
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
  boost::scoped_ptr<gridpack::math::Matrix> A(make_test_matrix(world, global_size));

  int lo, hi;
  A->localRowRange(lo, hi);

  BOOST_CHECK_GT(hi, lo);

  BOOST_CHECK_EQUAL(A->localRows(), local_size);
  BOOST_CHECK_EQUAL(A->rows(), global_size);
  BOOST_CHECK_EQUAL(A->cols(), global_size);
  gridpack::math::Matrix::StorageType tmptype(A->storageType());
  BOOST_CHECK_EQUAL(tmptype, the_storage_type);

  gridpack::parallel::Communicator self(world.divide(1));
  boost::scoped_ptr<gridpack::math::Matrix> B(make_test_matrix(self, global_size));

  BOOST_CHECK_EQUAL(B->localRows(), local_size);
  BOOST_CHECK_EQUAL(B->rows(), local_size);
  BOOST_CHECK_EQUAL(B->cols(), local_size);

}

BOOST_AUTO_TEST_CASE( storage )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::math::Matrix> 
    A(make_and_fill_test_matrix(world, 3, global_size));

  gridpack::math::Matrix::StorageType tmptype(A->storageType());
  BOOST_CHECK_EQUAL(tmptype, the_storage_type);

  switch (A->storageType()) {
  case (gridpack::math::Matrix::Dense):
    tmptype = gridpack::math::Matrix::Sparse;
    break;
  case (gridpack::math::Matrix::Sparse):
    tmptype = gridpack::math::Matrix::Dense;
    break;
  }
  boost::scoped_ptr<gridpack::math::Matrix> 
    B(gridpack::math::storageType(*A, tmptype));
  BOOST_CHECK_EQUAL(B->storageType(), tmptype);
  // norms should be identical
  BOOST_CHECK_CLOSE(A->norm2(), B->norm2(), delta);
}

BOOST_AUTO_TEST_CASE( set_and_get )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::math::Matrix> A(make_test_matrix(world, global_size));

  int lo, hi;
  A->localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(i));
    A->setElement(i, i, x);
  }
  A->ready();

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(i));
    gridpack::ComplexType y;
    A->getElement(i, i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
  }
}

BOOST_AUTO_TEST_CASE( bad_get )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::math::Matrix> A(make_test_matrix(world, global_size));

  int lo, hi;
  A->localRowRange(lo, hi);

  gridpack::ComplexType x(1.0);

  for (int i = lo; i < hi; ++i) {
    x = static_cast<gridpack::ComplexType>(i);
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
  boost::scoped_ptr<gridpack::math::Matrix> A(make_test_matrix(world, global_size));

  gridpack::ComplexType x(1.0);

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

  gridpack::math::Matrix A(world, local_size, global_size, the_storage_type);

  int lo, hi;
  A.localRowRange(lo, hi);

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
    A.setElements(n, &iidx[startidx], &jidx[startidx], &x[startidx]);
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
    A.getElements(n, &iidx[startidx], &jidx[startidx], &y[startidx]);

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
  A.localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(i));
    A.addElement(i, i, x);
    A.addElement(i, i, x);
  }
  A.ready();

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(2*i));
    gridpack::ComplexType y;
    A.getElement(i, i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
  }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(MatrixOperationsTest)

BOOST_AUTO_TEST_CASE( clone )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::math::Matrix> A(make_test_matrix(world, global_size));

  const int bw(1);
  int lo, hi;
  A->localRowRange(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(i));
    int jmin(std::max(i-bw, 0)), jmax(std::min(i+bw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->setElement(i, j, x);
    }
  }
  A->ready();

  boost::scoped_ptr<gridpack::math::Matrix> Aclone(A->clone());

  BOOST_CHECK_EQUAL(A->rows(), Aclone->rows());
  BOOST_CHECK_EQUAL(A->localRows(), Aclone->localRows());
  BOOST_CHECK_EQUAL(A->cols(), Aclone->cols());

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x;
    gridpack::ComplexType y;
    int jmin(std::max(i-bw, 0)), jmax(std::min(i+bw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->getElement(i, j, x);
      Aclone->getElement(i, j, y);
      BOOST_CHECK_CLOSE(real(x), real(y), delta);
      BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( add )
{
  
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::math::Matrix> A(make_test_matrix(world, global_size));

  int lo, hi;
  A->localRowRange(lo, hi);
 
  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(static_cast<double>(i));
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->setElement(i, j, x);
    }
  }
  A->ready();

  boost::scoped_ptr<gridpack::math::Matrix> B(A->clone());
  boost::scoped_ptr<gridpack::math::Matrix> C(gridpack::math::add(*A, *B));

  B->add(*A);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x;
    gridpack::ComplexType y;
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->getElement(i, j, x);
      C->getElement(i, j, y);
      BOOST_CHECK_CLOSE(2.0*real(x), real(y), delta);
      BOOST_CHECK_CLOSE(2.0*abs(x), abs(y), delta);
      B->getElement(i, j, y);
      BOOST_CHECK_CLOSE(2.0*real(x), real(y), delta);
      BOOST_CHECK_CLOSE(2.0*abs(x), abs(y), delta);
    }
  }

}

BOOST_AUTO_TEST_CASE( identity )
{
  gridpack::parallel::Communicator world;
  int global_size;
  boost::scoped_ptr<gridpack::math::Matrix> A(make_test_matrix(world, global_size));

  int lo, hi;
  A->localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(1.0);
    A->setElement(i, i, x);
  }
  A->ready();

  A->identity();

  boost::scoped_ptr<gridpack::math::Matrix> B(gridpack::math::identity(*A));
  
  boost::scoped_ptr<gridpack::math::Matrix> 
    C(make_and_fill_test_matrix(world, 3, global_size));
  C->identity();

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType x(1.0);
    gridpack::ComplexType y;
    A->getElement(i, i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
    B->getElement(i, i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
    C->getElement(i, i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
  }

}

BOOST_AUTO_TEST_CASE( scale )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::math::Matrix> 
    A(make_and_fill_test_matrix(world, 3, global_size));

  gridpack::ComplexType z(2.0);
  A->scale(z);
  
  int lo, hi;
  A->localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      gridpack::ComplexType 
        x(static_cast<gridpack::ComplexType>(i)*z), y;
      A->getElement(i, j, y);
      BOOST_CHECK_CLOSE(real(x), real(y), delta);
      BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( Transpose )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::math::Matrix> 
    A(make_and_fill_test_matrix(world, 3, global_size));

  boost::scoped_ptr<gridpack::math::Matrix> B(gridpack::math::transpose(*A));

  int lo, hi;
  A->localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      gridpack::ComplexType 
        x(static_cast<gridpack::ComplexType>(j)), y;
      B->getElement(i, j, y);
      BOOST_CHECK_CLOSE(real(x), real(y), delta);
      BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
    }
  }
}

BOOST_AUTO_TEST_CASE( ColumnDiagonalOps )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::math::Matrix> 
    A(make_and_fill_test_matrix(world, 3, global_size));
  int icolumn(global_size/2);

  boost::scoped_ptr<gridpack::math::Vector>  
    cvector(gridpack::math::column(*A, icolumn)),
    dvector(gridpack::math::diagonal(*A));

  int lo, hi;
  cvector->localIndexRange(lo, hi);

  for (int i = -1; i <= 1; ++i) {
    int idx(icolumn+i);
    if (lo <= idx && idx < hi) {
      gridpack::ComplexType 
        x(static_cast<gridpack::ComplexType>(idx));
      gridpack::ComplexType y;
      cvector->getElement(idx, y);
      BOOST_CHECK_CLOSE(real(x), real(y), delta);
      if (idx == icolumn) {
        dvector->getElement(idx, y);
        BOOST_CHECK_CLOSE(real(x), real(y), delta);
      }
    }
    (cvector->communicator()).barrier();
  }

  boost::scoped_ptr<gridpack::math::Matrix> 
    B(gridpack::math::diagonal(*dvector, the_storage_type));
  dvector->print();

  // norms of the diagonal matrix and original vector should be very
  // close
  double Bnorm(B->norm2());
  double vnorm(dvector->norm2());
  BOOST_CHECK_CLOSE(Bnorm, vnorm, delta);

  // make the diagonal matrix back into a vector and see that it has
  // not changed

  boost::scoped_ptr<gridpack::math::Vector>  
    bvector(gridpack::math::diagonal(*B));

  bvector->scale(-1.0);
  bvector->add(*dvector);
  vnorm = abs(bvector->norm2());
  
  // norm should be really really small
  BOOST_CHECK(vnorm < delta*delta);
}

BOOST_AUTO_TEST_CASE( AddDiagonal )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::math::Matrix> 
    A(make_test_matrix(world, global_size));
  A->identity();

  A->print();

  boost::scoped_ptr<gridpack::math::Vector>  
    v(new gridpack::math::Vector(A->communicator(), A->localRows()));
  v->fill(1.0);

  A->addDiagonal(*v);
  A->print();

  boost::scoped_ptr<gridpack::math::Vector>  
    d(diagonal(*A));
  d->print();

  double norm(d->norm1()/static_cast<double>(d->size()));

  BOOST_CHECK_CLOSE(norm, 2.0, delta);

  norm = d->norm2()/sqrt(static_cast<double>(d->size()));

  BOOST_CHECK_CLOSE(norm, 2.0, delta);
}
  
  

BOOST_AUTO_TEST_CASE( MatrixVectorMultiply )
{
  static const int bandwidth(3);
  static const gridpack::ComplexType scale(2.0);
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::math::Matrix> 
    A(make_and_fill_test_matrix(world, bandwidth, global_size)),
    T(transpose(*A));

  A->print();
  T->print();

  boost::scoped_ptr<gridpack::math::Vector>  
    xvector(new gridpack::math::Vector(A->communicator(), A->localRows())),
    yvector, zvector;

  xvector->fill(scale);
  yvector.reset(multiply(*A, *xvector));
  zvector.reset(transposeMultiply(*T, *xvector));
            
  int lo, hi;
  xvector->localIndexRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    int bw(bandwidth);
    if (i == 0 || i == global_size - 1) bw--;
    gridpack::ComplexType 
      x(static_cast<gridpack::ComplexType>(i*bw)*scale), y, z;
    yvector->getElement(i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);

    zvector->getElement(i, z);
    BOOST_CHECK_CLOSE(real(x), real(z), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(z), delta);
  }
}



BOOST_AUTO_TEST_CASE( MultiplyDiagonalTest )
{
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::math::Matrix> 
    A(make_and_fill_test_matrix(world, 3, global_size));
  boost::scoped_ptr<gridpack::math::Vector>
    dscale(new gridpack::math::Vector(A->communicator(), A->localRows()));
  gridpack::ComplexType z(2.0);
  dscale->fill(z);

  A->multiplyDiagonal(*dscale);

  std::cout << "From MultiplyDiagonalTest" << std::endl;
  dscale->print();
  A->print();

  int lo, hi;
  A->localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType 
      x(static_cast<gridpack::ComplexType>(i)*z), y;
    A->getElement(i, i, y);
    BOOST_CHECK_CLOSE(real(x), real(y), delta);
    BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
  }
}

BOOST_AUTO_TEST_CASE( MultiplyIdentity )
{
  static const int bandwidth(3);
  int global_size;
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::math::Matrix> 
    A(make_and_fill_test_matrix(world, bandwidth, global_size)),
    B(new gridpack::math::Matrix(A->communicator(), A->localRows(), A->localCols(), 
                                 gridpack::math::Matrix::Sparse));
  B->identity();
  boost::scoped_ptr<gridpack::math::Matrix> 
    C(gridpack::math::multiply(*A, *B));

  int lo, hi;
  A->localRowRange(lo, hi);

  for (int i = lo; i < hi; ++i) {
    int jmin(std::max(i-1, 0)), jmax(std::min(i+1,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      gridpack::ComplexType x, y;
      A->getElement(i, j, x);
      C->getElement(i, j, y);
      BOOST_CHECK_CLOSE(real(x), real(y), delta);
      BOOST_CHECK_CLOSE(abs(x), abs(y), delta);
    }
  }
}

BOOST_AUTO_TEST_CASE ( MatrixMatrixMultiply )
{
  gridpack::parallel::Communicator world;

  static const gridpack::ComplexType avalues[] =
    { 1.0,  2.0,  6.0,
      3.0,  4.0,  5.0 };
  static const gridpack::ComplexType bvalues[] =
    {  7.0, -1.0,
       0.0,  1.0,
      -3.0,  4.0 };
  static const gridpack::ComplexType cvalues[] =
    { -11.0, 25.0,
        6.0, 21.0 }; 
  boost::scoped_ptr<gridpack::math::Matrix> 
    A(new gridpack::math::Matrix(world, 2, 3)),
    B(new gridpack::math::Matrix(world, 3, 2,
                                 gridpack::math::Matrix::Sparse));

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

  k = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      iidx[k] = i + world.rank()*3;
      jidx[k] = j + world.rank()*2;
      k++;
    }
  }

  B->setElements(3*2, &iidx[0], &jidx[0], &bvalues[0]);
  B->ready();
  B->print();

  boost::scoped_ptr<gridpack::math::Matrix>
    C(gridpack::math::multiply(*A, *B));

  C->print();

  BOOST_CHECK_EQUAL(C->rows(), 2*world.size());
  BOOST_CHECK_EQUAL(C->cols(), 2*world.size());

  gridpack::ComplexType x, y;
  C->getElement(2*world.rank(), 2*world.rank(), x);
  y = -11.0;
  BOOST_CHECK_CLOSE(real(x), real(y), delta);
  BOOST_CHECK_CLOSE(abs(x), abs(y), delta);

  C->getElement(2*world.rank()+1, 2*world.rank()+1, x);
  y = 21.0;
  BOOST_CHECK_CLOSE(real(x), real(y), delta);
  BOOST_CHECK_CLOSE(abs(x), abs(y), delta);

}

BOOST_AUTO_TEST_CASE( NonSquareTranspose )
{
  gridpack::parallel::Communicator world;

  static const gridpack::ComplexType avalues[] =
    { 1.0,  2.0,  3.0,
      4.0,  5.0,  6.0 };
  boost::scoped_ptr<gridpack::math::Matrix> 
    A(new gridpack::math::Matrix(world, 2, 3, the_storage_type));

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

  boost::scoped_ptr<gridpack::math::Matrix> 
    B(gridpack::math::transpose(*A));
  B->print();

  // FIXME: check B contents

  boost::scoped_ptr<gridpack::math::Matrix> 
    C(new gridpack::math::Matrix(world, 3, 2, the_storage_type));
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
      
  static const gridpack::ComplexType avalues[] =
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

  size_t nproc(world.size());
  size_t me(world.rank());

  boost::scoped_ptr<gridpack::math::Matrix> A;
  switch (the_storage_type) {
  case gridpack::math::Matrix::Dense:
    A.reset(new gridpack::math::Matrix(world, isize, jsize, the_storage_type));
    break;
  case gridpack::math::Matrix::Sparse:
    A.reset(new gridpack::math::Matrix(world, isize, jsize, jsize));
    break;
  default:
    throw gridpack::Exception("Unknown Matrix storage type");
  }
  int lo, hi;
  A->localRowRange(lo, hi);
  
  for (int k = 0; k < n; ++k) {
    A->setElement(iidx[k] + isize*me, 
                  jidx[k] + jsize*me, 
                  avalues[k]);
  }
  A->ready();
  A->print();

  boost::scoped_ptr<gridpack::math::Matrix> 
    B(gridpack::math::transpose(*A));
  B->print();

  // FIXME: check B contents

  boost::scoped_ptr<gridpack::math::Matrix> C;
  switch (the_storage_type) {
  case gridpack::math::Matrix::Dense:
    C.reset(new gridpack::math::Matrix(world, jsize, isize, the_storage_type));
    break;
  case gridpack::math::Matrix::Sparse:
    C.reset(new gridpack::math::Matrix(world, jsize, isize, isize));
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
  boost::scoped_ptr<gridpack::math::Matrix> 
    A(make_test_matrix(world, global_size));

  int lo, hi;
  A->localRowRange(lo, hi);
 
  int halfbw(1);

  for (int i = lo; i < hi; ++i) {
    gridpack::ComplexType 
      x(static_cast<double>(i), static_cast<double>(i));
    int jmin(std::max(i-halfbw, 0)), jmax(std::min(i+halfbw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      A->setElement(i, j, x);
    }
  }
  A->ready();

  boost::scoped_ptr<gridpack::math::Matrix> 
    Areal(real(*A)),
    Aimag(imaginary(*A)),
    Aconj(conjugate(*A));

  for (int i = lo; i < hi; ++i) {
    int jmin(std::max(i-halfbw, 0)), jmax(std::min(i+halfbw,global_size-1));
    for (int j = jmin; j <= jmax; ++j) {
      gridpack::ComplexType a, areal, aimag, aconj;      
      A->getElement(i, j, a);
      Areal->getElement(i, j, areal);
      Aimag->getElement(i, j, aimag);
      Aconj->getElement(i, j, aconj);
      BOOST_CHECK_CLOSE(real(a), abs(areal), delta);
      BOOST_CHECK_CLOSE(imag(a), abs(aimag), delta);
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
  boost::scoped_ptr<gridpack::math::Matrix> 
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
  boost::scoped_ptr<gridpack::math::Matrix> 
    A(make_and_fill_test_matrix(world, bandwidth, global_size));

  A->print();

  std::string out(print_prefix);
  if (A->processor_size() > 1) {
    out += "binary_matrix_parallel.out";
  } else {
    out += "binary_matrix_serial.out";
  }
  A->saveBinary(out.c_str());

  boost::scoped_ptr<gridpack::math::Matrix> 
    B(new gridpack::math::Matrix(A->communicator(), 
                                 A->localRows(), A->localCols(),
                                 the_storage_type));
  B->loadBinary(out.c_str());

  B->scale(-1.0);
  B->add(*A);

  BOOST_CHECK_CLOSE(B->norm2(), 0.0, delta);

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
  gridpack::parallel::Communicator world;
  gridpack::math::Initialize();
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  gridpack::math::Finalize();
  return result;
}
