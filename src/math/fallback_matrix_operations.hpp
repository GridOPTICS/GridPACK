// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   fallback_matrix_operations.hpp
 * @author William A. Perkins
 * @date   2015-06-16 07:31:12 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _fallback_matrix_operations_hpp_
#define _fallback_matrix_operations_hpp_

#include "matrix.hpp"

namespace gridpack {
namespace math {
namespace fallback {

/// Get the diagonal from a Matrix and put in new Vector (fallback)
/** 
 * @e Collective.
 * 
 * @param A
 * @param x 
 * 
 * @return pointer to new, allocated Vector containing the diagonal elements of @c A
 */
template <typename T, typename I>
void diagonal(const MatrixT<T, I>& A, VectorT<T, I>& d)
{
  I lo, hi;
  A.localRowRange(lo, hi);
  for (I i = lo; i < hi; ++i) {
    T x;
    A.getElement(i, i, x);
    d.setElement(i, x);
  }
  d.ready();
}

/// Make a diagonal Matrix from a Vector
/** 
 * @e Collective.
 *
 * This 
 * 
 * @param x 
 * 
 * @return 
 */
template <typename T, typename I>
MatrixT<T, I> *diagonal(const VectorT<T, I>& x, 
                        const MatrixStorageType& stype = Sparse)
{
  MatrixT<T, I> *result;
  result = new MatrixT<T, I>(x.communicator(), x.localSize(), x.localSize(), stype);
  
  I lo, hi;
  x.localIndexRange(lo, hi);

  for (I i = lo; i < hi; ++i) {
    T v;
    x.getElement(i, v);
    result->setElement(i, i, v);
  }
  result->ready();
  return result;
}



/// Get a column from the Matrix and put in new Vector (fallback)
/** 
 * 
 * 
 * @param A 
 * @param cidx 
 * @param c
 * 
 */
template <typename T, typename I>
void
column(const MatrixT<T, I>& A, const I& cidx, VectorT<T, I>& c)
{
  I lo, hi;
  A.localRowRange(lo, hi);
  for (I i = lo; i < hi; ++i) {
    T x;
    A.getElement(i, cidx, x);
    c.setElement(i, x);
  }
  c.ready();
}

/// Multiply this matrix diagonal by the specified vector (fallback)
/** 
 * @e Collective.
 *
 * This is element by element multiplication
 * 
 * @param x factor by which all diagonal elements in the matrix are multiplied
 */
template <typename T, typename I>
void 
multiplyDiagonal(MatrixT<T, I>& A, const VectorT<T, I>& x)
{
  I lrows(A.localRows());
  std::vector<T> dnew(lrows);
  I lo, hi;
  x.localIndexRange(lo, hi);
  for (I i = lo; i < hi; ++i) {
    T v;
    x.getElement(i, v);
    dnew[i-lo] = v;
    A.getElement(i, i, v);
    dnew[i-lo] *= v;
  }

  for (I i = lo; i < hi; ++i) {
    T v(dnew[i-lo]);
    A.setElement(i, i, v);
  }

  A.ready();
}



/// Add the specified vector to the diagonal of this matrix (fallback)
/** 
 * @c Collective.
 *
 * @param A 
 * @param x 
 */
template <typename T, typename I>
void 
addDiagonal(MatrixT<T, I>& A, const VectorT<T, I>& x)
{
  I lo, hi;
  x.localIndexRange(lo, hi);
  for (I i = lo; i < hi; ++i) {
    T v;
    x.getElement(i, v);
    A.addElement(i, i, v);
  }
  A.ready();
}


// -------------------------------------------------------------
// denseMatrixMultiply
// -------------------------------------------------------------
template <typename T, typename I>
void
denseMatrixMultiply(const MatrixT<T, I>& A, const MatrixT<T, I>& B, MatrixT<T, I>& C)
{
  BOOST_ASSERT(A.rows() == C.rows());
  BOOST_ASSERT(B.cols() == C.cols());
  boost::scoped_ptr< VectorT<T, I> > 
    bc(new VectorT<T, I>(B.communicator(), B.localRows())), 
    rc(new VectorT<T, I>(C.communicator(), C.localRows()));
  int lo, hi;
  rc->localIndexRange(lo, hi);
  for (int j = 0; j < B.cols(); ++j) {
    math::column(B, j, *bc);
    math::multiply(A, *bc, *rc);
    for (int i = lo; i < hi; ++i) {
      typename VectorT<T, I>::TheType v;
      rc->getElement(i, v);
      C.setElement(i, j, v);
    }
  }
  C.ready();
}

template <typename T, typename I>
MatrixT<T, I> *
denseMatrixMultiply(const MatrixT<T, I>& A, const MatrixT<T, I>& B)
{
  BOOST_ASSERT(A.cols() == B.rows());
  BOOST_ASSERT(A.localCols() == B.localRows());
  MatrixT<T, I> *C(MatrixT<T, I>::createDense(A.communicator(),
                                              A.rows(), B.cols(),
                                              A.localRows(), B.localCols()));
  denseMatrixMultiply(A, B, *C);
  return C;
}


} // namespace math
} // namespace gridpack
} // namespace fallback


#endif
