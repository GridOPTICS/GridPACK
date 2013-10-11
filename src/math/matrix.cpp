// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   matrix.cpp
 * @author William A. Perkins
 * @date   2013-10-10 07:33:09 d3g096
 * 
 * @brief  Generic part of Matrix implementation
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/assert.hpp>
#include "matrix.hpp"
#include "gridpack/utilities/exception.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
// Matrix constructor
// -------------------------------------------------------------
Matrix::Matrix(MatrixImplementation *impl)
  : parallel::WrappedDistributed(impl), 
    utility::Uncopyable(),
    p_matrix_impl(impl)
{
  BOOST_ASSERT(p_matrix_impl);
}

Matrix::~Matrix(void)
{
  
}

// -------------------------------------------------------------
// Matrix::p_check_compatible
// -------------------------------------------------------------
void
Matrix::p_check_compatible(const Matrix& A) const
{
  if (this->communicator() != A.communicator()) {
    throw gridpack::Exception("incompatible: communicators do not match");
  }

  if ((this->rows() != A.rows()) || (this->cols() != A.cols())) {
    throw gridpack::Exception("incompatible: sizes do not match");
  }
  return;
}

// -------------------------------------------------------------
// Matrix operations
// -------------------------------------------------------------

// -------------------------------------------------------------
// add
// -------------------------------------------------------------
void
add(const Matrix& A, const Matrix& B, Matrix& result) 
{
  result.equate(A);
  result.add(B);
}


Matrix *
add(const Matrix& A, const Matrix& B) 
{
  Matrix *result = A.clone();
  add(A, B, *result);
  return result;
}


// -------------------------------------------------------------
// column
// -------------------------------------------------------------
Vector *
column(const Matrix& A, const int& cidx)
{
  Vector *colv(new Vector(A.communicator(), A.localRows()));
  column(A, cidx, *colv);
  return colv;
}

// -------------------------------------------------------------
// diagonal
// -------------------------------------------------------------
Vector *
diagonal(const Matrix& A)
{
  Vector *colv(new Vector(A.communicator(), A.localRows()));
  diagonal(A, *colv);
  return colv;
}

// -------------------------------------------------------------
// multiply
// -------------------------------------------------------------
Vector *
multiply(const Matrix& A, const Vector& x)
{
  Vector *result(new Vector(x.communicator(), x.localSize()));
  multiply(A, x, *result);
  return result;
}

// -------------------------------------------------------------
// identity
// -------------------------------------------------------------
Matrix *
identity(const Matrix& A)
{
  Matrix *result(A.clone());
  result->identity();
  return result;
}


} // namespace math
} // namespace gridpack
