// -------------------------------------------------------------
/**
 * @file   matrix.cpp
 * @author William A. Perkins
 * @date   2013-06-05 08:45:50 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/assert.hpp>
#include "gridpack/math/matrix.hpp"
#include "gridpack/utilities/exception.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
// Matrix constructor
// -------------------------------------------------------------
Matrix::Matrix(MatrixImplementation *impl)
  : parallel::Distributed(impl->communicator()), 
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


/** 
 * collective
 * 
 * @param A 
 * @param B 
 * 
 * @return pointer to new Matrix containing A+B
 */
Matrix *
add(const Matrix& A, const Matrix& B) 
{
  Matrix *result = A.clone();
  add(A, B, *result);
  return result;
}



} // namespace math
} // namespace gridpack
