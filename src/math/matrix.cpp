// -------------------------------------------------------------
/**
 * @file   matrix.cpp
 * @author William A. Perkins
 * @date   2013-05-15 13:39:36 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/math/matrix.hpp"

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
  
}

Matrix::~Matrix(void)
{
  
}

} // namespace math
} // namespace gridpack
