// -------------------------------------------------------------
/**
 * @file   matrix.cpp
 * @author William A. Perkins
 * @date   2013-05-17 13:00:26 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/assert.hpp>
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
  BOOST_ASSERT(p_matrix_impl);
}

Matrix::~Matrix(void)
{
  
}

} // namespace math
} // namespace gridpack
