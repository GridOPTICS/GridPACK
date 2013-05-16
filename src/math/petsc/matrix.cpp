/**
 * @file   matrix.cpp
 * @author William A. Perkins
 * @date   2013-05-16 11:18:43 d3g096
 * 
 * @brief  PETSc specific part of Matrix
 * 
 * 
 */

#include "boost/assert.hpp"
#include "gridpack/math/matrix.hpp"
#include "gridpack/math/petsc/petsc_matrix_implementation.hpp"


namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class Matrix
// -------------------------------------------------------------

// -------------------------------------------------------------
// Matrix:: constructors / destructor
// -------------------------------------------------------------
Matrix::Matrix(const parallel::Communicator& comm,
               const int& rows,
               const int& cols,
               const StorageType& storage_type)
  : parallel::Distributed(comm), utility::Uncopyable(),
    p_matrix_impl()
{
  switch (storage_type) {
  case Sparse:
    p_matrix_impl.reset(new PETScMatrixImplementation(this->communicator(),
                                                      rows, cols, false));
    break;
  case Dense:
    p_matrix_impl.reset(new PETScMatrixImplementation(this->communicator(),
                                                      rows, cols, true));
    break;
  default:
    BOOST_ASSERT(false);
  }
  BOOST_ASSERT(p_matrix_impl);
}


} // namespace math
} // namespace gridpack
