/**
 * @file   matrix.cpp
 * @author William A. Perkins
 * @date   2013-06-04 14:32:45 d3g096
 * 
 * @brief  PETSc specific part of Matrix
 * 
 * 
 */

#include "boost/assert.hpp"
#include "gridpack/math/matrix.hpp"
#include "gridpack/math/petsc/petsc_exception.hpp"
#include "gridpack/math/petsc/petsc_matrix_implementation.hpp"
#include "gridpack/math/petsc/petsc_matrix_extractor.hpp"


namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class Matrix
// -------------------------------------------------------------

// -------------------------------------------------------------
// Matrix:: constructors / destructor
// -------------------------------------------------------------
Matrix::Matrix(const parallel::Communicator& comm,
               const int& local_rows,
               const int& cols,
               const StorageType& storage_type)
  : parallel::Distributed(comm), utility::Uncopyable(),
    p_matrix_impl()
{
  switch (storage_type) {
  case Sparse:
    p_matrix_impl.reset(new PETScMatrixImplementation(this->communicator(),
                                                      local_rows, cols, false));
    break;
  case Dense:
    p_matrix_impl.reset(new PETScMatrixImplementation(this->communicator(),
                                                      local_rows, cols, true));
    break;
  default:
    BOOST_ASSERT(false);
  }
  BOOST_ASSERT(p_matrix_impl);
}

// -------------------------------------------------------------
// Matrix::add
// -------------------------------------------------------------
void
Matrix::add(const Matrix& B)
{
  Mat *pA(NULL), *pB(NULL);
  { 
    PETScMatrixExtractor extract;
    this->accept(extract);
    pA = extract.matrix();
  }
  {
    PETScMatrixExtractor extract;
    B.accept(extract);
    pB = extract.matrix();
  }

  PetscErrorCode ierr(0);
  try {
    PetscScalar one(1.0);
    ierr = MatAXPY(*pA, one, *pB, DIFFERENT_NONZERO_PATTERN);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

} // namespace math
} // namespace gridpack
