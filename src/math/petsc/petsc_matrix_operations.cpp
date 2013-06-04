// -------------------------------------------------------------
// file: petsc_matrix_operations.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April 17, 2013 by William A. Perkins
// Last Change: 2013-06-04 14:09:27 d3g096
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include <boost/assert.hpp>
#include "gridpack/math/matrix.hpp"
#include "gridpack/math/petsc/petsc_exception.hpp"
#include "gridpack/math/petsc/petsc_matrix_implementation.hpp"
#include "gridpack/math/petsc/petsc_matrix_extractor.hpp"
#include "gridpack/math/petsc/petsc_vector_implementation.hpp"
#include "gridpack/math/petsc/petsc_vector_extractor.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
// Matrix::add
// -------------------------------------------------------------
/** 
 * collective
 * 
 * @param A 
 * @param B 
 * 
 * @return 
 */
Matrix *
add(const Matrix& A, const Matrix& B) 
{
  PetscErrorCode ierr;
  if (A.communicator() != B.communicator()) {
    throw gridpack::Exception("Matrix::add: communicators do not match");
  }

  if ((A.rows() != B.rows()) || (A.cols() != B.cols())) {
    throw gridpack::Exception("Matrix::add: sizes do not match");
  }

  Matrix *result = A.clone();

  Mat *pA(NULL), *pB(NULL);
  { 
    PETScMatrixExtractor extract;
    result->accept(extract);
    pA = extract.matrix();
  }
  {
    PETScMatrixExtractor extract;
    B.accept(extract);
    pB = extract.matrix();
  }

  BOOST_ASSERT(pA != NULL);
  BOOST_ASSERT(pB != NULL);

  try {
    PetscScalar one(1.0);
    ierr = MatAXPY(*pA, one, *pB, DIFFERENT_NONZERO_PATTERN);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }

  return result;
}


} // namespace utility
} // namespace gridpack
