// -------------------------------------------------------------
// file: petsc_matrix_operations.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April 17, 2013 by William A. Perkins
// Last Change: 2013-06-05 09:50:13 d3g096
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
// transpose
// -------------------------------------------------------------
Matrix *
transpose(const Matrix& A)
{
  Matrix *result = A.clone();
  Mat *pA = PETScMatrix(*result);
  PetscErrorCode ierr(0);
  try {
    PetscScalar one(1.0);
    ierr = MatTranspose(*pA, MAT_REUSE_MATRIX, pA); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }

  return result;
}


} // namespace math
} // namespace gridpack
