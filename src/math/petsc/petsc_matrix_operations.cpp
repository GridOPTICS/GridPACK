// -------------------------------------------------------------
// file: petsc_matrix_operations.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April 17, 2013 by William A. Perkins
// Last Change: 2013-05-09 08:17:15 d3g096
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include <cassert>Error handling 
#include "matrix.hpp"
#include "petsc_matrix_implementation.hpp"
#include "petsc_matrix_extractor.hpp"
#include "petsc_vector_extractor.hpp"

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

  // FIXME: check A size == B size or throw

  // FIXME: check A parallel env == B parallel env or throw

  // FIXME: get a MPI_Comm from A parallel env
  MPI_Comm comm;

  Mat *pA(NULL), *pB(NULL);
  { 
    PETScMatrixExtractor extract;
    A.accept(extract);
    pA = extract.matrix();
  }
  {
    PETScMatrixExtractor extract;
    B.accept(extract);
    pB = extract.matrix();
  }

  assert(pA != NULL);
  assert(pB != NULL);

  Mat presult;
  ierr = MatDuplicate(*pA, MAT_COPY_VALUES, &presult);
  ierr = MatAXPY(*presult, 1.0, *pB, DIFFERENT_NONZERO_PATTERN);
  
  PETScMatrixImplementation *pMatrix(new PETScMatrixImplementation(A.comm(), presult));

  Matrix *result(new Matrix(pMatrix));
  return result;

}


} // namespace utility
} // namespace gridpack
