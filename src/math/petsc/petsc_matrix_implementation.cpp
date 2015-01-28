// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_implementation.cpp
 * @author William A. Perkins
 * @date   2015-01-28 11:01:01 d3g096
 * 
 * @brief  PETSc-specific matrix implementation
 * 
 * 
 */
// -------------------------------------------------------------

#include "matrix.hpp"
#include "implementation_visitor.hpp"
#include "petsc/petsc_matrix_implementation.hpp"
#include "petsc/petsc_exception.hpp"
#include "petsc/petsc_matrix_extractor.hpp"


namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScMatrixImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// PETScMatrixImplementation:: constructors / destructor
// -------------------------------------------------------------
/** 
 * 
 * 
 * @param comm 
 * @param local_rows rows to make local to this processor
 * @param cols 
 * @param dense 
 */
PETScMatrixImplementation::PETScMatrixImplementation(const parallel::Communicator& comm,
                                                     const PETScMatrixImplementation::IdxType& local_rows, 
                                                     const PETScMatrixImplementation::IdxType& local_cols,
                                                     const bool& dense)
  : MatrixImplementation(comm),
    p_mwrap(comm, local_rows, local_cols, dense)
{
}

PETScMatrixImplementation::PETScMatrixImplementation(const parallel::Communicator& comm,
                                                     const PETScMatrixImplementation::IdxType& local_rows, 
                                                     const PETScMatrixImplementation::IdxType& local_cols,
                                                     const PETScMatrixImplementation::IdxType& max_nonzero_per_row)
  : MatrixImplementation(comm),
    p_mwrap(comm, local_rows, local_cols, max_nonzero_per_row)
{
}

PETScMatrixImplementation::PETScMatrixImplementation(const parallel::Communicator& comm,
                                                     const PETScMatrixImplementation::IdxType& local_rows, 
                                                     const PETScMatrixImplementation::IdxType& local_cols,
                                                     const PETScMatrixImplementation::IdxType *nonzero_by_row)
  : MatrixImplementation(comm),
    p_mwrap(comm, local_rows, local_cols, nonzero_by_row)
{
}

PETScMatrixImplementation::PETScMatrixImplementation(Mat& m, const bool& copyMat)
  : MatrixImplementation(PetscMatrixWrapper::getCommunicator(m)),
    p_mwrap(m, copyMat)
{
}

PETScMatrixImplementation::~PETScMatrixImplementation(void)
{
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_setElement
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_setElement(const PETScMatrixImplementation::IdxType& i, 
                                        const PETScMatrixImplementation::IdxType& j, 
                                        const PETScMatrixImplementation::TheType& x)
{
  PetscErrorCode ierr(0);
  try {
    Mat *mat = p_mwrap.getMatrix();
    ierr = MatSetValue(*mat, i, j, x, INSERT_VALUES); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_set_elements
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_setElements(const PETScMatrixImplementation::IdxType& n, 
                                         const PETScMatrixImplementation::IdxType *i, 
                                         const PETScMatrixImplementation::IdxType *j, 
                                         const PETScMatrixImplementation::TheType *x)
{
  // FIXME: There's probably a better way
  for (PETScMatrixImplementation::IdxType k = 0; k < n; k++) {
    this->p_setElement(i[k], j[k], x[k]);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_add_element
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_addElement(const PETScMatrixImplementation::IdxType& i, 
                                        const PETScMatrixImplementation::IdxType& j, 
                                        const PETScMatrixImplementation::TheType& x)
{
  PetscErrorCode ierr(0);
  try {
    Mat *mat = p_mwrap.getMatrix();
    ierr = MatSetValue(*mat, i, j, x, ADD_VALUES); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_add_elements
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_addElements(const PETScMatrixImplementation::IdxType& n, 
                                         const PETScMatrixImplementation::IdxType *i, 
                                         const PETScMatrixImplementation::IdxType *j, 
                                         const PETScMatrixImplementation::TheType *x)
{
  // FIXME: There's probably a better way
  for (PETScMatrixImplementation::IdxType k = 0; k < n; k++) {
    this->p_addElement(i[k], j[k], x[k]);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_get_element
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_getElement(const PETScMatrixImplementation::IdxType& i, 
                                        const PETScMatrixImplementation::IdxType& j, 
                                        PETScMatrixImplementation::TheType& x) const
{
  PetscErrorCode ierr(0);
  try {
    const Mat *mat = p_mwrap.getMatrix();
    static const PETScMatrixImplementation::IdxType one(1);
    PetscInt iidx[one] = { i };
    PetscInt jidx[one] = { j };
    ierr = MatGetValues(*mat, one, &iidx[0], one, &jidx[0], &x); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_get_elements
// -------------------------------------------------------------
void
PETScMatrixImplementation::p_getElements(const PETScMatrixImplementation::IdxType& n,
                                         const PETScMatrixImplementation::IdxType *i, 
                                         const PETScMatrixImplementation::IdxType *j, 
                                         PETScMatrixImplementation::TheType *x) const
{
  // FIXME: There is a better way
  for (int k = 0; k < n; k++) {
    this->p_getElement(i[k], j[k], x[k]);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_clone
// -------------------------------------------------------------
MatrixImplementation *
PETScMatrixImplementation::p_clone(void) const
{
  PETScMatrixImplementation *result =
    new PETScMatrixImplementation(const_cast<Mat&>(*(p_mwrap.getMatrix())), true);
  return result;
}

} // namespace math
} // namespace gridpack
