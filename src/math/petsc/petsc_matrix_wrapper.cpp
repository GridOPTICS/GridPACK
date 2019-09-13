// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_wrapper.cpp
 * @author William A. Perkins
 * @date   2019-09-13 12:53:47 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created January 28, 2015 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------

#include "mpi.h"
#include "petsc/petsc_exception.hpp"
#include "petsc_matrix_wrapper.hpp"
#include "implementation_visitor.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PetscMatrixWrapper
// -------------------------------------------------------------

// -------------------------------------------------------------
// PETScMatrixImplementation::p_getCommunicator
// -------------------------------------------------------------
parallel::Communicator
PetscMatrixWrapper::getCommunicator(const Mat& m)
{
  MPI_Comm comm;
  PetscObjectGetComm((PetscObject)m, &comm);
  parallel::Communicator result(comm);
  return result;
}

// -------------------------------------------------------------
// PetscMatrixWrapper:: constructors / destructor
// -------------------------------------------------------------
/** 
 * 
 * 
 * @param comm 
 * @param local_rows 
 * @param local_cols 
 * @param dense 
 */
PetscMatrixWrapper::PetscMatrixWrapper(const parallel::Communicator& comm,
                                       const PetscInt& local_rows, const PetscInt& local_cols,
                                       const bool& dense)
  : ImplementationVisitable(),
    p_matrix(), p_matrixWrapped(false), p_destroyWrapped(true)
{
  p_build_matrix(comm, local_rows, local_cols);
  if (dense) {
    p_set_dense_matrix();
  } else {
    p_set_sparse_matrix();
  }
}

PetscMatrixWrapper::PetscMatrixWrapper(const parallel::Communicator& comm,
                                       const PetscInt& local_rows, const PetscInt& local_cols,
                                       const PetscInt& max_nonzero_per_row)
  : ImplementationVisitable(),
    p_matrix(), p_matrixWrapped(false), p_destroyWrapped(true)
{
  p_build_matrix(comm, local_rows, local_cols);
  p_set_sparse_matrix(max_nonzero_per_row);
}

PetscMatrixWrapper::PetscMatrixWrapper(const parallel::Communicator& comm,
                                       const PetscInt& local_rows, const PetscInt& local_cols,
                                       const PetscInt *nonzeros_by_row)
  : ImplementationVisitable(),
    p_matrix(), p_matrixWrapped(false), p_destroyWrapped(true)
{
  p_build_matrix(comm, local_rows, local_cols);
  p_set_sparse_matrix(nonzeros_by_row);
}

PetscMatrixWrapper::PetscMatrixWrapper(Mat& m, const bool& copyMat, const bool& destroyMat)
  : ImplementationVisitable(),
    p_matrix(), p_matrixWrapped(false),
    p_destroyWrapped(p_matrixWrapped ? destroyMat : true)
{
  PetscErrorCode ierr;
  try {

    if (copyMat) {
      ierr = MatDuplicate(m, MAT_COPY_VALUES, &p_matrix); CHKERRXX(ierr);
    } else {
      p_matrix = m;
      p_matrixWrapped = true;
    }

  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }

}

PetscMatrixWrapper::~PetscMatrixWrapper(void)
{
  PetscErrorCode ierr(0);
  if (!p_matrixWrapped || (p_matrixWrapped && p_destroyWrapped)) {
    try  {
      PetscBool ok;
      ierr = PetscInitialized(&ok); CHKERRXX(ierr);
      if (ok) {
        ierr = MatDestroy(&p_matrix); CHKERRXX(ierr);
      }
    } catch (...) {
      // just eat it
    }
  }
}

// -------------------------------------------------------------
// PetscMatrixWrapper::p_build_matrix
// -------------------------------------------------------------
void
PetscMatrixWrapper::p_build_matrix(const parallel::Communicator& comm,
                                   const PetscInt& local_rows, 
                                   const PetscInt& local_cols)
{
  PetscErrorCode ierr(0);
  try {

    // If any ownership arguments are specifed, *all* ownership
    // arguments need to be specified.  AND, if the matrix is square,
    // the local rows and cols need to be the same.

    PetscInt lrows(local_rows), grows(PETSC_DECIDE);
    ierr = PetscSplitOwnership(comm, &lrows, &grows); CHKERRXX(ierr);

    PetscInt lcols(local_cols), gcols(PETSC_DECIDE);
    ierr = PetscSplitOwnership(comm, &lcols, &gcols); CHKERRXX(ierr);
    
    ierr = MatCreate(comm, &p_matrix); CHKERRXX(ierr);
    ierr = MatSetSizes(p_matrix, lrows, lcols, grows, gcols); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscMatrixWrapper::p_set_dense_matrix
// -------------------------------------------------------------
void 
PetscMatrixWrapper::p_set_dense_matrix(void)
{
  PetscErrorCode ierr(0);
  try {
    parallel::Communicator comm(getCommunicator(p_matrix));
    if (comm.size() == 1) {
      ierr = MatSetType(p_matrix, MATSEQDENSE); CHKERRXX(ierr);
      ierr = MatSeqDenseSetPreallocation(p_matrix, PETSC_NULL); CHKERRXX(ierr);
    } else {
      ierr = MatSetType(p_matrix, MATDENSE); CHKERRXX(ierr);
      ierr = MatMPIDenseSetPreallocation(p_matrix, PETSC_NULL); CHKERRXX(ierr);
    }
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscMatrixWrapper::p_set_sparse_matrix
// -------------------------------------------------------------
void 
PetscMatrixWrapper::p_set_sparse_matrix(void)
{
  PetscErrorCode ierr(0);
  try {
    parallel::Communicator comm(getCommunicator(p_matrix));
    if (comm.size() == 1) {
      ierr = MatSetType(p_matrix, MATSEQAIJ); CHKERRXX(ierr);
    } else {
      ierr = MatSetType(p_matrix, MATMPIAIJ); CHKERRXX(ierr);
    }
    ierr = MatSetFromOptions(p_matrix); CHKERRXX(ierr);
    ierr = MatSetUp(p_matrix); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

void 
PetscMatrixWrapper::p_set_sparse_matrix(const PetscInt& max_nz_per_row)
{
  PetscErrorCode ierr(0);
  const PetscInt diagonal_non_zero_guess(max_nz_per_row);
  PetscInt offdiagonal_non_zero_guess(static_cast<PetscInt>(diagonal_non_zero_guess));
  offdiagonal_non_zero_guess = std::max(offdiagonal_non_zero_guess, 10);

  try {
    parallel::Communicator comm(getCommunicator(p_matrix));
    if (comm.size() == 1) {
      ierr = MatSetType(p_matrix, MATSEQAIJ); CHKERRXX(ierr);
      ierr = MatSeqAIJSetPreallocation(p_matrix, 
                                       diagonal_non_zero_guess + offdiagonal_non_zero_guess,
                                       PETSC_NULL); CHKERRXX(ierr);
    } else {
      ierr = MatSetType(p_matrix, MATMPIAIJ); CHKERRXX(ierr);
      ierr = MatMPIAIJSetPreallocation(p_matrix, 
                                       diagonal_non_zero_guess,
                                       PETSC_NULL,
                                       offdiagonal_non_zero_guess, 
                                       PETSC_NULL); CHKERRXX(ierr);
    }
    ierr = MatSetFromOptions(p_matrix); CHKERRXX(ierr);
    ierr = MatSetUp(p_matrix); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

void 
PetscMatrixWrapper::p_set_sparse_matrix(const PetscInt *nz_by_row)
{
  std::vector<PetscInt> diagnz;
  PetscInt lrows(this->localRows());
  diagnz.reserve(lrows);
  std::copy(nz_by_row, nz_by_row+lrows, 
            std::back_inserter(diagnz));

  PetscErrorCode ierr(0);
  try {
    parallel::Communicator comm(getCommunicator(p_matrix));
    if (comm.size() == 1) {
      ierr = MatSetType(p_matrix, MATSEQAIJ); CHKERRXX(ierr);
      ierr = MatSeqAIJSetPreallocation(p_matrix, 
                                       PETSC_DECIDE,
                                       &diagnz[0]); CHKERRXX(ierr);
    } else {
      std::vector<PetscInt> offdiagnz;
      offdiagnz.reserve(lrows);
      std::copy(nz_by_row, nz_by_row+lrows, 
                std::back_inserter(offdiagnz));
      ierr = MatSetType(p_matrix, MATMPIAIJ); CHKERRXX(ierr);
      ierr = MatMPIAIJSetPreallocation(p_matrix, 
                                       PETSC_DECIDE,
                                       &diagnz[0],
                                       PETSC_DECIDE, 
                                       &offdiagnz[0]); CHKERRXX(ierr);
    }
    ierr = MatSetFromOptions(p_matrix); CHKERRXX(ierr);
    ierr = MatSetUp(p_matrix); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscMatrixWrapper::localRowRange
// -------------------------------------------------------------
void
PetscMatrixWrapper::localRowRange(PetscInt& lo, PetscInt& hi) const
{
  PetscErrorCode ierr(0);
  try {
    PetscInt plo, phi;
    ierr = MatGetOwnershipRange(p_matrix, &plo, &phi); CHKERRXX(ierr);
    lo = plo;
    hi = phi;
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscMatrixWrapper::rows
// -------------------------------------------------------------
PetscInt 
PetscMatrixWrapper::rows(void) const
{
  PetscErrorCode ierr(0);
  PetscInt result(0);
  try {
    PetscInt rows;
    ierr = MatGetSize(p_matrix, &rows, PETSC_NULL); CHKERRXX(ierr);
    result = rows;
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
  return result;
}

// -------------------------------------------------------------
// PetscMatrixWrapper::localRows
// -------------------------------------------------------------
PetscInt
PetscMatrixWrapper::localRows(void) const
{
  PetscErrorCode ierr(0);
  PetscInt result(0);
  try {
    PetscInt rows;
    ierr = MatGetLocalSize(p_matrix, &rows, PETSC_NULL); CHKERRXX(ierr);
    result = rows;
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
  return result;
}  

// -------------------------------------------------------------
// PetscMatrixWrapper::cols
// -------------------------------------------------------------
PetscInt
PetscMatrixWrapper::cols(void) const
{
  PetscErrorCode ierr(0);
  PetscInt result(0);
  try {
    PetscInt cols;
    ierr = MatGetSize(p_matrix, PETSC_NULL, &cols); CHKERRXX(ierr);
    result = cols;
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
  return result;
}  

// -------------------------------------------------------------
// PetscMatrixWrapper::localCols
// -------------------------------------------------------------
PetscInt
PetscMatrixWrapper::localCols(void) const
{
  PetscErrorCode ierr(0);
  PetscInt result(0);
  try {
    PetscInt cols;
    ierr = MatGetLocalSize(p_matrix, PETSC_NULL, &cols); CHKERRXX(ierr);
    result = cols;
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
  return result;
}  

// -------------------------------------------------------------
// PetscMatrixWrapper::real
// -------------------------------------------------------------
void
PetscMatrixWrapper::real(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = MatRealPart(p_matrix); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscMatrixWrapper::imaginary
// -------------------------------------------------------------
void
PetscMatrixWrapper::imaginary(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = MatImaginaryPart(p_matrix); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscMatrixWrapper::p_conjugate
// -------------------------------------------------------------
void
PetscMatrixWrapper::conjugate(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = MatConjugate(p_matrix); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PETScMatrixImplementation::p_norm2
// -------------------------------------------------------------
double
PetscMatrixWrapper::norm2(void) const
{
  PetscErrorCode ierr(0);
  double result(0.0);
  try {
    PetscReal norm;
    ierr = MatNorm(p_matrix, NORM_FROBENIUS, &norm); CHKERRXX(ierr);
    result = norm;
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
  return result;
}

// -------------------------------------------------------------
// PetscMatrixWrapper::zero
// -------------------------------------------------------------
void
PetscMatrixWrapper::zero(void)
{
  Mat *mat = getMatrix();
  PetscErrorCode ierr(0);
  try {
    ierr = MatZeroEntries(*mat); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscMatrixWrapper::p_ready
// -------------------------------------------------------------
void
PetscMatrixWrapper::ready(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = MatAssemblyBegin(p_matrix, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
    ierr = MatAssemblyEnd(p_matrix, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
    if (false) {
      MatInfo info;
      ierr = MatGetInfo(p_matrix,MAT_LOCAL,&info);
      parallel::Communicator comm(getCommunicator(p_matrix));
      std::cerr << comm.rank() << ": Matrix::ready(): "
                << "size = (" << this->rows() << "x" << this->cols() << "), "
                << "#assembled = " << info.assemblies << ", "
                << "#mallocs = " << info.mallocs << ", "
                << "#nz allocated = " << info.nz_allocated << ", "
                << "#nz used = " << info.nz_used << ", "
                << std::endl;
    }
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// petsc_print_matrix
// -------------------------------------------------------------
static void
petsc_print_matrix(const Mat mat, const char* filename, PetscViewerFormat format)
{
  PetscErrorCode ierr;
  try {
    PetscViewer viewer;
    MPI_Comm comm = PetscObjectComm((PetscObject)mat);
    int me, nproc;
    ierr = MPI_Comm_rank(comm, &me);
    ierr = MPI_Comm_size(comm, &nproc);
    if (filename != NULL) {
      ierr = PetscViewerASCIIOpen(comm, filename, &viewer); CHKERRXX(ierr);
    } else {
      ierr = PetscViewerASCIIGetStdout(comm, &viewer); CHKERRXX(ierr);
    }
#if PETSC_VERSION_GE(3,7,0)
    ierr = PetscViewerPushFormat(viewer, format);
#else
    ierr = PetscViewerSetFormat(viewer, format);
#endif
    CHKERRXX(ierr);
    PetscInt grow, gcol, lrow, lcol;
    switch (format) {
    case PETSC_VIEWER_DEFAULT:
      ierr = MatGetSize(mat, &grow, &gcol); CHKERRXX(ierr);
      ierr = MatGetLocalSize(mat, &lrow, &lcol); CHKERRXX(ierr);
#if PETSC_VERSION_LT(3,7,0)
      ierr = PetscViewerASCIISynchronizedAllow(viewer, PETSC_TRUE); CHKERRXX(ierr);
#else
      ierr = PetscViewerASCIIPushSynchronized(viewer); CHKERRXX(ierr);
#endif
      ierr = PetscViewerASCIIPrintf(viewer,             
                                    "# Matrix distribution\n"); CHKERRXX(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,             
                                    "# proc   rows     cols\n");  CHKERRXX(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,             
                                    "# ---- -------- --------\n"); CHKERRXX(ierr);
      ierr = PetscViewerASCIISynchronizedPrintf(viewer, "# %4d %8d %8d\n",
                                                me, lrow, lcol);  CHKERRXX(ierr);
      ierr = PetscViewerFlush(viewer); CHKERRXX(ierr);
      ierr = MPI_Barrier(comm);
      ierr = PetscViewerASCIIPrintf(viewer,             
                                    "# ---- -------- --------\n");  CHKERRXX(ierr);
      ierr = PetscViewerASCIIPrintf(viewer, "# %4d %8d %8d\n", 
                                    nproc, grow, gcol);  CHKERRXX(ierr);
#if PETSC_VERSION_GE(3,7,0)
      ierr = PetscViewerASCIIPopSynchronized(viewer); CHKERRXX(ierr);
#endif
    default:
      break;
    }
    ierr = MatView(mat, viewer); CHKERRXX(ierr);
#if PETSC_VERSION_GE(3,7,0)
    ierr = PetscViewerPopFormat(viewer); CHKERRXX(ierr);
#endif
    
    if (filename != NULL) {
      ierr = PetscViewerDestroy(&viewer); CHKERRXX(ierr);
    } else {
      // FIXME: causes a SEGV?
      // ierr = PetscViewerDestroy(&viewer); CHKERRXX(ierr);
    }
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}


// -------------------------------------------------------------
// PetscMatrixWrapper::print
// -------------------------------------------------------------
void
PetscMatrixWrapper::print(const char* filename) const
{
  petsc_print_matrix(p_matrix, filename, PETSC_VIEWER_DEFAULT);
}

// -------------------------------------------------------------
// PetscMatrixWrapper::save
// -------------------------------------------------------------
void
PetscMatrixWrapper::save(const char* filename) const
{
  petsc_print_matrix(p_matrix, filename, PETSC_VIEWER_ASCII_MATLAB);
}

// -------------------------------------------------------------
// PetscMatrixWrapper::loadBinary
// -------------------------------------------------------------
void
PetscMatrixWrapper::loadBinary(const char* filename)
{
  PetscErrorCode ierr;
  try {
    PetscViewer viewer;
    parallel::Communicator comm(getCommunicator(p_matrix));
    ierr = PetscViewerBinaryOpen(comm,
                                 filename,
                                 FILE_MODE_READ,
                                 &viewer); CHKERRXX(ierr);
#if PETSC_VERSION_GE(3,7,0)
    ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_NATIVE);
#else
    ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_NATIVE);
#endif
    CHKERRXX(ierr);    
    ierr = MatLoad(p_matrix, viewer); CHKERRXX(ierr);
#if PETSC_VERSION_GE(3,7,0)
    ierr = PetscViewerPopFormat(viewer); CHKERRXX(ierr);
#endif
    ierr = PetscViewerDestroy(&viewer); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// PetscMatrixWrapper::saveBinary
// -------------------------------------------------------------
void
PetscMatrixWrapper::saveBinary(const char* filename) const
{
  PetscErrorCode ierr;
  try {
    PetscViewer viewer;
    parallel::Communicator comm(getCommunicator(p_matrix));
    ierr = PetscViewerBinaryOpen(comm,
                                 filename,
                                 FILE_MODE_WRITE,
                                 &viewer); CHKERRXX(ierr);
#if PETSC_VERSION_GE(3,7,0)
    ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_NATIVE);
#else
    ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_NATIVE); 
#endif
    CHKERRXX(ierr);
    ierr = MatView(p_matrix, viewer); CHKERRXX(ierr);
#if PETSC_VERSION_GE(3,7,0)
    ierr = PetscViewerPopFormat(viewer); CHKERRXX(ierr);
#endif
    ierr = PetscViewerDestroy(&viewer); CHKERRXX(ierr);
  } catch (const PETSC_EXCEPTION_TYPE& e) {
    throw PETScException(ierr, e);
  }
}


// -------------------------------------------------------------
// PetscMatrixWrapper::p_accept
// -------------------------------------------------------------
void 
PetscMatrixWrapper::p_accept(ImplementationVisitor& visitor)
{
  visitor.visit(*this);
}

void 
PetscMatrixWrapper::p_accept(ConstImplementationVisitor& visitor) const
{
  visitor.visit(*this);
}


} // namespace math
} // namespace gridpack
