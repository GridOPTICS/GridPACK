// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   petsc_linear_matrx_solver_impl.cpp
 * @author William A. Perkins
 * @date   2013-10-24 09:07:08 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include "petsc_linear_matrix_solver_impl.hpp"
#include "petsc_matrix_implementation.hpp"
#include "petsc_matrix_extractor.hpp"
#include "petsc_exception.hpp"
#include "petsc/petsc_configuration.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PetscLinearMatrixSolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// static members
// -------------------------------------------------------------
MatOrderingType 
PetscLinearMatrixSolverImplementation::p_supportedOrderingType[] =  {
  MATORDERINGNATURAL,
  MATORDERINGND,
  MATORDERING1WD,
  MATORDERINGRCM,
  MATORDERINGQMD,
  MATORDERINGROWLENGTH,
  MATORDERINGAMD
};

MatSolverPackage 
PetscLinearMatrixSolverImplementation::p_supportedSolverPackage[] = {
  MATSOLVERSUPERLU_DIST,
  MATSOLVERSUPERLU,
  MATSOLVERMUMPS,
  MATSOLVERPETSC
};

// -------------------------------------------------------------
// PetscLinearMatrixSolverImplementation:: constructors / destructor
// -------------------------------------------------------------
PetscLinearMatrixSolverImplementation::PetscLinearMatrixSolverImplementation(const Matrix& A)
  : LinearMatrixSolverImplementation(A),
    p_factored(false),
    p_orderingType(MATORDERINGND),
    p_solverPackage(MATSOLVERSUPERLU_DIST),
    p_factorType(MAT_FACTOR_LU),
    p_fill(5),
    p_iterations(2)
{
  // FIXME: maybe enforce the following: A is square, A uses sparse storage
}

PetscLinearMatrixSolverImplementation::~PetscLinearMatrixSolverImplementation(void)
{
  PetscErrorCode ierr;
  try  {
    PetscBool ok;
    ierr = PetscInitialized(&ok);
    if (ok) {
      ierr = MatDestroy(&p_Fmat);
    }
  } catch (...) {
    // just eat it
  }
}

// -------------------------------------------------------------
// PetscLinearMatrixSolverImplementation::p_configure
// -------------------------------------------------------------
void
PetscLinearMatrixSolverImplementation::p_configure(utility::Configuration::Cursor *props)
{
  std::string prefix(petscProcessOptions(this->communicator(), props));

  PetscErrorCode ierr(0);

  std::string mstr;
  size_t n;
  bool found;

  mstr = props->get("Ordering", MATORDERINGND);
  boost::to_lower(mstr);

  n = sizeof(p_supportedOrderingType)/sizeof(MatOrderingType);
  found = false;
  for (int i = 0; i < n; ++i) {
    if (mstr == p_supportedOrderingType[i]) {
      p_orderingType = p_supportedOrderingType[i];
      found = true;
      break;
    }
  }

  if (!found) {
    std::string msg = 
      boost::str(boost::format("%s PETSc configuration: unrecognized \"Ordering\": \"%s\"") %
                 this->configurationKey() % p_orderingType);
    throw Exception(msg);
  }

  mstr = props->get("Package", MATSOLVERSUPERLU_DIST);
  boost::to_lower(mstr);

  n = sizeof(p_supportedSolverPackage)/sizeof(MatSolverPackage);
  found = false;
  for (int i = 0; i < n; ++i) {
    if (mstr == p_supportedSolverPackage[i]) {
      p_solverPackage = p_supportedSolverPackage[i];
      found = true;
      break;
    }
  }

  if (!found) {
    std::string msg = 
      boost::str(boost::format("%s PETSc PETSc configuration: unrecognized \"Package\": \"%s\"") %
                 this->configurationKey() % p_solverPackage);
    throw Exception(msg);
  }

  // FIXME: I cannot make this test work. Not sure why. It would be
  // nice to be able to find out if the package works before we
  // actually try it.

  // PetscBool supported(PETSC_TRUE);
  // try {
  //   Mat *A(PETScMatrix(*p_A));
  //   ierr  = MatGetFactorAvailable(*A, p_solverPackage, p_factorType, &supported); CHKERRXX(ierr);
  // } catch (const PETSc::Exception& e) {
  //   throw PETScException(ierr, e);
  // }

  // if (!supported) {
  //   std::string msg = 
  //     boost::str(boost::format("%s PETSc configuration: unsupported \"Package\": \"%s\" (not built in PETSc?)") %
  //                              this->configurationKey() % p_solverPackage);
  //   throw Exception(msg);
  // }

  p_iterations = props->get("Iterations", p_iterations);
  if (p_iterations <= 0) {
    std::string msg = 
      boost::str(boost::format("%s PETSc configuration: bad \"Iterations\": %d") %
                               this->configurationKey() % p_iterations);
    throw Exception(msg);
  }    

  p_fill = props->get("Fill", p_fill);
  if (p_fill <= 0) {
    std::string msg = 
      boost::str(boost::format("%s PETSc configuration: bad \"Fill\": %d") %
                               this->configurationKey() % p_fill);
    throw Exception(msg);
  }    
}


// -------------------------------------------------------------
// PetscLinearMatrixSolverImplementation::p_factor
// -------------------------------------------------------------
void
PetscLinearMatrixSolverImplementation::p_factor(void) const
{
  PetscErrorCode ierr(0);
  
  try {
    Mat *A(PETScMatrix(*p_A));
    MatFactorInfo  info;
    IS perm, iperm;

    ierr = MatGetOrdering(*A, MATORDERINGND, &perm, &iperm); CHKERRXX(ierr);
    ierr = MatGetFactor(*A, p_solverPackage, p_factorType, &p_Fmat);CHKERRXX(ierr);
    info.fill = p_fill;

    ierr = MatLUFactorSymbolic(p_Fmat, *A, perm, iperm, &info); CHKERRXX(ierr);

    for (int i = 0; i < p_iterations; ++i) {
      ierr = MatLUFactorNumeric(p_Fmat, *A, &info); CHKERRXX(ierr);
    }

    ierr = ISDestroy(&perm); CHKERRXX(ierr);
    ierr = ISDestroy(&iperm); CHKERRXX(ierr);

  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  p_factored = true;
}

// -------------------------------------------------------------
// PetscLinearMatrixSolverImplementation::p_solve
// -------------------------------------------------------------
Matrix *
PetscLinearMatrixSolverImplementation::p_solve(const Matrix& B) const
{
  PetscErrorCode ierr(0);
  Mat X;

  // FIXME: enforce B is dense, or copy and convert
  const Mat *Bmat(PETScMatrix(B));

  try {
    if (!p_factored) {
      p_factor();
    }
    ierr = MatDuplicate(*Bmat, MAT_DO_NOT_COPY_VALUES, &X); CHKERRXX(ierr);
    ierr = MatMatSolve(p_Fmat, *Bmat, X); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }

  PETScMatrixImplementation *ximpl = 
    new PETScMatrixImplementation(this->communicator(), X);
  Matrix *result = new Matrix(ximpl);

  try {
    ierr = MatDestroy(&X); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }

  return result;
}

} // namespace math
} // namespace gridpack
