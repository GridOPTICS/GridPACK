// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   petsc_linear_matrx_solver_impl.hpp
 * @author William A. Perkins
 * @date   2019-12-03 08:18:32 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_linear_matrx_solver_impl_hpp_
#define _petsc_linear_matrx_solver_impl_hpp_

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <petscconf.h>
#include <petscmat.h>
#include "linear_matrix_solver_implementation.hpp"
#include "petsc_configurable.hpp"
#include "petsc_matrix_implementation.hpp"
#include "petsc_matrix_extractor.hpp"
#include "petsc_exception.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PetscLinearMatrixSolverImplementation
// -------------------------------------------------------------
template <typename T, typename I = int>
class PetscLinearMatrixSolverImplementation 
  : public LinearMatrixSolverImplementation<T, I>,
    private PETScConfigurable
{
public:

  typedef typename BaseLinearMatrixSolverInterface<T, I>::MatrixType MatrixType;

  /// Default constructor.
  PetscLinearMatrixSolverImplementation(const MatrixType& A)
    : LinearMatrixSolverImplementation<T, I>(A),
      PETScConfigurable(this->communicator()),
      p_factored(false),
      p_orderingType(MATORDERINGND),
#if defined(PETSC_HAVE_SUPERLU_DIST)
      p_solverPackage(MATSOLVERSUPERLU_DIST),
#elif defined(PETSC_HAVE_MUMPS)
      p_solverPackage(MATSOLVERMUMPS),
#else
      p_solverPackage(MATSOLVERPETSC),
#endif    
      p_factorType(MAT_FACTOR_LU),
      p_fill(5), p_pivot(false)
  {
    // FIXME: maybe enforce the following: A is square, A uses sparse storage
  }

  /// Destructor
  ~PetscLinearMatrixSolverImplementation(void)
  {
    PetscErrorCode ierr(0);
    try  {
      PetscBool ok;
      ierr = PetscInitialized(&ok);
      if (ok && p_factored) {
        ierr = MatDestroy(&p_Fmat);
      }
    } catch (...) {
      // just eat it
    }
  }

protected:

  /// Choose a matrix solver type based on PETSc version
  typedef 
#if PETSC_VERSION_LT(3,9,0)
  MatSolverPackage
#else
  MatSolverType
#endif
  ThePetscMatSolverType;

  /// The underlying PETSc factored coefficient matrix
  mutable Mat p_Fmat;

  /// Is p_Fmat ready?
  mutable bool p_factored;

  /// List of supported matrix ordering
  static MatOrderingType p_supportedOrderingType[];

  /// Number of supported matrix orderings
  static int p_nSupportedOrderingTypes;

  /// PETSc matrix ordering type
  MatOrderingType p_orderingType;

  /// List of supported solver packages
  static ThePetscMatSolverType p_supportedSolverPackage[];

  /// Number of supported solver packages
  static int p_nSupportedSolverPackages;

  /// PETSC solver package for factorization
  ThePetscMatSolverType p_solverPackage;

  /// PETSc factorization method to use
  MatFactorType p_factorType;

  /// Fill levels to use in decomposition
  int p_fill;

  /// Flag to enable pivoting
  bool p_pivot;

  /// Do what is necessary to build this instance
  void p_build(const std::string& option_prefix)
  {
    // Yep. Empty.
  }

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props)
  {
    std::string mstr;
    size_t n;
    bool found;

    if (props) {
      mstr = props->get("Ordering", MATORDERINGND);
    } else {
      mstr = MATORDERINGND;
    }
    boost::to_lower(mstr);

    n = p_nSupportedOrderingTypes;
    found = false;
    for (size_t i = 0; i < n; ++i) {
      if (mstr == p_supportedOrderingType[i]) {
        p_orderingType = p_supportedOrderingType[i];
        found = true;
        break;
      }
    }

    if (!found) {
      std::string msg = 
        boost::str(boost::format("%s PETSc configuration: unrecognized \"Ordering\": \"%s\"") %
                   this->configurationKey() % mstr);
      throw Exception(msg);
    }

    if (props) {
      mstr = props->get("Package", 
#if defined(PETSC_HAVE_SUPERLU_DIST)
                        MATSOLVERSUPERLU_DIST
#elif defined(PETSC_HAVE_MUMPS)
                        MATSOLVERMUMPS
#else
                        MATSOLVERPETSC
#endif
                        );
      boost::to_lower(mstr);

      n = p_nSupportedSolverPackages;
      found = false;
      for (size_t i = 0; i < n; ++i) {
        if (mstr == p_supportedSolverPackage[i]) {
          p_solverPackage = p_supportedSolverPackage[i];
          found = true;
          break;
        }
      }

      if (!found) {
        std::string msg = 
          boost::str(boost::format("%s PETSc PETSc configuration: unrecognized \"Package\": \"%s\"") %
                     this->configurationKey() % mstr);
        throw Exception(msg);
      }
      p_fill = props->get("Fill", p_fill);
      p_pivot = props->get("Pivot", p_pivot);
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
  
    if (p_fill <= 0) {
      std::string msg = 
        boost::str(boost::format("%s PETSc configuration: bad \"Fill\": %d") %
                   this->configurationKey() % p_fill);
      throw Exception(msg);
    }    

    this->build(props);
  }

  /// Factor the coefficient matrix
  void p_factor(void) const
  {
    PetscErrorCode ierr(0);
  
    try {
      Mat *A(PETScMatrix(*LinearMatrixSolverImplementation<T, I>::p_A));
      MatFactorInfo  info;
      IS perm, iperm;

      ierr = MatGetOrdering(*A, p_orderingType, &perm, &iperm); CHKERRXX(ierr);
      ierr = MatGetFactor(*A, p_solverPackage, p_factorType, &p_Fmat);CHKERRXX(ierr);
      info.fill = p_fill;
      info.dtcol = (p_pivot ? 1 : 0);

      ierr = MatLUFactorSymbolic(p_Fmat, *A, perm, iperm, &info); CHKERRXX(ierr);
      ierr = MatLUFactorNumeric(p_Fmat, *A, &info); CHKERRXX(ierr);

      ierr = ISDestroy(&perm); CHKERRXX(ierr);
      ierr = ISDestroy(&iperm); CHKERRXX(ierr);

    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }
    p_factored = true;
  }

  /// Solve w/ the specified RHS Matrix (specialized)
  MatrixType *p_solve(const MatrixType& B) const
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
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }

    PETScMatrixImplementation<T, I> *ximpl = 
      new PETScMatrixImplementation<T, I>(X, true);
    MatrixT<T, I> *result = new MatrixT<T, I>(ximpl);

    try {
      ierr = MatDestroy(&X); CHKERRXX(ierr);
    } catch (const PETSC_EXCEPTION_TYPE& e) {
      throw PETScException(ierr, e);
    }

    return result;
  }

};

template <typename T, typename I>
MatOrderingType
PetscLinearMatrixSolverImplementation<T, I>::p_supportedOrderingType[] = {
  MATORDERINGNATURAL,
  MATORDERINGND,
  MATORDERING1WD,
  MATORDERINGRCM,
  MATORDERINGQMD,
  MATORDERINGROWLENGTH
#if PETSC_VERSION_GT(3,5,0)
  ,
  MATORDERINGWBM,
  MATORDERINGSPECTRAL,
  MATORDERINGAMD
#endif
};

template <typename T, typename I>
int
PetscLinearMatrixSolverImplementation<T, I>::p_nSupportedOrderingTypes =
  sizeof(p_supportedOrderingType)/sizeof(MatOrderingType);

template <typename T, typename I>
typename PetscLinearMatrixSolverImplementation<T, I>::ThePetscMatSolverType
PetscLinearMatrixSolverImplementation<T, I>::p_supportedSolverPackage[] = {
  MATSOLVERSUPERLU_DIST,
  MATSOLVERSUPERLU,
  MATSOLVERMUMPS,
  MATSOLVERPETSC
};

template <typename T, typename I>
int
PetscLinearMatrixSolverImplementation<T, I>::p_nSupportedSolverPackages = 
  sizeof(p_supportedSolverPackage)/sizeof(PetscLinearMatrixSolverImplementation<T, I>::ThePetscMatSolverType);

} // namespace math
} // namespace gridpack


#endif
