// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_exception.hpp
 * @author William A. Perkins
 * @date   2014-09-29 09:47:51 d3g096
 * 
 * @brief 
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_exception_hpp_
#define _petsc_exception_hpp_

#include <string>
#include <exception>
#include <petscsys.h>

// You gotta love PETSc consistency 

#if PETSC_VERSION_(3,4,0)
#undef PETSC_VERSION_RELEASE
#define PETSC_VERSION_RELEASE 0
#endif

// With PETSc version 3.5, PETSc::Exception was no longer defined.  It
// was replaced with std::runtime_error. 

#if PETSC_VERSION_LT(3,5,0)
#include <petscsys.hh>
#define PETSC_EXCEPTION_TYPE PETSc::Exception
#else
#define PETSC_EXCEPTION_TYPE std::runtime_error
#endif 

#include "gridpack/utilities/exception.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScException
// -------------------------------------------------------------
class PETScException : public gridpack::Exception {
public:

  /// Construct with a PETSc error code.
  explicit PETScException(const PetscErrorCode ierr);

  /// Construct from a PETSc exception
  explicit PETScException(const PETSC_EXCEPTION_TYPE& e);

  /// Construct from a PETSc exception and error code
  PETScException(const PetscErrorCode& ierr, const PETSC_EXCEPTION_TYPE& e);

  /// Copy constructor
  PETScException(const PETScException& old);

  /// Destructor
  ~PETScException(void) throw();

protected:
  
  /// PETSc error code causing the problem
  PetscErrorCode petsc_err_;

  
};





} // namespace math
} // namespace gridpack

#endif
