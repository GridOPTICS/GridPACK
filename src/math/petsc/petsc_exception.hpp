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
 * @date   2019-08-13 09:12:24 d3g096
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
#include <petscversion.h>

// You gotta love PETSc consistency 

#if PETSC_VERSION_(3,4,0)
#undef PETSC_VERSION_RELEASE
#define PETSC_VERSION_RELEASE 0
#endif

// Default type for PETSc exceptions
#define PETSC_EXCEPTION_TYPE std::runtime_error

// If PETSc is compiled with the C++, use its C++ exception
// facility. Otherwise, we need to make our own.

#ifdef PETSC_CLANGUAGE_CXX
// With PETSc version 3.5, PETSc::Exception was no longer defined.  It
// was replaced with std::runtime_error. 

#if PETSC_VERSION_LT(3,5,0)
#include <petscsys.hh>
#undef PETSC_EXCEPTION_TYPE
#define PETSC_EXCEPTION_TYPE PETSc::Exception
#define GRIDPACK_USES_PETSC_EXCEPTION 1
#endif
#endif

#include "gridpack/utilities/exception.hpp"

namespace gridpack {
namespace math {

#ifndef CHKERRXX

// If PETSc is not compiled as C++, the CHKERRXX macro is not defined.
// This does two things: (1) starts PETSc error handling, and (2)
// throws an exception.  GridPACK's PETSc interface depends heavily on
// it. So, make a replacement, that hopefully does the same thing.

extern void throw_petsc_exception(int line, const char *file, int ierr);

#define CHKERRXX(ierr)  do { if (PetscUnlikely(ierr)) gridpack::math::throw_petsc_exception(__LINE__, __FILE__, ierr); } while (0);

#endif



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
