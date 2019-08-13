// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_exception.cpp
 * @author William A. Perkins
 * @date   2019-08-13 08:56:30 d3g096
 * 
 * @brief  Implementation of PETScException
 * 
 * 
 */
// -------------------------------------------------------------

#include <sstream>
#include "petsc_exception.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
// throw_petsc_exception
//
// This needs to throw the exception #define'd by
// PETSC_EXCEPTION_TYPE, not PETScException to be compatible with the
// way CHRERRXX (when defined by PETSc) is used.
// -------------------------------------------------------------
void
throw_petsc_exception(int line, const char *file, int ierr)
{
  std::ostringstream msg;
  const char *buf;
  PetscErrorMessage(ierr, &buf, PETSC_NULL);
  msg << file << ": " << line << ": "
      << "PETSc Error (" << ierr << "): " << buf;
  throw PETSC_EXCEPTION_TYPE(msg.str());
}

// -------------------------------------------------------------
//  class PETScException
// -------------------------------------------------------------

// -------------------------------------------------------------
// PETScException:: constructors / destructor
// -------------------------------------------------------------
PETScException::PETScException(const PetscErrorCode ierr)
  : gridpack::Exception(), petsc_err_(ierr)
{
  std::ostringstream msg;
  msg << "PETSc error (" << petsc_err_ << ")";
  const char *buf;
  PetscErrorMessage(petsc_err_, &buf, PETSC_NULL);
  msg << ": " << buf;
  message_ = msg.str();
}          

PETScException::PETScException(const PETSC_EXCEPTION_TYPE& e)
  : gridpack::Exception(), petsc_err_(0)
{
  std::ostringstream msg;

  msg << "PETSc error: " 
#ifdef GRIDPACK_USES_PETSC_EXCEPTION
      << e.msg()
#else
      << e.what()
#endif
    ;
  message_ = msg.str();
}          

PETScException::PETScException(const PetscErrorCode& ierr, const PETSC_EXCEPTION_TYPE& e)
  : gridpack::Exception(), petsc_err_(0)
{
  std::ostringstream msg;
  msg << "PETSc error (" << petsc_err_ << "): " 
#if GRIDPACK_USES_PETSC_EXCEPTION
      << e.msg()
#else
      << e.what()
#endif
    ;
  message_ = msg.str();
}          

PETScException::PETScException(const PETScException& old)
  : gridpack::Exception(old), petsc_err_(old.petsc_err_)
{
}

PETScException::~PETScException(void) throw()
{
}


} // namespace math
} // namespace gridpack
