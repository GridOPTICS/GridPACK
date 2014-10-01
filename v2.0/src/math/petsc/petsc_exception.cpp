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
 * @date   2014-09-12 13:40:09 d3g096
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
#if PETSC_VERSION_LT(3,5,0)
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
#if PETSC_VERSION_LT(3,5,0)
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
