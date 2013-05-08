/**
 * @file   petsc_exception.cpp
 * @author William A. Perkins
 * @date   2013-05-08 13:45:22 d3g096
 * 
 * @brief  Implementation of PETScException
 * 
 * 
 */

#include <sstream>
#include "gridpack/math/petsc/petsc_exception.hpp"

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

PETScException::PETScException(const PETSc::Exception& e)
  : gridpack::Exception(), petsc_err_(0)
{
  std::ostringstream msg;
  msg << "PETSc error: " << e.msg();
  message_ = msg.str();
}          

PETScException::PETScException(const PetscErrorCode& ierr, const PETSc::Exception& e)
  : gridpack::Exception(), petsc_err_(0)
{
  std::ostringstream msg;
  msg << "PETSc error (" << petsc_err_ << "): " << e.msg();
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
