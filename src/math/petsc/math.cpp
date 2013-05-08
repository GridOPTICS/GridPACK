/**
 * @file   math.cpp
 * @author William A. Perkins
 * @date   2013-05-08 14:06:21 d3g096
 * 
 * @brief  
 * 
 * 
 */

#include <petscsys.h>
#include "gridpack/math/petsc/petsc_exception.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
// Initialize
// -------------------------------------------------------------
/// Does whatever is necessary to start up the PETSc library
void
Initialize(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = PetscInitializeNoArguments(); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

// -------------------------------------------------------------
// Initialized
// -------------------------------------------------------------
bool
Initialized(void)
{
  PetscErrorCode ierr(0);
  PetscBool result;
  try {
    ierr = PetscInitialized(&result); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  return result;
}

// -------------------------------------------------------------
// Finalize
// -------------------------------------------------------------
/// Does whatever is necessary to shut down the PETSc library
void
Finalize(void)
{
  PetscErrorCode ierr(0);
  try {
    ierr = PetscFinalize(); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

} // namespace math
} // namespace gridpack
