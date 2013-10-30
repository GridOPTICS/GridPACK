// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   math.cpp
 * @author William A. Perkins
 * @date   2013-10-09 13:22:31 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <petscsys.h>
#include "gridpack/math/math.hpp"
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
  if (Initialized()) return;
  PetscErrorCode ierr(0);
  try {
    ierr = PetscInitializeNoArguments(); CHKERRXX(ierr);
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD,
                                  "gridpack.petscrc",
                                  PETSC_FALSE); CHKERRXX(ierr);
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
  if (!Initialized()) return;
  PetscErrorCode ierr(0);
  try {
    ierr = PetscFinalize(); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
}

} // namespace math
} // namespace gridpack
