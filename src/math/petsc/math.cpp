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
 * @date   2014-01-15 10:06:10 d3g096
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
  PetscBool flg;
  try {
// Turn this on to enable PETSc logging.
#if 0
    int argc = 2;
    char **args;
    args = new char*[2];
    args[0] = new char[32];
    args[1] = new char[32];
    sprintf(args[0],"powerflow.x");
    sprintf(args[1],"-log_summary");
    ierr = PetscInitialize(&argc,&args,NULL,NULL); CHKERRXX(ierr);
    delete [] args[1];
    delete [] args[0];
    delete [] args;
#else
    ierr = PetscInitializeNoArguments(); CHKERRXX(ierr);
#endif
    PetscOptionsHasName(NULL,"-log_summary",&flg);
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD,
                                  "gridpack.petscrc",
                                  PETSC_FALSE); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  // Print out some information on processor configuration
  int mpi_err, me, nproc;
  mpi_err = MPI_Comm_rank(PETSC_COMM_WORLD, &me);
  mpi_err = MPI_Comm_size(PETSC_COMM_WORLD, &nproc);
  if (me == 0) {
    printf("\nGridPACK math module configured on %d processors\n",nproc);
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
