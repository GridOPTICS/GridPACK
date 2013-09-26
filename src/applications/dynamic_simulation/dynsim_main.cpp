// -------------------------------------------------------------
/**
 * @file   dynsim_main.cpp
 * @author Shuangshuang Jin
 * @date   September 19, 2013
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/math/math.hpp"
#include "gridpack/applications/dynamic_simulation/dynsim_app.hpp"

// Calling program for the dynamis simulation applications

main(int *argc, char **argv)
{
  // Initialize MPI libraries
  int ierr = MPI_init(&argc, &argv);
  // Intialize Math libraries
  gridpack::math::initialize();

  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

  gridpack::dynamic_simulation::DynSimApp app;
  app.execute();

  GA_Terminate();

  // Terminate Math libraries
  gridpack::math::Finalize();
  // Clean up MPI libraries
  ierr = MPI_Finalize();
}

