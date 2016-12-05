/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hw_main.cpp
 * @author Bruce Palmer
 * @date   2014-12-09 14:40:37 d3g096
 * 
 * @brief  
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "hw_app.hpp"

// Calling program for the hello_world applications

int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
#if 0
  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);

  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);
#endif

  gridpack::hello_world::HWApp app;
  app.execute(argc, argv);

#if 0
  GA_Terminate();

  // Terminate Math libraries
  gridpack::math::Finalize();
  // Clean up MPI libraries
  ierr = MPI_Finalize();
  return ierr;
#endif
  return 0;
}
