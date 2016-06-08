/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_main.cpp
 * @author Bruce Palmer
 * @date   2014-12-09 14:39:50 d3g096
 * 
 * @brief  
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "pf_app.hpp"

// Calling program for the powerflow applications

int
main(int argc, char **argv)
{

  gridpack::parallel::Environment env(argc,argv);
#if 0
  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);
  GA_Initialize();

  int stack = 8000000, heap = 8000000;
  MA_init(C_DBL, stack, heap);
#endif

  // Initialize Math libraries
  gridpack::math::Initialize();

  gridpack::powerflow::PFApp app;
  app.execute(argc, argv);

  // Terminate Math libraries
  gridpack::math::Finalize();

  GA_Terminate();
#if 0
  // Clean up MPI libraries
  ierr = MPI_Finalize();
  return ierr;
#else
  return 0;
#endif
}
