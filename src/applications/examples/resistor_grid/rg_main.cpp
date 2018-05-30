/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   rg_main.cpp
 * @author Bruce Palmer
 * @date   2014-12-09 14:41:19 d3g096
 * 
 * @brief  
 */
// -------------------------------------------------------------

#include "gridpack/include/gridpack.hpp"
#include "rg_app.hpp"

// Calling program for the resistor_grid applications

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

  // Initialize Math libraries
  gridpack::math::Initialize();

  gridpack::resistor_grid::RGApp app;
  app.execute(argc, argv);

  // Terminate Math libraries
  gridpack::math::Finalize();
#if 0
  GA_Terminate();

  // Clean up MPI libraries
  ierr = MPI_Finalize();
#endif
  return 0;
}
