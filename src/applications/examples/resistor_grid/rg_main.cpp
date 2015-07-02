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

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/math/math.hpp"
#include "rg_app.hpp"

// Calling program for the resistor_grid applications

int
main(int argc, char **argv)
{
  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);

  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

  // Initialize Math libraries
  gridpack::math::Initialize();

  gridpack::resistor_grid::RGApp app;
  app.execute(argc, argv);

  GA_Terminate();

  // Terminate Math libraries
  gridpack::math::Finalize();
  // Clean up MPI libraries
  ierr = MPI_Finalize();
  return ierr;
}
