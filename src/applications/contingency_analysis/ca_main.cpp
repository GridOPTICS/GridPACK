/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ca_main.cpp
 * @author Yousu Chen 
 * @date   January 20, 2014
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/math/math.hpp"
#include "gridpack/applications/contingency_analysis/ca_app.hpp"

// Calling program for the contingency_analysis applications

main(int argc, char **argv)
{
  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);
  // Intialize Math libraries
  gridpack::math::Initialize();

  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

  gridpack::contingency_analysis::CAApp app;
  app.execute(argc, argv);

  GA_Terminate();

  // Terminate Math libraries
  gridpack::math::Finalize();
  // Clean up MPI libraries
  ierr = MPI_Finalize();
}

