/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   rtpr_main.cpp
 * @author Bruce Palemr
 * @date   2019-10-09 13:25:25 d3g293
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/math/math.hpp"
#include "rtpr_driver.hpp"

// Calling program for the contingency_analysis applications

int
main(int argc, char **argv)
{
  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);

  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

  // Intialize Math libraries
  gridpack::math::Initialize(&argc,&argv);

  {
    gridpack::rtpr::RTPRDriver driver;
    driver.execute(argc, argv);
  }

  // Terminate Math libraries
  gridpack::math::Finalize();

  GA_Terminate();
  // Clean up MPI libraries
  ierr = MPI_Finalize();
  //return ierr;
}

