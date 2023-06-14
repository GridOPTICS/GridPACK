/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wind_main.cpp
 *
 * @brief Main driver for Wind DSA application
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "wind_driver.hpp"

// Calling program for dynamic security analysis

const char *help = "Dynamic Security Assessment Application for Wind Energy";

int
main(int argc, char **argv)
{
  // Initialize libraries (parallel and math)
  gridpack::Environment env(argc,argv,help);

  {
    gridpack::contingency_analysis::WindDriver driver;
    driver.execute2(argc, argv);
  }

  return 0;
}

