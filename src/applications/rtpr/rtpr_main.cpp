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
#include "gridpack/include/gridpack.hpp"
#include "rtpr_driver.hpp"

#ifdef USE_GOSS
const char *help = "Real Time Path Rating Application\n"
                   "Runtime input options:\n"
                   "-uri: URI parameter for GOSS server\n"
                   "-username: user name for GOSS server\n"
                   "-passwd: password for GOSS server\n"
                   "-input: topic name for initial input\n";
#else
const char *help = "Real Time Path Rating Application\n";              
#endif

// Calling program for the contingency_analysis applications

int
main(int argc, char **argv)
{
  gridpack::Environment env(argc, argv, help);

  {
    gridpack::rtpr::RTPRDriver driver;
    driver.execute(argc, argv, env);
  }
}

