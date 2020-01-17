
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hadrec_main.cpp
 * @author Bruce Palmer
 * @date   2016-07-14 14:23:07 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"

const char* help = "HADREC Test application";

int main(int argc, char **argv)
{
  // Initialize libraries (parallel and math)
  gridpack::Environment env(argc,argv,help);

  if (1) {
    gridpack::utility::CoarseTimer *timer =
      gridpack::utility::CoarseTimer::instance();
    gridpack::parallel::Communicator world;

    // read configuration file
    gridpack::utility::Configuration *config =
      gridpack::utility::Configuration::configuration();
    if (argc >= 2 && argv[1] != NULL) {
      char inputfile[256];
      sprintf(inputfile,"%s",argv[1]);
      config->open(inputfile,world);
    } else {
      config->open("input.xml",world);
    }
  }

  return 0;
}

