/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_main.cpp
 * @author Bruce Palmer
 * @date   2016-07-14 14:23:28 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/voltage_stability/vs_app_module.hpp"

// Calling program for the state estimation application

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

  if (1) {
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

    // setup and run state estimation calculation
    boost::shared_ptr<gridpack::voltage_stability::VSNetwork>
      se_network(new gridpack::voltage_stability::VSNetwork(world));

    gridpack::voltage_stability::VSAppModule vs_app;
    /*
    vs_app.readNetwork(se_network,config);
    vs_app.initialize();
    vs_app.readMeasurements();
    vs_app.solve();
    vs_app.write();
    */
    vs_app.solvepowerflow();
    printf("debug 1. Solve power flow and transfer complete. Writing FVSI \n");
    vs_app.write_FVSI();
  }

  GA_Terminate();

  // Terminate Math libraries
  gridpack::math::Finalize();
  // Clean up MPI libraries
  ierr = MPI_Finalize();
}

