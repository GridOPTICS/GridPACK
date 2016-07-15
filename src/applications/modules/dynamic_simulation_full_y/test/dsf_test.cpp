/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_main.cpp
 * @author Shuangshuang Jin
 * @date   2016-07-14 14:30:38 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "gridpack/applications/modules/dynamic_simulation_full_y/dsf_app_module.hpp"

// test program for dynamic simulation module

int
main(int argc, char **argv)
{
  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);

  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

  // Intialize Math libraries
  gridpack::math::Initialize();

  if (1) {
    gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
    int t_total = timer->createCategory("Dynamic Simulation: Total Application");
    timer->start(t_total);

    gridpack::parallel::Communicator world;

    // read configuration file 
    int t_config = timer->createCategory("Dynamic Simulation: Config");
    timer->start(t_config);
    gridpack::utility::Configuration *config =
      gridpack::utility::Configuration::configuration();
    if (argc >= 2 && argv[1] != NULL) { 
      char inputfile[256]; 
      sprintf(inputfile,"%s",argv[1]);
      config->open(inputfile,world);
    } else {
      config->open("input.xml",world);
    }
    timer->stop(t_config);

    // setup and run dynamic simulation calculation
    gridpack::utility::Configuration::CursorPtr cursor;
   
    // setup and run dynamic simulation calculation
    boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork>
      ds_network(new gridpack::dynamic_simulation::DSFullNetwork(world));
    gridpack::dynamic_simulation::DSFullApp ds_app;

    // read in network
    cursor = config->getCursor("Configuration.Dynamic_simulation");
    std::string basefile;
    cursor->get("networkConfiguration",&basefile);
    ds_app.readNetwork(ds_network,config,basefile.c_str());

    // read in information on generators
    ds_app.readGenerators();

    // read in faults from input file
    cursor = config->getCursor("Configuration.Dynamic_simulation");
    std::vector<gridpack::dynamic_simulation::DSFullBranch::Event> faults;
    faults = ds_app.getFaults(cursor);

    // run dynamic simulation
    ds_app.initialize();
    ds_app.solve(faults[0]);
    // ds_app.write();
    timer->stop(t_total);
    timer->dump();
  }

  GA_Terminate();

  // Terminate Math libraries
  gridpack::math::Finalize();
  // Clean up MPI libraries
  ierr = MPI_Finalize();
}

