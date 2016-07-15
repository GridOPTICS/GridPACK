/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_test.cpp
 * @author Bruce Palmer
 * @date   2016-07-14 14:23:29 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/math/math.hpp"
#include "ds_app_module.hpp"

// Calling program for the contingency_analysis applications

int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc,argv);
  gridpack::math::Initialize();

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
    boost::shared_ptr<gridpack::dynamic_simulation::DSNetwork>
      ds_network(new gridpack::dynamic_simulation::DSNetwork(world));

    gridpack::dynamic_simulation::DSAppModule ds_app;
    ds_app.readNetwork(ds_network,config);
    ds_app.readGenerators();
    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config->getCursor("Configuration.Dynamic_simulation");
    std::vector<gridpack::dynamic_simulation::DSBranch::Event> faults;
    faults = ds_app.getFaults(cursor);
    ds_app.initialize();
    ds_app.setGeneratorWatch();
    ds_app.solve(faults[0]);
    ds_app.write();
  }

  // Terminate Math libraries
  gridpack::math::Finalize();
  return 0;
}

