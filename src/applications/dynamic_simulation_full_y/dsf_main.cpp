/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsf_main.cpp
 * @date   2023-11-08 09:22:59 d3g096
 *
 * @brief
 */
// -------------------------------------------------------------

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/environment/environment.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "gridpack/applications/modules/dynamic_simulation_full_y/dsf_app_module.hpp"
#include <vector>


// Calling program for the dynamis simulation applications

int
main(int argc, char **argv)
{
  gridpack::Environment env(argc, argv, NULL, 200000, 200000);

  gridpack::NoPrint *noprint_ins = gridpack::NoPrint::instance();
  noprint_ins->setStatus(false);

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

    gridpack::dynamic_simulation::DSFullApp ds_app;

    // Run powerflow and initialize DS network. Using solvePowerFlowBeforeDynSimu
    // encapsulates pf_network inside its own scope so it is destroyed cleanly
    // before DS simulation begins, avoiding GA handle lifecycle issues.
    if (argc >= 2 && argv[1] != NULL) {
      ds_app.solvePowerFlowBeforeDynSimu(argv[1]);
    } else {
      ds_app.solvePowerFlowBeforeDynSimu("input.xml");
    }

    ds_app.readGenerators();
    ds_app.readSequenceData();
    ds_app.initialize();
    ds_app.setGeneratorWatch();

    // read in faults from input file
    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config->getCursor("Configuration.Dynamic_simulation");
    std::vector<gridpack::dynamic_simulation::Event> faults;
    faults = ds_app.getEvents(cursor);

    ds_app.solvePreInitialize(faults[0]);

    while(!ds_app.isDynSimuDone()){
      ds_app.executeOneSimuStep( );
    }

    //ds_app.write();
    timer->stop(t_total);
    timer->dump();
  }

}



