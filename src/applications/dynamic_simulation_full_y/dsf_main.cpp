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

    // setup and run powerflow calculation
    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config->getCursor("Configuration.Powerflow");
    bool useNonLinear = false;
    useNonLinear = cursor->get("UseNonLinear", useNonLinear);

    boost::shared_ptr<gridpack::powerflow::PFNetwork>
      pf_network(new gridpack::powerflow::PFNetwork(world));

    gridpack::powerflow::PFAppModule pf_app;
    pf_app.readNetwork(pf_network, config);
    pf_app.initialize();
    if (useNonLinear) {
      pf_app.nl_solve();
    } else {
      pf_app.solve();
    }
    pf_app.write();
    pf_app.saveData();
   
    // setup and run dynamic simulation calculation
    boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork>
      ds_network(new gridpack::dynamic_simulation::DSFullNetwork(world));
    gridpack::dynamic_simulation::DSFullApp ds_app;
    pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
      gridpack::dynamic_simulation::DSFullBranch>(ds_network);

    // transfer results from PF calculation to DS calculation
    ds_app.transferPFtoDS(pf_network, ds_network); 

    // run dynamic simulation
    ds_app.setNetwork(ds_network, config);
    //ds_app.readNetwork(ds_network,config);
    ds_app.readGenerators();
    ds_app.readSequenceData();
    //printf("ds_app.initialize:\n");
    ds_app.initialize();
    ds_app.setGeneratorWatch();
    //printf("gen ID:	mac_ang_s0	mac_spd_s0	pmech	pelect\n");
    //printf("Step	time:	bus_id	mac_ang_s1	mac_spd_s1\n");
    //printf("ds_app.solve:\n");
    //ds_app.solve(faults[0]);

    // read in faults from input file
    //gridpack::utility::Configuration::CursorPtr cursor;
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

