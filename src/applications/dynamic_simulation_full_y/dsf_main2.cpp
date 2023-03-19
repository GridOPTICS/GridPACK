/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "gridpack/applications/modules/dynamic_simulation_full_y/dsf_app_module.hpp"
#include <vector>

void run_dynamics(int argc, char **argv)
{
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
	
  ds_app.setup();

  ds_app.run(0.5);

  //ds_app.write();
  timer->stop(t_total);
  timer->dump();
}

// Calling program for the dynamis simulation applications

int main(int argc, char **argv)
{
  gridpack::NoPrint *noprint_ins = gridpack::NoPrint::instance();
  noprint_ins->setStatus(false);
  
  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);

  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

  // Intialize Math libraries
  gridpack::math::Initialize(&argc,&argv);

  run_dynamics(argc,argv);
  
  GA_Terminate();

  // Terminate Math libraries
  gridpack::math::Finalize();
  // Clean up MPI libraries
  ierr = MPI_Finalize();
}

