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
  gridpack::parallel::Communicator world;
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Dynamic Simulation: Total Application");
  timer->start(t_total);

  std::string inputfile;
  if (argc >= 2 && argv[1] != NULL) {
    inputfile = argv[1];
  } else {
    inputfile = "input.xml";
  }
  gridpack::dynamic_simulation::DSFullApp ds_app;
  ds_app.solvePowerFlowBeforeDynSimu(inputfile.c_str());
  ds_app.readGenerators();
  ds_app.readSequenceData();
  //printf("ds_app.initialize:\n");
  ds_app.initialize();
  ds_app.setGeneratorWatch();
	
  ds_app.setup();
  int ngen = ds_app.numGenerators();
  int nload = ds_app.numLoads();
  int nline = ds_app.numLines();
  if (world.rank() == 0) {
    std::cout << " Number of generators: "<<ngen<<std::endl;
    std::cout << " Number of loads:      "<<nload<<std::endl;
    std::cout << " Number of lines:      "<<nline<<std::endl;
  }

  ds_app.run();

  //ds_app.write();
  timer->stop(t_total);
  timer->dump();
}

// Calling program for the dynamis simulation applications

int main(int argc, char **argv)
{
  gridpack::Environment env(argc, argv, NULL, 200000, 200000);

  gridpack::NoPrint *noprint_ins = gridpack::NoPrint::instance();
  noprint_ins->setStatus(false);
  
  run_dynamics(argc,argv);
  
}

