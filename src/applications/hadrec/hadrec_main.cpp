
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
#include "gridpack/applications/modules/hadrec/hadrec_app_module.hpp"
#include "gridpack/applications/modules/dynamic_simulation_full_y/dsf_app_module.hpp"
#include <vector>


const char* help = "HADREC Test application";

int main(int argc, char **argv)
{
  // Initialize libraries (parallel and math)
  gridpack::Environment env(argc,argv,help);
  
  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);

  GA_Initialize();
  int stack = 200000, heap = 200000;
  MA_init(C_DBL, stack, heap);

  // Intialize Math libraries
  gridpack::math::Initialize(&argc,&argv);

  if (1) {
	  
/*
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
*/
	gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
    int t_total = timer->createCategory("Dynamic Simulation: Total Application");
	
	boost::shared_ptr<gridpack::hadrec::HADRECAppModule> hadrec_app_sptr (new gridpack::hadrec::HADRECAppModule() );
	hadrec_app_sptr->solvePowerFlowBeforeDynSimu(argc, argv);
	hadrec_app_sptr->transferPFtoDS();
	hadrec_app_sptr->initializeDynSimu();
	
	std::vector<gridpack::dynamic_simulation::Event> action_list;
	action_list.clear();
	
	while(!hadrec_app_sptr->isDynSimuDone()){
		hadrec_app_sptr->executeDynSimuOneStep(action_list);
	}
	
	timer->stop(t_total);
    timer->dump();

  }
  
  GA_Terminate();

  // Terminate Math libraries
  gridpack::math::Finalize();
  // Clean up MPI libraries
  ierr = MPI_Finalize();

  return 0;
}

