
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

  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Dynamic Simulation: Total Application");

  //allocate the memory for the HADRECAppModule
  boost::shared_ptr<gridpack::hadrec::HADRECAppModule>
    hadrec_app_sptr (new gridpack::hadrec::HADRECAppModule() );

  // solve power flow
  std::string file;
  if (argc > 1) {
    file = argv[1];
  } else {
    file = "input.xml";
  }

  hadrec_app_sptr->solvePowerFlowBeforeDynSimu(const_cast<char *>(file.c_str()));

  // transfer power flow results to dynamic simulation
  hadrec_app_sptr->transferPFtoDS();

  // initialize dynamic simulation
  hadrec_app_sptr->initializeDynSimu();

  gridpack::hadrec::HADRECAction loadshedact;
  loadshedact.actiontype = 0;
  loadshedact.bus_number = 5;
  loadshedact.componentID = "1";
  loadshedact.percentage = -0.2;

  int isteps = 0;
  bool bApplyAct = true; //false;  // whether apply the action in the simulation steps
  std::vector<double> ob_vals;
  int idxtmp;
  
  // test getOblist
  std::vector<int> obs_genBus;
  std::vector<std::string> obs_genIDs;
  std::vector<int> obs_loadBus;
  std::vector<std::string> obs_loadIDs;
  std::vector<int> obs_vBus;
  hadrec_app_sptr->getObservationLists(obs_genBus, obs_genIDs,
      obs_loadBus, obs_loadIDs, obs_vBus);
  
  printf("-----------renke debug, getObservationLists------------\n");
  printf("-----------ob gen bus list, ");
  for (idxtmp=0; idxtmp<obs_genBus.size(); idxtmp++){
	  printf(" %d, ", obs_genBus[idxtmp]);
  }
  printf(" \n");
  
  printf("-----------ob gen ID list, ");
  for (idxtmp=0; idxtmp<obs_genIDs.size(); idxtmp++){
	  printf(" %s, ", obs_genIDs[idxtmp].c_str());
  }
  printf(" \n");
  
  printf("-----------ob bus list, ");
  for (idxtmp=0; idxtmp<obs_vBus.size(); idxtmp++){
	  printf(" %d, ", obs_vBus[idxtmp]);
  }
  printf(" \n");

  while(!hadrec_app_sptr->isDynSimuDone()){
    // if the dynamic simulation is not done (hit the end time)
    if ( bApplyAct && (isteps == 2500 || isteps == 3000 ||
          isteps == 3500 || isteps == 4000 ) ){
      //apply action
      hadrec_app_sptr->applyAction(loadshedact);
      //printf("----renke debug load shed, isteps: %d \n", isteps);
    }
    //execute one dynamic simulation step
    hadrec_app_sptr->executeDynSimuOneStep();
	
	ob_vals.clear();
	ob_vals = hadrec_app_sptr->getObservations();
	
    printf("observations, ");
    for (idxtmp=0; idxtmp<ob_vals.size(); idxtmp++) {
       printf(" %16.12f, ", ob_vals[idxtmp]);
    }
    printf(" \n");

    isteps++;
  }

  //timer->stop(t_total);
  //timer->dump();


  //start the reload and second time dynamic simulation here
  // transfer power flow results to dynamic simulation
  bool btest_2dynasimu = true; //true;
  if (btest_2dynasimu) {
    
    //hadrec_app_sptr->solvePowerFlowBeforeDynSimu(argc, argv);
    
    printf("---------------renke debug, hadrec main, second dyn starts----------------\n");
    hadrec_app_sptr->transferPFtoDS();

    // initialize dynamic simulation
    hadrec_app_sptr->initializeDynSimu();

    isteps = 0;
    //bApplyAct = true;  // whether apply the action in the simulation steps

    while(!hadrec_app_sptr->isDynSimuDone()){
      // if the dynamic simulation is not done (hit the end time)
      if ( bApplyAct && (isteps == 2500 || isteps == 3000 ||
            isteps == 3500 || isteps == 4000 ) ){
        //apply action
        hadrec_app_sptr->applyAction(loadshedact);
        //printf("----renke debug load shed, isteps: %d \n", isteps);
      }
      //execute one dynamic simulation step
      hadrec_app_sptr->executeDynSimuOneStep();
	  
	  ob_vals.clear();
	  ob_vals = hadrec_app_sptr->getObservations();
	
      printf("observations, ");
      for (idxtmp=0; idxtmp<ob_vals.size(); idxtmp++) {
         printf(" %16.12f, ", ob_vals[idxtmp]);
      }
      printf(" \n");
	
      isteps++;
    }
  }

  timer->stop(t_total);
  timer->dump();

  return 0;
}

