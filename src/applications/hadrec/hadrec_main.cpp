
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hadrec_main.cpp
 * @author Bruce Palmer
 * @date   2020-04-23 13:26:55 d3g096
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
#include <string>


const char* help = "HADREC Test application";

int main(int argc, char **argv)
{
  gridpack::parallel::Communicator world;
  int me = world.rank();
  // Initialize libraries (parallel and math)
  gridpack::NoPrint *noprint_ins = gridpack::NoPrint::instance();
  //  noprint_ins->setStatus(true);

  //bool bnoprint = gridpack::NoPrint::instance()->status();
  //printf ("------------- hadrec_main function test 1  bnoprint: %d \n", bnoprint);

  gridpack::Environment env(argc,argv,help);

  //bnoprint = gridpack::NoPrint::instance()->status();
  //printf ("------------- hadrec_main function test 2  bnoprint: %d \n", bnoprint);

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

  //bnoprint = gridpack::NoPrint::instance()->status();
  //printf ("------------- hadrec_main function test 3  bnoprint: %d \n", bnoprint);

  //hadrec_app_sptr->solvePowerFlowBeforeDynSimu(const_cast<char *>(file.c_str()), 2);

  hadrec_app_sptr->solvePowerFlowBeforeDynSimu(const_cast<char *>(file.c_str()));

  // transfer power flow results to dynamic simulation
  hadrec_app_sptr->transferPFtoDS();

  std::vector<gridpack::dynamic_simulation::Event> BusFaults;
  BusFaults.clear();

  // initialize dynamic simulation
  hadrec_app_sptr->initializeDynSimu(BusFaults);

  bool debugoutput = false; // whether print out debug staffs
  double lp, lq, pg, qg;
  int busno = 5;
  bool btmp;

  //-----test get load and get generator function-----------
  //if (debugoutput){

  busno = 5;
  btmp = hadrec_app_sptr->getBusTotalLoadPower(busno, lp, lq);
  if (me==0) printf("------------test hadrec_main, load at bus %d, has P: %f, Q: %f\n", busno, lp, lq);
  busno = 7;
  btmp = hadrec_app_sptr->getBusTotalLoadPower(busno, lp, lq);
  if (me==0) printf("------------test hadrec_main, load at bus %d, has P: %f, Q: %f\n", busno, lp, lq);
  busno = 9;
  btmp = hadrec_app_sptr->getBusTotalLoadPower(busno, lp, lq);
  if (me==0) printf("------------test hadrec_main, load at bus %d, has P: %f, Q: %f\n", busno, lp, lq);

  busno = 1;
  std::string genid = "1 ";
  btmp = hadrec_app_sptr->getGeneratorPower(busno, genid, pg, qg);
  if (me==0) printf("------------test hadrec_main, %d find generator at bus %d, has P: %f, Q: %f\n", btmp, busno, pg, qg);

  busno = 2;
  genid = "1 ";
  btmp = hadrec_app_sptr->getGeneratorPower(busno, genid, pg, qg);
  if (me==0) printf("------------test hadrec_main, %d find generator at bus %d, has P: %f, Q: %f\n", btmp, busno, pg, qg);

  busno = 3;
  genid = "1 ";
  btmp = hadrec_app_sptr->getGeneratorPower(busno, genid, pg, qg);
  if (me==0) printf("------------test hadrec_main, %d find generator at bus %d, has P: %f, Q: %f\n", btmp, busno, pg, qg);
  //}


  gridpack::hadrec::HADRECAction loadshedact;
  loadshedact.actiontype = 0;
  loadshedact.bus_number = 5;
  loadshedact.componentID = "1";
  loadshedact.percentage = -0.2;

  gridpack::hadrec::HADRECAction loadshedact1;
  loadshedact1.actiontype = 0;
  loadshedact1.bus_number = 7;
  loadshedact1.componentID = "1";
  loadshedact1.percentage = -0.2;

  gridpack::hadrec::HADRECAction linetrip;
  linetrip.actiontype = 1;
  //linetrip.bus_number = 6;
  linetrip.brch_from_bus_number = 6;
  linetrip.brch_to_bus_number = 7;
  linetrip.branch_ckt = "1 ";


  int isteps = 0;
  bool bApplyAct_LoadShedding = false;  // whether apply the load shedding action in the simulation steps
  bool bApplyAct_LineTripping = true;  // whether apply the line tripping action in the simulation steps
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

  if (debugoutput && me == 0){
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

    printf("-----------ob load bus list, ");
    for (idxtmp=0; idxtmp<obs_loadBus.size(); idxtmp++){
      printf(" %d, ", obs_loadBus[idxtmp]);
    }
    printf(" \n");

    printf("-----------ob load ID list, ");
    for (idxtmp=0; idxtmp<obs_loadIDs.size(); idxtmp++){
      printf(" %s, ", obs_loadIDs[idxtmp].c_str());
    }
    printf(" \n");

    printf("-----------ob bus list, ");
    for (idxtmp=0; idxtmp<obs_vBus.size(); idxtmp++){
      printf(" %d, ", obs_vBus[idxtmp]);
    }
    printf(" \n");
  }

  std::vector<int> zone_id;
  std::vector<double> tmp_p;
  std::vector<double> tmp_q;

  hadrec_app_sptr->getZoneLoads(tmp_p, tmp_q, zone_id);
  if (me == 0) {
    printf("\n-------------------get zone load information, total zones: %d \n\n", zone_id.size());
    for (idxtmp=0; idxtmp<zone_id.size(); idxtmp++){

      printf(" zone number: %d, total load p: %f, total load q: %f,\n", zone_id[idxtmp], tmp_p[idxtmp], tmp_q[idxtmp]);

    }
  }

  hadrec_app_sptr->getZoneGeneratorPower(tmp_p, tmp_q, zone_id);
  if (me == 0) {
    printf("\n-------------------get zone generation information, total zones: %d \n\n", zone_id.size());
    for (idxtmp=0; idxtmp<zone_id.size(); idxtmp++){

      printf(" zone number: %d, total generation p: %f, total generation q: %f,\n", zone_id[idxtmp], tmp_p[idxtmp], tmp_q[idxtmp]);

    }
  }

  while(!hadrec_app_sptr->isDynSimuDone()){
    // if the dynamic simulation is not done (hit the end time)
    if ( bApplyAct_LoadShedding && (isteps == 2500 || isteps == 3000 ||
          isteps == 3500 || isteps == 4000 ||
          isteps == 4500 || isteps == 5000 || isteps == 5500 ) ){
      //apply action
      hadrec_app_sptr->applyAction(loadshedact);
      hadrec_app_sptr->applyAction(loadshedact1);
      //printf("----renke debug load shed, isteps: %d \n", isteps);
    }

    if ( bApplyAct_LineTripping && isteps == 400){
      if (me == 0) printf("----renke debug line trip, isteps: %d \n", isteps);
      hadrec_app_sptr->applyAction(linetrip);
    }

    //execute one dynamic simulation step
    hadrec_app_sptr->executeDynSimuOneStep();

    ob_vals.clear();
    ob_vals = hadrec_app_sptr->getObservations();

    if (debugoutput && me == 0) {
      printf("observations, ");
      for (idxtmp=0; idxtmp<ob_vals.size(); idxtmp++) {
        printf(" %16.12f, ", ob_vals[idxtmp]);
      }
      printf(" \n");
    }

    isteps++;
  }

  if (me == 0) printf("\n----------------finished first round of dynamic simulation----\n ");
  //timer->stop(t_total);
  //timer->dump();


  //start the reload and second time dynamic simulation here
  // transfer power flow results to dynamic simulation
  bool btest_2dynasimu = false;
  if (btest_2dynasimu) {

    gridpack::dynamic_simulation::Event busfault;
    busfault.start = 1.0;
    busfault.end = 1.2;
    busfault.step = 0.005;
    busfault.isBus = true;
    busfault.bus_idx = 7;

    BusFaults.clear();
    BusFaults.push_back(busfault);

    //hadrec_app_sptr->solvePowerFlowBeforeDynSimu(argc, argv);
    //hadrec_app_sptr->solvePowerFlowBeforeDynSimu(const_cast<char *>(file.c_str()), 1);
    hadrec_app_sptr->solvePowerFlowBeforeDynSimu(const_cast<char *>(file.c_str()));

    if (me == 0) printf("\n---------------renke debug, hadrec main, second dyn starts----------------\n");
    hadrec_app_sptr->transferPFtoDS();

    // initialize dynamic simulation
    hadrec_app_sptr->initializeDynSimu(BusFaults);

    busno = 5;
    hadrec_app_sptr->getBusTotalLoadPower(busno, lp, lq);
    if (me == 0) printf("------------test hadrec_main, load at bus %d, has P: %f, Q: %f\n", busno, lp, lq);
    busno = 7;
    hadrec_app_sptr->getBusTotalLoadPower(busno, lp, lq);
    if (me == 0) printf("------------test hadrec_main, load at bus %d, has P: %f, Q: %f\n", busno, lp, lq);
    busno = 9;
    hadrec_app_sptr->getBusTotalLoadPower(busno, lp, lq);
    if (me == 0) printf("------------test hadrec_main, load at bus %d, has P: %f, Q: %f\n", busno, lp, lq);

    busno = 1;
    std::string genid = "1";
    btmp = hadrec_app_sptr->getGeneratorPower(busno, genid, pg, qg);
    if (me == 0) printf("------------test hadrec_main, %d find generator at bus %d, has P: %f, Q: %f\n", btmp, busno, pg, qg);

    busno = 2;
    genid = "1";
    btmp = hadrec_app_sptr->getGeneratorPower(busno, genid, pg, qg);
    if (me == 0) printf("------------test hadrec_main, %d find generator at bus %d, has P: %f, Q: %f\n", btmp, busno, pg, qg);

    busno = 3;
    genid = "1";
    btmp = hadrec_app_sptr->getGeneratorPower(busno, genid, pg, qg);
    if (me == 0) printf("------------test hadrec_main, %d find generator at bus %d, has P: %f, Q: %f\n", btmp, busno, pg, qg);

    isteps = 0;
    //bApplyAct_LoadShedding = true;  // whether apply the action in the simulation steps

    while(!hadrec_app_sptr->isDynSimuDone()){
      // if the dynamic simulation is not done (hit the end time)
      if ( bApplyAct_LoadShedding && (isteps == 2500 || isteps == 3000 ||
            isteps == 3500 || isteps == 4000 ) ){
        //apply action
        hadrec_app_sptr->applyAction(loadshedact);
        hadrec_app_sptr->applyAction(loadshedact1);
        //printf("----renke debug load shed, isteps: %d \n", isteps);
      }
      //execute one dynamic simulation step
      hadrec_app_sptr->executeDynSimuOneStep();

      ob_vals.clear();
      ob_vals = hadrec_app_sptr->getObservations();

      if (debugoutput && me == 0){
        printf("observations, ");
        for (idxtmp=0; idxtmp<ob_vals.size(); idxtmp++) {
          printf(" %16.12f, ", ob_vals[idxtmp]);
        }
        printf(" \n");
      }

      isteps++;
    }
  }

  //printf(" ----------------finished \n ");
  // Make sure we could do it again with a new instance if we wanted to
  hadrec_app_sptr.reset(new gridpack::hadrec::HADRECAppModule());
  hadrec_app_sptr.reset();

  timer->stop(t_total);
  //timer->dump();

  return 0;
}

