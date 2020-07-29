/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hadrec_app.cpp
 * @author Bruce Palmer
 * @date   2020-04-23 13:19:33 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "hadrec_app_module.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/math/math.hpp"
#include <vector>

/**
 * Basic constructor
 */
gridpack::hadrec::HADRECAppModule::HADRECAppModule(void)
  : config_sptr(new gridpack::utility::Configuration())
{
	bconfig_sptr_set = false;
}

/**
 * Basic destructor
 */
gridpack::hadrec::HADRECAppModule::~HADRECAppModule(void)
{
}

/**
 * solve power flow before run dynamic simulation 
 */
void gridpack::hadrec::HADRECAppModule::solvePowerFlowBeforeDynSimu(
    const char *inputfile, int pfcase_idx){
	
	gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
    //t_total = timer->createCategory("Dynamic Simulation: Total Application");
    //timer->start(t_total);

    gridpack::parallel::Communicator world;
	//world.reset(new gridpack::parallel::Communicator() );

    // read configuration file 
    t_config = timer->createCategory("Dynamic Simulation: Config");
    timer->start(t_config);
	
    //gridpack::utility::Configuration *config =
    //  gridpack::utility::Configuration::configuration();
	  
    //config_sptr.reset(gridpack::utility::Configuration::configuration());
	//config_sptr.reset( new gridpack::utility::Configuration() )
	
	if (!bconfig_sptr_set){
		if (inputfile) {
			config_sptr->open(inputfile, world);
		} else {
			config_sptr->open("input.xml",world);
		}
		bconfig_sptr_set = true;
	}
    timer->stop(t_config);

    // setup and run powerflow calculation
    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config_sptr->getCursor("Configuration.Powerflow");
    bool useNonLinear = false;
    useNonLinear = cursor->get("UseNonLinear", useNonLinear);
	
    pf_network.reset(new gridpack::powerflow::PFNetwork(world));

    pf_app_sptr.reset(new gridpack::powerflow::PFAppModule()) ;
	
	bool bnoprint = gridpack::NoPrint::instance()->status();

	if (bnoprint) { //if no print flag set to be true, suppress any intermedium output from pf app module
		pf_app_sptr->suppressOutput(true);
	}
	
    pf_app_sptr->readNetwork(pf_network, &(*config_sptr), pfcase_idx);
    pf_app_sptr->initialize();
    if (useNonLinear) {
      pf_app_sptr->nl_solve();
    } else {
      pf_app_sptr->solve();
    }
	
	if (!bnoprint) { //if no print flag set to be true, do not print the power flow solution
		pf_app_sptr->write();
	}
	
    pf_app_sptr->saveData();
	
    // setup dynamic simulation network
   
    ds_network.reset(new gridpack::dynamic_simulation::DSFullNetwork(world));
	//pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
    //  gridpack::dynamic_simulation::DSFullBranch>(ds_network);
	
}

/**
 * transfer data from power flow to dynamic simulation 
 */
void gridpack::hadrec::HADRECAppModule::transferPFtoDS(){
	  
    //ds_network.reset(new gridpack::dynamic_simulation::DSFullNetwork(world));
	pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
      gridpack::dynamic_simulation::DSFullBranch>(ds_network);
	
    int numBus = pf_network->numBuses();
    int i;
    gridpack::component::DataCollection *pfData;
    gridpack::component::DataCollection *dsData;
    double rval;
    for (i=0; i<numBus; i++) {
      pfData = pf_network->getBusData(i).get();
      dsData = ds_network->getBusData(i).get();
      pfData->getValue("BUS_PF_VMAG",&rval);
      dsData->setValue(BUS_VOLTAGE_MAG,rval);
      ///printf("Step0 bus%d mag = %f\n", i+1, rval);
      pfData->getValue("BUS_PF_VANG",&rval);
      dsData->setValue(BUS_VOLTAGE_ANG,rval);
      int ngen = 0;
      if (pfData->getValue(GENERATOR_NUMBER, &ngen)) {
        int j;
        for (j=0; j<ngen; j++) {
          pfData->getValue("GENERATOR_PF_PGEN",&rval,j);
          dsData->setValue(GENERATOR_PG,rval,j);
          //printf("save PGEN: %f\n", rval);
          pfData->getValue("GENERATOR_PF_QGEN",&rval,j);
          dsData->setValue(GENERATOR_QG,rval,j);
          //printf("save QGEN: %f\n", rval);
        }
      }
  }	
}

/**
 * do initialization only for dynamics simulation
 * @param faults lists of faults that might be simulated
 * @param dscase_idx index pointing to dyr parameter file that should be used if
 *                   a list of files is supplied in input deck
 */
void gridpack::hadrec::HADRECAppModule::initializeDynSimu
(std::vector<gridpack::dynamic_simulation::Event> faults,
 int dscase_idx){ 
// the definition of struct Event is at
// /src/applications/modules/dynamic_simulation_full_y/dsf_components.hpp

  ds_app_sptr.reset(new gridpack::dynamic_simulation::DSFullApp());

  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config_sptr->getCursor("Configuration.Dynamic_simulation");
  
  if (faults.empty()){
	faults = ds_app_sptr->getFaults(cursor);
  }

  // run dynamic simulation
  ds_app_sptr->setNetwork(ds_network, &(*config_sptr));
  //ds_app_sptr->readNetwork(ds_network,config);
  ds_app_sptr->readGenerators(dscase_idx);
  //printf("ds_app_sptr->initialize:\n");
  ds_app_sptr->initialize();
  ds_app_sptr->setGeneratorWatch();
  //printf("gen ID:	mac_ang_s0	mac_spd_s0	pmech	pelect\n");
  //printf("Step	time:	bus_id	mac_ang_s1	mac_spd_s1\n");
  //printf("ds_app_sptr->solve:\n");
  //ds_app_sptr->solve(faults[0]);
  ds_app_sptr->setObservations(cursor);
  p_obs_genBus.clear();
  p_obs_genIDs.clear();
  p_obs_loadBus.clear();
  p_obs_loadIDs.clear();
  p_obs_vBus.clear();
  //p_obs_vals.clear();
  ds_app_sptr->getObservationLists(p_obs_genBus, p_obs_genIDs,
      p_obs_loadBus, p_obs_loadIDs, p_obs_vBus);

  ds_app_sptr->solvePreInitialize(faults[0]);
}

/**
 * do a fully initialization before running dynamics simulation
 * @param case_idx index pointing to network configuration and dyr parameter
 * file that should be used if a list of files is supplied in input deck
 */
void gridpack::hadrec::HADRECAppModule::fullInitializationBeforeDynSimuSteps(
    const char *inputfile,
    const std::vector<gridpack::dynamic_simulation::Event>& BusFaults,
    int pfcase_idx, int dscase_idx){
	
	solvePowerFlowBeforeDynSimu(inputfile, pfcase_idx);
	transferPFtoDS();

	initializeDynSimu(BusFaults, dscase_idx);
		
}

/**
 * Execute only one simulation time step 
 */
void gridpack::hadrec::HADRECAppModule::executeDynSimuOneStep(){
		
	ds_app_sptr->executeOneSimuStep();
	
/*
   std::vector<double> rSpd, rAng, vMag, vAng;

   ds_app_sptr->getObservations(vMag, vAng, rSpd, rAng);
   int rank = ds_network->communicator().rank();
   if (rank == 0 && vMag.size() > 0) {
     int i;
     printf("Voltage observation, ");
     for (i=0; i<vMag.size(); i++) {
       printf(" %16.12f, %16.12f,",
           vMag[i],vAng[i]);
     }
     printf(" \n");
   }
   if (rank == 0 && rSpd.size() > 0) {
     int i;
     printf("Generator observation, ");
     for (i=0; i<rSpd.size(); i++) {
       printf(" %16.12f,  %16.12f, ",
           rSpd[i],rAng[i]);
     }
     printf(" \n");
   }
   
*/
}

/**
 * return observations after each simulation time step
 */

std::vector<double> gridpack::hadrec::HADRECAppModule::getObservations(){
	
	std::vector<double> rSpd, rAng, vMag, vAng, fOnline;
	ds_app_sptr->getObservations(vMag, vAng, rSpd, rAng, fOnline);
	
	std::vector<double> obs_vals;
	obs_vals.clear();
	
	int i;
	for (i=0; i<rSpd.size(); i++){
		obs_vals.push_back(rSpd[i]);
	}
	
	for (i=0; i<rAng.size(); i++){
		obs_vals.push_back(rAng[i]);
	}

	for (i=0; i<vMag.size(); i++){
		obs_vals.push_back(vMag[i]);
	}

	for (i=0; i<vAng.size(); i++){
		obs_vals.push_back(vAng[i]);
	}	

	for (i=0; i<fOnline.size(); i++){
		obs_vals.push_back(fOnline[i]);
	}	
	
	return obs_vals;
	
}

/**
 * return observations list
 */

void gridpack::hadrec::HADRECAppModule::getObservationLists(
    std::vector<int> &genBuses, std::vector<std::string> &genIDs,
    std::vector<int> &loadBuses, std::vector<std::string> &loadIDs,
    std::vector<int> &busIDs){
	
	ds_app_sptr->getObservationLists(genBuses, genIDs, loadBuses,
       loadIDs, busIDs);
	
}

/**
 * Check whether the dynamic simulation is done
 */
bool gridpack::hadrec::HADRECAppModule::isDynSimuDone( ){
	
	return ds_app_sptr->isDynSimuDone();

}

/**
 * apply actions
 */
void gridpack::hadrec::HADRECAppModule::applyAction(gridpack::hadrec::HADRECAction control_action){
	if ( control_action.actiontype == 0 ){
		return ds_app_sptr->applyLoadShedding(control_action.bus_number, control_action.componentID, control_action.percentage );
	}

}
