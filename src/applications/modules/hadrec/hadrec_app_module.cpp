/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hadrec_app.cpp
 * @author Bruce Palmer
 * @date   2018-06-20 11:07:20 d3g096
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
{
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
void gridpack::hadrec::HADRECAppModule::solvePowerFlowBeforeDynSimu(int argc, char **argv){
	
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
	  
	config_sptr.reset(gridpack::utility::Configuration::configuration());
	
    if (argc >= 2 && argv[1] != NULL) { 
      char inputfile[256]; 
      sprintf(inputfile,"%s",argv[1]);
      config_sptr->open(inputfile, world);
    } else {
      config_sptr->open("input.xml",world);
    }
    timer->stop(t_config);

    // setup and run powerflow calculation
    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config_sptr->getCursor("Configuration.Powerflow");
    bool useNonLinear = false;
    useNonLinear = cursor->get("UseNonLinear", useNonLinear);

    pf_network.reset(new gridpack::powerflow::PFNetwork(world));

    pf_app_sptr.reset(new gridpack::powerflow::PFAppModule()) ;
    pf_app_sptr->readNetwork(pf_network, &(*config_sptr) );
    pf_app_sptr->initialize();
    if (useNonLinear) {
      pf_app_sptr->nl_solve();
    } else {
      pf_app_sptr->solve();
    }
    pf_app_sptr->write();
    pf_app_sptr->saveData();
	
    // setup dynamic simulation network
   
    ds_network.reset(new gridpack::dynamic_simulation::DSFullNetwork(world));
	pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
      gridpack::dynamic_simulation::DSFullBranch>(ds_network);
	
}

/**
 * transfer data from power flow to dynamic simulation 
 */
void gridpack::hadrec::HADRECAppModule::transferPFtoDS(){
	  
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
 */
void gridpack::hadrec::HADRECAppModule::initializeDynSimu(){
	
	ds_app_sptr.reset(new gridpack::dynamic_simulation::DSFullApp());
	
	gridpack::utility::Configuration::CursorPtr cursor;
	cursor = config_sptr->getCursor("Configuration.Dynamic_simulation");
    std::vector<gridpack::dynamic_simulation::Event> faults;
    faults = ds_app_sptr->getFaults(cursor);

    // run dynamic simulation
    ds_app_sptr->setNetwork(ds_network, &(*config_sptr));
    //ds_app_sptr->readNetwork(ds_network,config);
    ds_app_sptr->readGenerators();
    //printf("ds_app_sptr->initialize:\n");
    ds_app_sptr->initialize();
    ds_app_sptr->setGeneratorWatch();
    //printf("gen ID:	mac_ang_s0	mac_spd_s0	pmech	pelect\n");
    //printf("Step	time:	bus_id	mac_ang_s1	mac_spd_s1\n");
    //printf("ds_app_sptr->solve:\n");
    //ds_app_sptr->solve(faults[0]);
	
	ds_app_sptr->solvePreInitialize(faults[0]);
	
	
}

/**
 * do a fully initialization before running dynamics simulation
 */
void gridpack::hadrec::HADRECAppModule::fullInitializationBeforeDynSimuSteps(int argc, char **argv){
	
	solvePowerFlowBeforeDynSimu(argc, argv);
	transferPFtoDS();
	initializeDynSimu();
		
}

/**
 * Execute only one simulation time step 
 */
void gridpack::hadrec::HADRECAppModule::executeDynSimuOneStep(){
		
	ds_app_sptr->executeOneSimuStep();
	
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
