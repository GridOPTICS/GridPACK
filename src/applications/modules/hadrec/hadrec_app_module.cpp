/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hadrec_app.cpp
 * @author Bruce Palmer
 * @date   2024-04-18 14:17:01 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "hadrec_app_module.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/math/math.hpp"
#include <gridpack/utilities/exception.hpp>
#include <vector>

/**
 * Basic constructor
 */
gridpack::hadrec::HADRECAppModule::HADRECAppModule(void)
  : config_sptr(new gridpack::utility::Configuration())
{
	bconfig_sptr_set = false;
	p_PFuseNonLinear = false;

	ds_app_sptr.reset(new gridpack::dynamic_simulation::DSFullApp());
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

    // read configuration file 
    t_config = timer->createCategory("Dynamic Simulation: Config");
    timer->start(t_config);
	
    //gridpack::utility::Configuration *config =
    //  gridpack::utility::Configuration::configuration();
	  
    //config_sptr.reset(gridpack::utility::Configuration::configuration());
	//config_sptr.reset( new gridpack::utility::Configuration() )
	
	if (!bconfig_sptr_set){
		if (inputfile) {
                  config_sptr->open(inputfile, this->communicator());
		} else {
                  config_sptr->open("input.xml", this->communicator());
		}
		bconfig_sptr_set = true;
	}
    timer->stop(t_config);

    // setup and run powerflow calculation
    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config_sptr->getCursor("Configuration.Powerflow");
    bool useNonLinear = false;
    useNonLinear = cursor->get("UseNonLinear", useNonLinear);
	
    pf_network.reset(new gridpack::powerflow::PFNetwork(this->communicator()));

    pf_app_sptr.reset(new gridpack::powerflow::PFAppModule()) ;
	
	bool bnoprint = gridpack::NoPrint::instance()->status();

	if (bnoprint) { //if no print flag set to be true, suppress any intermedium output from pf app module
		pf_app_sptr->suppressOutput(true);
	}
	
    pf_app_sptr->readNetwork(pf_network, &(*config_sptr), pfcase_idx);
    pf_app_sptr->initialize();
    pf_analytics.reset(new gridpack::analysis::NetworkAnalytics
                       <gridpack::powerflow::PFNetwork>(pf_network));
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
   
    ds_network.reset(new gridpack::dynamic_simulation::DSFullNetwork(this->communicator()));
	//pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
    //  gridpack::dynamic_simulation::DSFullBranch>(ds_network);
}

/**
 * solve power flow before run dynamic simulation, return the flag whether the power flow is solved successfully
 */
bool gridpack::hadrec::HADRECAppModule::solvePowerFlowBeforeDynSimu_withFlag(
    const char *inputfile, int pfcase_idx){
	
	gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
    //t_total = timer->createCategory("Dynamic Simulation: Total Application");
    //timer->start(t_total);

    // read configuration file 
    t_config = timer->createCategory("Dynamic Simulation: Config");
    timer->start(t_config);
	
    //gridpack::utility::Configuration *config =
    //  gridpack::utility::Configuration::configuration();
	  
    //config_sptr.reset(gridpack::utility::Configuration::configuration());
	//config_sptr.reset( new gridpack::utility::Configuration() )
	
	if (!bconfig_sptr_set){
		if (inputfile) {
                  config_sptr->open(inputfile, this->communicator());
		} else {
                  config_sptr->open("input.xml",this->communicator());
		}
		bconfig_sptr_set = true;
	}
    timer->stop(t_config);

    // setup and run powerflow calculation
    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config_sptr->getCursor("Configuration.Powerflow");
    bool useNonLinear = false;
    useNonLinear = cursor->get("UseNonLinear", useNonLinear);
	
    pf_network.reset(new gridpack::powerflow::PFNetwork(this->communicator()));

    pf_app_sptr.reset(new gridpack::powerflow::PFAppModule()) ;
	
	bool bnoprint = gridpack::NoPrint::instance()->status();

	if (bnoprint) { //if no print flag set to be true, suppress any intermedium output from pf app module
		pf_app_sptr->suppressOutput(true);
	}
	
    pf_app_sptr->readNetwork(pf_network, &(*config_sptr), pfcase_idx);
    pf_app_sptr->initialize();
    pf_analytics.reset(new gridpack::analysis::NetworkAnalytics<gridpack::powerflow::PFNetwork>(pf_network));
	
	bool pf_succ_flag;
    if (useNonLinear) {
      pf_succ_flag = pf_app_sptr->nl_solve();
    } else {
      pf_succ_flag = pf_app_sptr->solve();
    }
	
	if (!bnoprint) { //if no print flag set to be true, do not print the power flow solution
		pf_app_sptr->write();
	}
	
    pf_app_sptr->saveData();
	
    // setup dynamic simulation network
   
    ds_network.reset(new gridpack::dynamic_simulation::DSFullNetwork(this->communicator()));
	//pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
    //  gridpack::dynamic_simulation::DSFullBranch>(ds_network);
	return pf_succ_flag;
	
}

/**
 * read in power flow data 
 */
void gridpack::hadrec::HADRECAppModule::readPowerFlowData(
    const char *inputfile, int pfcase_idx){
			
	gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
    //t_total = timer->createCategory("Dynamic Simulation: Total Application");
    //timer->start(t_total);

    // read configuration file 
    t_config = timer->createCategory("Dynamic Simulation: Config");
    timer->start(t_config);
	
    //gridpack::utility::Configuration *config =
    //  gridpack::utility::Configuration::configuration();
	  
    //config_sptr.reset(gridpack::utility::Configuration::configuration());
	//config_sptr.reset( new gridpack::utility::Configuration() )
	
	if (!bconfig_sptr_set){
		if (inputfile) {
			config_sptr->open(inputfile, this->communicator());
		} else {
			config_sptr->open("input.xml",this->communicator());
		}
		bconfig_sptr_set = true;
	}
    timer->stop(t_config);

    // setup and run powerflow calculation
    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = config_sptr->getCursor("Configuration.Powerflow");
    //bool useNonLinear = false;
    p_PFuseNonLinear = cursor->get("UseNonLinear", p_PFuseNonLinear);
	
    pf_network.reset(new gridpack::powerflow::PFNetwork(this->communicator()));

    pf_app_sptr.reset(new gridpack::powerflow::PFAppModule()) ;
	
	bool bnoprint = gridpack::NoPrint::instance()->status();

	if (bnoprint) { //if no print flag set to be true, suppress any intermedium output from pf app module
		pf_app_sptr->suppressOutput(true);
	}
	
    pf_app_sptr->readNetwork(pf_network, &(*config_sptr), pfcase_idx);
	
	// setup dynamic simulation network
   
    ds_network.reset(new gridpack::dynamic_simulation::DSFullNetwork(this->communicator()));
	//pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
    //  gridpack::dynamic_simulation::DSFullBranch>(ds_network);
}

/**
 * solve power flow
 */
bool gridpack::hadrec::HADRECAppModule::solvePowerFlow(){
	
	pf_app_sptr->initialize();
        pf_analytics.reset(new gridpack::analysis::NetworkAnalytics<gridpack::powerflow::PFNetwork>(pf_network));
	
	bool pf_succ_flag;
    if (p_PFuseNonLinear) {
      pf_succ_flag = pf_app_sptr->nl_solve();
    } else {
      pf_succ_flag = pf_app_sptr->solve();
    }
	
	bool bnoprint = gridpack::NoPrint::instance()->status();
	
	if (!bnoprint) { //if no print flag set to be true, do not print the power flow solution
		pf_app_sptr->write();
	}
	
    pf_app_sptr->saveDataAlsotoOrg();
	
	return pf_succ_flag;
	
}

/**
 * transfer data from power flow to dynamic simulation 
 */
void gridpack::hadrec::HADRECAppModule::transferPFtoDS()
{	  
  pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
		    gridpack::dynamic_simulation::DSFullBranch>(ds_network);

  ds_app_sptr->transferPFtoDS(pf_network,ds_network);
  ds_analytics.reset(new gridpack::analysis::NetworkAnalytics
                     <gridpack::dynamic_simulation::DSFullNetwork>(ds_network));
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

  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config_sptr->getCursor("Configuration.Dynamic_simulation");
  
  if (faults.empty()){
	faults = ds_app_sptr->getEvents(cursor);
  }

  // run dynamic simulation
  ds_app_sptr->setNetwork(ds_network, &(*config_sptr));
  //ds_app_sptr->readNetwork(ds_network,config);
  
  //printf("ds_app_sptr->readGenerators:\n");
  ds_app_sptr->readGenerators(dscase_idx);
  //printf("ds_app_sptr->initialize:\n");
  ds_app_sptr->initialize();
  ds_analytics.reset(new gridpack::analysis::NetworkAnalytics
                     <gridpack::dynamic_simulation::DSFullNetwork>(ds_network));
  //printf("ds_app_sptr->setGeneratorWatch:\n");
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
	
	std::vector<double> rSpd, rAng, genP, genQ, vMag, vAng, fOnline; 
	//ds_app_sptr->getObservations(vMag, vAng, rSpd, rAng, genP, genQ, fOnline);
	
	std::vector<double> busfreq;
	ds_app_sptr->getObservations_withBusFreq(vMag, vAng, rSpd, rAng, genP, genQ, fOnline, busfreq);
	
	std::vector<double> obs_vals;
	obs_vals.clear();
	
	int i;
	for (i=0; i<rSpd.size(); i++){
		obs_vals.push_back(rSpd[i]);
	}
	
	for (i=0; i<rAng.size(); i++){
		obs_vals.push_back(rAng[i]);
	}
	
	for (i=0; i<genP.size(); i++){
		obs_vals.push_back(genP[i]);
	}
	
	for (i=0; i<genQ.size(); i++){
		obs_vals.push_back(genQ[i]);
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
	
	for (i=0; i<busfreq.size(); i++){
		obs_vals.push_back(busfreq[i]);
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
 * return observations list with bus frequency as observations
 */

void gridpack::hadrec::HADRECAppModule::getObservationLists_withBusFreq(
    std::vector<int> &genBuses, std::vector<std::string> &genIDs,
    std::vector<int> &loadBuses, std::vector<std::string> &loadIDs,
    std::vector<int> &busIDs, std::vector<int> &busfreqIDs){
	
	ds_app_sptr->getObservationLists_withBusFreq(genBuses, genIDs, loadBuses,
       loadIDs, busIDs, busfreqIDs);
	
}

/**
 * Check whether the dynamic simulation is done
 */
bool gridpack::hadrec::HADRECAppModule::isDynSimuDone( ){
	
	return ds_app_sptr->isDynSimuDone();

}

/**
 * Return values for total active and reactive load power on bus
 * @param bus_id original bus index
 * @param lp active load power
 * @param lq reactive load power
 * @return false if bus is not found on this processor
 */
bool gridpack::hadrec::HADRECAppModule::getBusTotalLoadPower(int bus_id,
    double &total_p, double &total_q)
{
	double lp = 0.0;
	double lq = 0.0;
	bool flag = false;
	flag = ds_app_sptr->getBusTotalLoadPower(bus_id, lp, lq);
    total_p = lp;
	total_q = lq;
	
	return flag;
}

/**
 * Return real and reactive power produced by requested generator
 * @param bus_id original index for bus hosting generator
 * @param gen_id 2-character identifier for generator
 * @param pg active power produced by generator
 * @param qg reactive power produced by generator
 * @return false if generator is not found on this processor
 */
bool gridpack::hadrec::HADRECAppModule::getGeneratorPower(int bus_id,
    std::string gen_id, double &pg, double &qg)
{ 
  	double pgtmp = 0.0;
	double qgtmp = 0.0;
	bool flag = false;
	flag = ds_app_sptr->getGeneratorPower(bus_id, gen_id, pgtmp, qgtmp);
    pg = pgtmp;
	qg = qgtmp;
	
	return flag;
}

/**
 * Return total active and reactive loads for each zone
 * @param load_p active load for all zones
 * @param load_q reactive load for all zones
 * @param zone_id label for each zone
 */
bool gridpack::hadrec::HADRECAppModule::getZoneLoads(std::vector<double> &load_p, std::vector<double> &load_q,
    std::vector<int> &zone_id) const
{
	ds_app_sptr->getZoneLoads(load_p, load_q, zone_id);
	
	return true;
}

/**
 * Return total active and reactive generator power for each zone
 * @param generator_p active generator power for all zones
 * @param generator_q reactive generator power for all zones
 * @param zone_id label for each zone
 */
bool gridpack::hadrec::HADRECAppModule::getZoneGeneratorPower(std::vector<double> &generator_p,
    std::vector<double> &generator_q, std::vector<int> &zone_id) const
{
	ds_app_sptr->getZoneGeneratorPower(generator_p, generator_q, zone_id);
	
	return true;
}

/**
 * apply actions
 */
void gridpack::hadrec::HADRECAppModule::applyAction(gridpack::hadrec::HADRECAction control_action){
	
	if ( control_action.actiontype == 0 ){
		return ds_app_sptr->applyLoadShedding(control_action.bus_number, control_action.componentID, control_action.percentage );
	}
	
	if ( control_action.actiontype == 2 ){
		return ds_app_sptr->applyGeneratorTripping(control_action.bus_number, control_action.componentID );
	}
	
	if ( control_action.actiontype == 3 ){
		return ds_app_sptr->applyConstYLoad_Change_P(control_action.bus_number, control_action.percentage );
	}
	
	if ( control_action.actiontype == 4 ){
		return ds_app_sptr->applyConstYLoad_Change_Q(control_action.bus_number, control_action.percentage );
	}
	
	if ( control_action.actiontype == 5 || control_action.actiontype == 6 || control_action.actiontype == 7 || control_action.actiontype == 8){
		return ds_app_sptr->applyGFIAdjustment(control_action.actiontype-5, control_action.bus_number, control_action.componentID, control_action.percentage );
	}
	
	if ( control_action.actiontype == 9 ){
		return ds_app_sptr->applyConstYLoadShedding(control_action.bus_number, control_action.percentage );
	}
	
	if ( control_action.actiontype == 1 ){
		if (control_action.brch_from_bus_number > 0 || control_action.brch_to_bus_number > 0){
			return ds_app_sptr->setLineTripAction(control_action.brch_from_bus_number, control_action.brch_to_bus_number, control_action.branch_ckt);
		}
		if (control_action.brch_from_bus_number == -1 && control_action.brch_to_bus_number == -1 && control_action.bus_number > 0){
		
			return ds_app_sptr->setLineTripAction(control_action.bus_number);
		}
	}

}

/**
 * set the wide area control signals of the PSS of a certain generator
 * input bus_number: generator bus number
 * input bus_number: generator gen ID
 * input wideAreaControlSignal:  wide area control signal for the PSS of the generator
 */
void gridpack::hadrec::HADRECAppModule::setWideAreaControlSignal(int bus_number, std::string genid, double wideAreaControlSignal){
	
	ds_app_sptr->setWideAreaControlSignal(bus_number, genid, wideAreaControlSignal);
	
}

/**
 * execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
 * the vector  vloadP and vloadQ
 */
void gridpack::hadrec::HADRECAppModule::scatterInjectionLoad(const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ){
	
	ds_app_sptr->scatterInjectionLoad(vbusNum, vloadP, vloadQ);
	
}

/**
* execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
* the vector  vloadP and vloadQ - new implemnetation by removing the contribution of the original constant Y load from y-maxtrix, 
* and model the entire load change as injection current
*/
void gridpack::hadrec::HADRECAppModule::scatterInjectionLoadNew(const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ){
	
	ds_app_sptr->scatterInjectionLoadNew(vbusNum, vloadP, vloadQ);
	
}

/**
* execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
* the vector  vloadP and vloadQ - this implemnetation keeps the Y load component of the bus still at the bus, 
* while only compenstates the difference
*/
void gridpack::hadrec::HADRECAppModule::scatterInjectionLoadNew_compensateY(const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ){
	
	ds_app_sptr->scatterInjectionLoadNew_compensateY(vbusNum, vloadP, vloadQ);
	
}

/**
 * execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
 * the vector  vloadP and vloadQ - new implemnetation by removing the contribution of the original constant Y load from y-maxtrix, 
 * and model the entire load change as injection current, also add a large parallel connecting impedance at the certain bus
*/
void gridpack::hadrec::HADRECAppModule::scatterInjectionLoadNew_Norton(const std::vector<int>& vbusNum, 
							const std::vector<double>& vloadP, const std::vector<double>& vloadQ, 
							const std::vector<double>& vimpedanceR, const std::vector<double>& vimpedanceI){
								
	ds_app_sptr->scatterInjectionLoadNew_Norton(vbusNum, vloadP, vloadQ, vimpedanceR, vimpedanceI);					
}

/**
* execute load scattering, the values of the STATIC load current at certain buses vbusNum will be changed to the values of 
* the vector  vCurR and vCurI - new implemnetation by removing the contribution of the original constant Y load from y-maxtrix, 
* and model the entire load change as injection current
*/
void gridpack::hadrec::HADRECAppModule::scatterInjectionLoadNewConstCur(const std::vector<int>& vbusNum, const std::vector<double>& vCurR, const std::vector<double>& vCurI){
	
	ds_app_sptr->scatterInjectionLoadNewConstCur(vbusNum, vCurR, vCurI);
	
}

/**
 * get the power flow solution for the specific bus, vmag and v angle
 * @param bus original number, bus solution vmag and v angle
 * @return false if location of bus is not found in
 * network
 */

bool gridpack::hadrec::HADRECAppModule::getPFSolutionSingleBus(
    int bus_number, double &bus_mag, double &bus_angle){
	
	double Vmag = 1.0;
	double Vangle = 0.0;
	
	bool ret = pf_app_sptr->getPFSolutionSingleBus(bus_number, Vmag, Vangle);
	
	bus_mag = Vmag;
	bus_angle = Vangle;
	
	return ret;
}


/**
 * Modify generator parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param gen_id two character token specifying generator on bus
 * @param genParam string representing dictionary name of data element
 *                to be modified
 * @param value new value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::hadrec::HADRECAppModule::modifyDataCollectionGenParam(
    int bus_id, std::string gen_id, std::string genParam, double value)
{
    return pf_app_sptr->modifyDataCollectionGenParam(bus_id, gen_id,genParam,value);
}
bool gridpack::hadrec::HADRECAppModule::modifyDataCollectionGenParam(
    int bus_id, std::string gen_id, std::string genParam, int value)
{
    return pf_app_sptr->modifyDataCollectionGenParam(bus_id, gen_id,genParam,value);
}

/**
 * Modify load parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param load_id two character token specifying load on bus
 * @param loadParam string representing dictionary name of data element
 *                to be modified
 * @param value new value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::hadrec::HADRECAppModule::modifyDataCollectionLoadParam(
    int bus_id, std::string load_id, std::string loadParam, double value)
{
    return pf_app_sptr->modifyDataCollectionLoadParam(bus_id, load_id,loadParam,value);
}
bool gridpack::hadrec::HADRECAppModule::modifyDataCollectionLoadParam(
    int bus_id, std::string load_id, std::string loadParam, int value)
{
    return pf_app_sptr->modifyDataCollectionLoadParam(bus_id, load_id,loadParam,value);
}

/**
 * Modify parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param busParam string representing dictionary name of data element
 *                to be modified
 * @param value new value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::hadrec::HADRECAppModule::modifyDataCollectionBusParam(
    int bus_id, std::string busParam, double value)
{
  return pf_app_sptr->modifyDataCollectionBusParam(bus_id, busParam,value);
}
bool gridpack::hadrec::HADRECAppModule::modifyDataCollectionBusParam(
    int bus_id, std::string busParam, int value)
{
  return pf_app_sptr->modifyDataCollectionBusParam(bus_id, busParam,value);
}

/**
 * Modify parameters in data collection for specified branch
 * @param bus1, bus2 bus IDs for from and to bus
 * @param ckt two character token specifying branch
 * @param branchParam string representing dictionary name of data element
 *                to be modified
 * @param value new value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::hadrec::HADRECAppModule::modifyDataCollectionBranchParam(
    int bus1, int bus2, std::string ckt,
    std::string branchParam, double value)
{
  return pf_app_sptr->modifyDataCollectionBranchParam(bus1,bus2,ckt,branchParam,value);
}
bool gridpack::hadrec::HADRECAppModule::modifyDataCollectionBranchParam(
    int bus1, int bus2, std::string ckt,
    std::string branchParam, int value)
{
  return pf_app_sptr->modifyDataCollectionBranchParam(bus1,bus2,ckt,branchParam,value);
}

/**
 * Get generator parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param gen_id two character token specifying generator on bus
 * @param genParam string representing dictionary name of data element
 *                to be modified
 * @param value value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::hadrec::HADRECAppModule::getDataCollectionGenParam(
    int bus_id, std::string gen_id,
    std::string genParam, double &value)
{ 
  double *pvalue = &value;
  return pf_app_sptr->getDataCollectionGenParam(bus_id, gen_id, genParam, pvalue);
}
bool gridpack::hadrec::HADRECAppModule::getDataCollectionGenParam(
    int bus_id, std::string gen_id,
    std::string genParam, int &value)
{
	int *pvalue = &value;
  return pf_app_sptr->getDataCollectionGenParam(bus_id, gen_id, genParam, pvalue);
}

/**
 * Get load parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param load_id two character token specifying load on bus
 * @param loadParam string representing dictionary name of data element
 *                to be modified
 * @param value value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::hadrec::HADRECAppModule::getDataCollectionLoadParam(
    int bus_id, std::string load_id,
    std::string loadParam, double &value)
{
  double *pvalue = &value;
  return pf_app_sptr->getDataCollectionLoadParam(bus_id, load_id, loadParam, pvalue);
}
bool gridpack::hadrec::HADRECAppModule::getDataCollectionLoadParam(
    int bus_id, std::string load_id,
    std::string loadParam, int &value)
{
	int *pvalue = &value;
  return pf_app_sptr->getDataCollectionLoadParam(bus_id, load_id, loadParam, pvalue);
}

/**
 * Get parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param busParam string representing dictionary name of data element
 *                to be modified
 * @param value value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::hadrec::HADRECAppModule::getDataCollectionBusParam(
    int bus_id, std::string busParam, double &value)
{
	double *pvalue = &value;
  return pf_app_sptr->getDataCollectionBusParam(bus_id, busParam, pvalue);
}
bool gridpack::hadrec::HADRECAppModule::getDataCollectionBusParam(
    int bus_id, std::string busParam, int &value)
{
	int *pvalue = &value;
  return pf_app_sptr->getDataCollectionBusParam(bus_id, busParam, pvalue);
}

/**
 * Get parameters in data collection for specified branch
 * @param bus1, bus2 bus IDs for from and to bus
 * @param ckt two character token specifying branch
 * @param branchParam string representing dictionary name of data element
 *                to be modified
 * @param value value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::hadrec::HADRECAppModule::getDataCollectionBranchParam(
    int bus1, int bus2, std::string ckt,
    std::string branchParam, double &value)
{
	double *pvalue = &value;
  return pf_app_sptr->getDataCollectionBranchParam(bus1, bus2, ckt, branchParam, pvalue);
}
bool gridpack::hadrec::HADRECAppModule::getDataCollectionBranchParam(
    int bus1, int bus2, std::string ckt,
    std::string branchParam, int &value)
{
	int *pvalue = &value;
  return pf_app_sptr->getDataCollectionBranchParam(bus1, bus2, ckt, branchParam, pvalue);
}

/**
 * Export final solved power flow to PSS/E formatted file, version 23
 * @param filename name of file to store network configuration
 */
void gridpack::hadrec::HADRECAppModule::exportPSSE23(std::string filename)
{
  pf_app_sptr->exportPSSE23(filename);
}

/**
 * Export final solved power flow to PSS/E formatted file, version 33
 * @param filename name of file to store network configuration
 */
void gridpack::hadrec::HADRECAppModule::exportPSSE33(std::string filename)
{
  pf_app_sptr->exportPSSE33(filename);
}

/**
 * Export final solved power flow to PSS/E formatted file, version 34
 * @param filename name of file to store network configuration
 */
void gridpack::hadrec::HADRECAppModule::exportPSSE34(std::string filename)
{
  pf_app_sptr->exportPSSE34(filename);
}

/**
 * Set the state of some device on the network
 * @param bus_id bus ID
 * @param dev_id two character identifier of device
 * @param type of device (GENERATOR, EXCITER, OR GOVERNOR)
 * @param name of the state variable (ANGLE,SPEED_DEV)
 * @param value new value of parameter
 * @return false if this device or
 parameter not found
 */
bool gridpack::hadrec::HADRECAppModule::setState(int bus_id,
    std::string dev_id, std::string device,
    std::string name, double value)
{
  return ds_app_sptr->setState(bus_id, dev_id, device, name, value);
}

/**
 * Get the state of some device on the network
 * @param  bus_id bus ID
 * @param  dev_id two character identifier of device
 * @param  type of device (GENERATOR, EXCITER, OR GOVERNOR)
 * @param  name of the state variable (ANGLE, SPEED_DEV)
 * @return value current value of parameter
 * @return false if this device or state variable not found
 */
bool gridpack::hadrec::HADRECAppModule::getState(int bus_id,
    std::string dev_id, std::string device,
    std::string name, double *value)
{
  return ds_app_sptr->getState(bus_id, dev_id, device, name, value);
}


int gridpack::hadrec::HADRECAppModule::totalBuses(void) const
{
  int result(0);

  if (ds_analytics) {
    result = ds_analytics->totalBuses();
  } else if (pf_analytics) {
    result = pf_analytics->totalBuses();
  } else {
    throw gridpack::Exception("HADRECAppModule::totalBuses(): network not defined");
  }

  return result;
}


int gridpack::hadrec::HADRECAppModule::totalBranches(void) const
{
  int result(0);

  if (ds_analytics) {
    result = ds_analytics->totalBranches();
  } else if (pf_analytics) {
    result = pf_analytics->totalBranches();
  } else {
    throw gridpack::Exception("HADRECAppModule::totalBranches(): network not defined");
  }

  return result;
}


std::vector<int>
gridpack::hadrec::HADRECAppModule::getConnectedBranches(int oidx) const
{
  std::vector<int> result;
  if (ds_analytics) {
    result = ds_analytics->getConnectedBranches(oidx);
  } else if (pf_analytics) {
    result = pf_analytics->getConnectedBranches(oidx);
  } else {
    throw gridpack::Exception("HADRECAppModule::totalBranches(): network not defined");
  }
  
  return result;
}

void
gridpack::hadrec::HADRECAppModule::getBranchEndpoints(const int& idx, int *fbus, int *tbus) const
{
  if (ds_analytics) {
    ds_analytics->getBranchEndpoints(idx, fbus, tbus);
  } else if (pf_analytics) {
    pf_analytics->getBranchEndpoints(idx, fbus, tbus);
  } else {
    throw gridpack::Exception("HADRECAppModule::getBranchEndpoints(): network not defined");
  }
}

int gridpack::hadrec::HADRECAppModule::numGenerators(void) const
{
  int result(0);

  if (ds_analytics) {
    result = ds_analytics->numGenerators();
  } else if (pf_analytics) {
    result = pf_analytics->numGenerators();
  } else {
    throw gridpack::Exception("HADRECAppModule::numGenerators(): network not defined");
  }

  return result;
}

int gridpack::hadrec::HADRECAppModule::numLoads(void) const
{
  int result(0);

  if (ds_analytics) {
    result = ds_analytics->numLoads();
  } else if (pf_analytics) {
    result = pf_analytics->numLoads();
  } else {
    throw gridpack::Exception("HADRECAppModule::numLoads(): network not defined");
  }

  return result;
}

int gridpack::hadrec::HADRECAppModule::numStorage(void) const
{
  int result(0);

  if (ds_analytics) {
    result = ds_analytics->numStorage();
  } else if (pf_analytics) {
    result = pf_analytics->numStorage();
  } else {
    throw gridpack::Exception("HADRECAppModule::numStorage(): network not defined");
  }

  return result;
}
