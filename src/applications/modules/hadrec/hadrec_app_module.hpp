/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hadrec_app_module.hpp
 * @author Bruce Palmer
 * @date   2020-04-16 10:52:34 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _hadrec_app_module_h_
#define _hadrec_app_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "gridpack/applications/modules/dynamic_simulation_full_y/dsf_app_module.hpp"

namespace gridpack {
namespace hadrec {

struct HADRECAction
{
	int actiontype = 0;  // 0: load shedding
	int bus_number = -1;
	std::string componentID = "1";
	double percentage = 0.0;
	//p_obs_vals.clear();
	
};


class HADRECAppModule
{
  public:
    /**
     * Basic constructor
     */
    HADRECAppModule(void);

    /**
     * Basic destructor
     */
    ~HADRECAppModule(void);
	
	/**
	 * solve power flow before run dynamic simulation 
	 */
	void solvePowerFlowBeforeDynSimu(const char *inputfile, int pfcase_idx = -1);
	
	/**
	* transfer data from power flow to dynamic simulation 
	*/
	void transferPFtoDS();
	
	/**
    * do initialization only for dynamics simulation
    * @param faults lists of faults that might be simulated
    * @param dscase_idx index pointing to dyr parameter file that should be used if
    *                   a list of files is supplied in input deck
     */
    void initializeDynSimu(std::vector<gridpack::dynamic_simulation::Event> faults,
        int dscase_idx=-1);
	
	/**
	* do a fully initialization before running dynamics simulation
   * @param case_idx index pointing to network configuration and dyr parameter
   * file that should be used if a list of files is supplied in input deck
	*/
	void fullInitializationBeforeDynSimuSteps(const char *inputfile,
           const std::vector<gridpack::dynamic_simulation::Event>& BusFaults,
           int pfcase_idx=-1, int dscase_idx=-1);
	
	/**
	* Execute only one simulation time step 
	*/
	void executeDynSimuOneStep();
	
	/**
	* Check whether the dynamic simulation is done
	*/
	bool isDynSimuDone( );
	
	/**
	* apply actions
	*/
	void applyAction(gridpack::hadrec::HADRECAction control_action);
	
	/**
	* return observations after each simulation time step
	*/
	std::vector<double> getObservations();
	
	/**
	* return observations list
	*/
	void getObservationLists(std::vector<int> &genBuses,
       std::vector<std::string> &genIDs, std::vector<int> &loadBuses,
       std::vector<std::string> &loadIDs, std::vector<int> &busIDs);
	

  private:
   boost::shared_ptr<gridpack::utility::Configuration> config_sptr;
	boost::shared_ptr<gridpack::powerflow::PFNetwork> pf_network;
	boost::shared_ptr<gridpack::powerflow::PFAppModule> pf_app_sptr;
	
	boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork> ds_network;
	boost::shared_ptr<gridpack::dynamic_simulation::DSFullApp> ds_app_sptr;
	
    int t_total;
	int t_config;

   // Observations
   std::vector<int> p_obs_genBus;
   std::vector<std::string> p_obs_genIDs;
   std::vector<int> p_obs_loadBus;
   std::vector<std::string> p_obs_loadIDs;
   std::vector<int> p_obs_vBus;
   //std::vector<double> p_obs_vals;
	
};

} // hadrec
} // gridpack
#endif
