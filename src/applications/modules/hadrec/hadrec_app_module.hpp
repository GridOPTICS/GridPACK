/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hadrec_app_module.hpp
 * @author Bruce Palmer
 * @date   2014-01-28 11:30:49 d3g096
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

// Calling program for powerflow application

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
	void solvePowerFlowBeforeDynSimu(int argc, char **argv);
	
	/**
	* transfer data from power flow to dynamic simulation 
	*/
	void transferPFtoDS();
	
	/**
     * do initialization only for dynamics simulation
     */
    void initializeDynSimu();
	
	/**
	* do a fully initialization before running dynamics simulation
	*/
	void fullInitializationBeforeDynSimuSteps(int argc, char **argv);
	
	/**
	* Execute only one simulation time step 
	*/
	void executeDynSimuOneStep(
		std::vector<gridpack::dynamic_simulation::Event> action_list);
	
	/**
	* Check whether the dynamic simulation is done
	*/
	bool isDynSimuDone( );
	

  private:
    boost::shared_ptr<gridpack::utility::Configuration> config_sptr;
	boost::shared_ptr<gridpack::powerflow::PFNetwork> pf_network;
	boost::shared_ptr<gridpack::powerflow::PFAppModule> pf_app_sptr;
	
	boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork> ds_network;
	boost::shared_ptr<gridpack::dynamic_simulation::DSFullApp> ds_app_sptr;
	
    int t_total;
	int t_config;
	

};

} // hadrec
} // gridpack
#endif
