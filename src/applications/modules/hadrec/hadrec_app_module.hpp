/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hadrec_app_module.hpp
 * @author Bruce Palmer
 * @date   2023-10-25 08:28:44 d3g096
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
	int actiontype = 0;  // 0: load shedding; 1: line tripping; 2: generator tripping; 3: Constant-Y load P change; 4: Constant-Y load Q change
						 // 5: GFI mp adjust; 6: GFI mq adjust; 7: GFI Pset adjust; 8: GFI Qset adjust
	int bus_number = -1;
	std::string componentID = "1";
	double percentage = 0.0;
	
	int brch_from_bus_number = -1;
	int brch_to_bus_number = -1;
	std::string branch_ckt = "1";
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
	 * solve power flow before run dynamic simulation, return the flag whether the power flow is solved successfully
	 */
	bool solvePowerFlowBeforeDynSimu_withFlag(const char *inputfile, int pfcase_idx = -1);
	
    /**
	 * read in power flow data 
	 */
	void readPowerFlowData(const char *inputfile, int pfcase_idx);
	
	/**
	 * solve power flow
	 */
	bool solvePowerFlow();
	
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
	 * set the wide area control signals of the PSS of a certain generator
	 * input bus_number: generator bus number
	 * input bus_number: generator gen ID
	 * input wideAreaControlSignal:  wide area control signal for the PSS of the generator
	 */
	void setWideAreaControlSignal(int bus_number, std::string genid, double wideAreaControlSignal);
	
	/**
	 * execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
	 * the vector  vloadP and vloadQ
	*/
	void scatterInjectionLoad(const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ);
	
	/**
	* execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
	* the vector  vloadP and vloadQ - new implemnetation by removing the contribution of the original constant Y load from y-maxtrix, 
	* and model the entire load change as injection current
	*/
 
	void scatterInjectionLoadNew(const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ);
	
	/**
	* execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
	* the vector  vloadP and vloadQ - this implemnetation keeps the Y load component of the bus still at the bus, 
	* while only compenstates the difference
	*/
 
	void scatterInjectionLoadNew_compensateY(const std::vector<int>& vbusNum, const std::vector<double>& vloadP, const std::vector<double>& vloadQ);
	
	/**
	 * execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
	 * the vector  vloadP and vloadQ - new implemnetation by removing the contribution of the original constant Y load from y-maxtrix, 
	 * and model the entire load change as injection current, also add a large parallel connecting impedance at the certain bus
	 */
	 void scatterInjectionLoadNew_Norton(const std::vector<int>& vbusNum, 
							const std::vector<double>& vloadP, const std::vector<double>& vloadQ, 
							const std::vector<double>& vimpedanceR, const std::vector<double>& vimpedanceI);
	
	/**
	* execute load scattering, the values of the STATIC load current at certain buses vbusNum will be changed to the values of 
	* the vector  vCurR and vCurI - new implemnetation by removing the contribution of the original constant Y load from y-maxtrix, 
	* and model the entire load change as injection current
	*/
	void scatterInjectionLoadNewConstCur(const std::vector<int>& vbusNum, const std::vector<double>& vCurR, const std::vector<double>& vCurI);
	
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
	   
	/**
	 * return observations list with bus frequency as observations
	 */

	void getObservationLists_withBusFreq(
		std::vector<int> &genBuses, std::vector<std::string> &genIDs,
		std::vector<int> &loadBuses, std::vector<std::string> &loadIDs,
		std::vector<int> &busIDs, std::vector<int> &busfreqIDs);
	   
	/**
     * Return values for total active and reactive load power on bus
     * @param bus_id original bus index
     * @param lp active load power
     * @param lq reactive load power
     * @return false if bus is not found on this processor
     */
    bool getBusTotalLoadPower(int bus_id, double &total_p, double &total_q);

    /**
     * Return real and reactive power produced by requested generator
     * @param bus_id original index for bus hosting generator
     * @param gen_id 2-character identifier for generator
     * @param pg active power produced by generator
     * @param qg reactive power produced by generator
     * @return false if generator is not found on this processor
     */
    bool getGeneratorPower(int bus_id, std::string gen_id, double &pg, double &qg);
	
	    /**
     * Return total active and reactive loads for each zone
     * @param load_p active load for all zones
     * @param load_q reactive load for all zones
     * @param zone_id label for each zone
     */
    bool getZoneLoads(std::vector<double> &load_p, std::vector<double> &load_q,
        std::vector<int> &zone_id) const;

    /**
     * Return total active and reactive generator power for each zone
     * @param generator_p active generator power for all zones
     * @param generator_q reactive generator power for all zones
     * @param zone_id label for each zone
     */
    bool getZoneGeneratorPower(std::vector<double> &generator_p,
        std::vector<double> &generator_q, std::vector<int> &zone_id) const;
		
	/**
	* get the power flow solution for the specific bus, vmag and v angle
	* @param bus original number, bus solution vmag and v angle
	* @return false if location of bus is not found in
	* network
	*/

	bool getPFSolutionSingleBus(int bus_number, double &bus_mag, double &bus_angle);
	
	/**
     * Modify generator parameters in data collection for specified bus
     * @param bus_id bus ID
     * @param gen_id two character token specifying generator on bus
     * @param genParam string representing dictionary name of data element
     *                to be modified
     * @param value new value of parameter
     * @return return false if parameter is not found
     */
    bool modifyDataCollectionGenParam(int bus_id, std::string gen_id,
        std::string genParam, double value);
    bool modifyDataCollectionGenParam(int bus_id, std::string gen_id,
        std::string genParam, int value);

    /**
     * Modify load parameters in data collection for specified bus
     * @param bus_id bus ID
     * @param load_id two character token specifying load on bus
     * @param loadParam string representing dictionary name of data element
     *                to be modified
     * @param value new value of parameter
     * @return return false if parameter is not found
     */
    bool modifyDataCollectionLoadParam(int bus_id, std::string load_id,
        std::string loadParam, double value);
    bool modifyDataCollectionLoadParam(int bus_id, std::string load_id,
        std::string loadParam, int value);

    /**
     * Modify parameters in data collection for specified bus
     * @param bus_id bus ID
     * @param busParam string representing dictionary name of data element
     *                to be modified
     * @param value new value of parameter
     * @return return false if parameter is not found
     */
    bool modifyDataCollectionBusParam(int bus_id,
        std::string busParam, double value);
    bool modifyDataCollectionBusParam(int bus_id,
        std::string busParam, int value);
		
	    /**
     * Modify parameters in data collection for specified branch
     * @param bus1, bus2 bus IDs for from and to bus
     * @param ckt two character token specifying branch
     * @param branchParam string representing dictionary name of data element
     *                to be modified
     * @param value new value of parameter
     * @return return false if parameter is not found
     */
    bool modifyDataCollectionBranchParam(int bus1, int bus2, std::string ckt,
        std::string branchParam, double value);
    bool modifyDataCollectionBranchParam(int bus1, int bus2, std::string ckt,
        std::string branchParam, int value);

	/**
     * Get generator parameters in data collection for specified bus
     * @param bus_id bus ID
     * @param gen_id two character token specifying generator on bus
     * @param genParam string representing dictionary name of data element
     *                to be modified
     * @param value value of parameter
     * @return return false if parameter is not found
     */
    bool getDataCollectionGenParam(int bus_id, std::string gen_id,
        std::string genParam, double &value);
    bool getDataCollectionGenParam(int bus_id, std::string gen_id,
        std::string genParam, int &value);

    /**
     * Get load parameters in data collection for specified bus
     * @param bus_id bus ID
     * @param load_id two character token specifying load on bus
     * @param loadParam string representing dictionary name of data element
     *                to be modified
     * @param value value of parameter
     * @return return false if parameter is not found
     */
    bool getDataCollectionLoadParam(int bus_id, std::string load_id,
        std::string loadParam, double &value);
    bool getDataCollectionLoadParam(int bus_id, std::string load_id,
        std::string loadParam, int &value);

    /**
     * Get parameters in data collection for specified bus
     * @param bus_id bus ID
     * @param busParam string representing dictionary name of data element
     *                to be modified
     * @param value value of parameter
     * @return return false if parameter is not found
     */
    bool getDataCollectionBusParam(int bus_id,
        std::string busParam, double &value);
    bool getDataCollectionBusParam(int bus_id,
        std::string busParam, int &value);

    /**
     * Get parameters in data collection for specified branch
     * @param bus1, bus2 bus IDs for from and to bus
     * @param ckt two character token specifying branch
     * @param branchParam string representing dictionary name of data element
     *                to be modified
     * @param value value of parameter
     * @return return false if parameter is not found
     */
    bool getDataCollectionBranchParam(int bus1, int bus2, std::string ckt,
        std::string branchParam, double &value);
    bool getDataCollectionBranchParam(int bus1, int bus2, std::string ckt,
        std::string branchParam, int &value);
	
	/**
     * Export final solved power flow to PSS/E formatted file, version 23
     * @param filename name of file to store network configuration
     */
	void exportPSSE23(std::string filename);

  /**
   * Export final solved power flow to PSS/E formatted file, version 33
   * @param filename name of file to store network configuration
   */
  void exportPSSE33(std::string filename);

  /**
   * Export final solved power flow to PSS/E formatted file, version 34
   * @param filename name of file to store network configuration
   */
  void exportPSSE34(std::string filename);

   /**
    * Set the state of some device on the network
    * @param bus_id bus ID
    * @param dev_id two character identifier of device
    * @param device type of device to be modified
    * @param name string labeling parameter to be modified
    * @param value new value of parameter
    * @return false if this device or parameter not found
    */
   bool setState(int bus_id, std::string dev_id, std::string device,
       std::string name, double value);

   /**
    * Get the state of some device on the network
    * @param bus_id bus ID
    * @param dev_id two character identifier of device
    * @param device type of device to be modified
    * @param name string labeling parameter to be modified
    * @param value current value of parameter
    * @return false if this device or parameter not found
    */
   bool getState(int bus_id, std::string dev_id, std::string device,
       std::string name, double *value);

  private:
   boost::shared_ptr<gridpack::utility::Configuration> config_sptr;
	boost::shared_ptr<gridpack::powerflow::PFNetwork> pf_network;
	boost::shared_ptr<gridpack::powerflow::PFAppModule> pf_app_sptr;
	
	boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork> ds_network;
	boost::shared_ptr<gridpack::dynamic_simulation::DSFullApp> ds_app_sptr;
	
   gridpack::parallel::Communicator p_comm;

    int t_total;
	int t_config;
	
	bool bconfig_sptr_set;
	bool p_PFuseNonLinear;

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
