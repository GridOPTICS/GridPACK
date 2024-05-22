/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsf_app_module.hpp
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#ifndef _dsf_app_module_h_
#define _dsf_app_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/parallel/global_vector.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "dsf_factory.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/analysis/network_analytics.hpp"
//#include "gridpack/applications/modules/hadrec/hadrec_app_module.hpp"


namespace gridpack {
namespace dynamic_simulation {

    // Calling program for dynamic simulation application

class DSFullApp
{
  public:
    /**
     * Basic constructor
     */
    DSFullApp(void);

    /**
     * Basic constructor with commmunicator argument
     * @param comm communicator that application object is restricted to
     */
    DSFullApp(gridpack::parallel::Communicator comm);

    /**
     * Basic destructor
     */
    ~DSFullApp(void);

  /**
   * Solve power flow and use it to initialize dynamic
   * simulation. This creates and correctly handles the power flow and
   * dynamic simulation networks internally. No need to create them
   * externally.
   */
  void solvePowerFlowBeforeDynSimu(const char *inputfile, const int& pf_idx = -1);

    /**
     * Read in and partition the dynamic simulation network. The input file is read
     * directly from the Dynamic_simulation block in the configuration file so no
     * external file names or parameters need to be passed to this routine
     * @param network pointer to a DSFullNetwork object. This should not have any
     * buses or branches defined on it.
     * @param config pointer to open configuration file
     * @param otherfile name of network configuration file if different from the
     * one in the input deck
     */
    void readNetwork(boost::shared_ptr<DSFullNetwork> &network,
        gridpack::utility::Configuration *config,
        const char *otherfile = NULL);

    /**
     * Assume that DSFullNetwork already exists and just cache an internal pointer
     * to it. This routine does not call the partition function. Also read in
     * simulation parameters from configuration file
     * @param network pointer to a complete DSFullNetwork object.
     * @param config pointer to open configuration file
     */
    void setNetwork(boost::shared_ptr<DSFullNetwork> &network,
        gridpack::utility::Configuration *config);

    /**
     * Read generator parameters. These will come from a separate file (most
     * likely). The name of this file comes from the input configuration file.
     * @param ds_idx index of dyr file if a list of dyr files are provided.
     */
    void readGenerators(int ds_idx = -1);

    /**
     * Read sequence data from a file.
     */
    void readSequenceData();

    /**
     * Set up exchange buffers and other internal parameters and initialize
     * network components using data from data collection
     */
    void initialize();

    /**
     * Reinitialize calculation from data collections
     */
    void reload();

  /**
     * Reset data structures
     */
    void reset();


    /**
     * Execute the time integration portion of the application
     */
    void solve(gridpack::dynamic_simulation::Event fault);
	
	/**
     * initialization before the time step integration starts 
     */
    void solvePreInitialize(gridpack::dynamic_simulation::Event fault);

     /**
	Setup before the dynamic simulation begins
     **/
  void setup();
	/**
     * Execute only one simulation time step 
     */
    void executeOneSimuStep();

     /**
	Run upto a given time
     **/
  void run(double tend);

  /**
     Run till end time
  **/
  void run();

     
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
	void scatterInjectionLoadNew(const std::vector<int>& vbusNum, 
								const std::vector<double>& vloadP, const std::vector<double>& vloadQ);
								
	/**
	 * execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
	 * the vector  vloadP and vloadQ - this implemnetation keeps the Y load component of the bus still at the bus, 
	 * while only compenstates the difference
	*/
	void scatterInjectionLoadNew_compensateY(const std::vector<int>& vbusNum, 
								const std::vector<double>& vloadP, const std::vector<double>& vloadQ);
								
	/**
	 * execute load scattering, the P and Q values of the STATIC load at certain buses vbusNum will be changed to the values of 
	 * the vector  vloadP and vloadQ - new implemnetation by removing the contribution of the original constant Y load from y-maxtrix, 
     * and model the entire load change as injection current, also add a large parallel connecting impedance at the certain bus
	*/
	void scatterInjectionLoadNew_Norton(const std::vector<int>& vbusNum, 
								const std::vector<double>& vloadP, const std::vector<double>& vloadQ, 
								const std::vector<double>& vimpedanceR, const std::vector<double>& vimpedanceI);
								
	/**
	 * execute load scattering with constant current load , the values of the STATIC load current at certain buses vbusNum will be changed to the values of 
     * the vector  vCurR and vCurI - new implemnetation by removing the contribution of the original constant Y load from y-maxtrix, 
     * and model the entire load change as injection current
    */
 
    void scatterInjectionLoadNewConstCur(const std::vector<int>& vbusNum, 
										const std::vector<double>& vCurR, const std::vector<double>& vCurI);
	
	/**
     * execute load shedding	 
     */
    void applyLoadShedding(int bus_number, std::string loadid, double percentage);
	
	/**
	 * execute constant Y load shedding at a curtain bus	 
	 * bus number
	 * percentage: float load shed percentage, for example -0.2 means shed 20%
	 */
	void applyConstYLoadShedding(int bus_number, double percentage );
	
	/**
	 * set the wide area control signals of the PSS of a certain generator
	 * input bus_number: generator bus number
	 * input bus_number: generator gen ID
	 * input wideAreaControlSignal:  wide area control signal for the PSS of the generator
	 */
	void setWideAreaControlSignal(int bus_number, std::string genid, double wideAreaControlSignal);
	
	/**
     * execute Grid Forming Inverter control parameters adjustment
	 * input controlTyp: 0: GFI mp adjust; 1: GFI mq adjust; 2: GFI Pset adjust; 3: GFI Qset adjust; others: invalid
	 * input bus_number: GFI bus number
	 * input bus_number: GFI gen ID
	 * input newParValScaletoOrg: GFI new parameter scale value to the very intial value at the begining of dynamic simulation
     */
    void applyGFIAdjustment(int controlType, int bus_number, std::string genid, double newParValScaletoOrg);
	
	/**
     * execute Constant Y load P change	 
     */
    void applyConstYLoad_Change_P(int bus_number, double loadPChangeMW);
	void clearConstYLoad_Change_P();
	
	/**
     * execute Constant Y load Q change	 
     */
    void applyConstYLoad_Change_Q(int bus_number, double loadPChangeMVAR);
	void clearConstYLoad_Change_Q();
	
	void setConstYLoadtoZero_P(int bus_number);
	void setConstYLoadtoZero_Q(int bus_number);
	
	/**
     * set constant Y load to impedancer and impedancei
    */
    void setConstYLoadImpedance(int bus_number, double impedancer, double impedancei);
	
    /**
     * execute generator tripping 
     */
    void applyGeneratorTripping(int bus_number, std::string genid);
	
	/**
	* set all the necessery flags for the two buses and one branch for the line needs to trip
	* this function is for single branch flags set-up, may need to be called
	* multiple times for multiple line tripping 
	*/
	void setLineTripAction(int brch_from_bus_number, int brch_to_bus_number, std::string branch_ckt);
	
	/**
	* set all the necessery flags for the two buses and one branch for the line needs to trip
	* this function will trip a branch, given a bus number, just find any one of the 
	* connected line(not transformer) with the bus, and trip that one
	* this function is for single branch flags set-up, may need to be called
	* multiple times for multiple line tripping 
	*/
	void setLineTripAction(int bus_number);
	
	/**
	* clear all the necessery flags for the all buses and branches for the lines needs to trip
	* this function is for all the branches' flags clear-up, just need to be called
	* once to clear all the tripping lines's flag
	*/
	void clearLineTripAction();
	
	/**
     * Check whether the dynamic simulation is done
     */
    bool isDynSimuDone();

    /**
     * Write out final results of dynamic simulation calculation to standard output
     */
    void write(const char *signal);

  /**
 * Read events starting at the config cursor  and form a list of events
 * @param  cursor to Events in input.xml file
 * @return a list of events
 * Note : The configuration needs to be already set
 : setNetwork method
 *
 */
    std::vector<gridpack::dynamic_simulation::Event>
      getEvents(gridpack::utility::Configuration::CursorPtr cursor);

      /**
   * Read events set in the config file and form a list of events
   * @return a list of events
   * Note : The configuration needs to be already set
          : with the setNetwork method
   *
   */
    std::vector<gridpack::dynamic_simulation::Event>
      getEvents();

  /**
   * Set an event for the dynamic simulation
   * @param event info
   */
  void setEvent(gridpack::dynamic_simulation::Event event);

    /**
     * Read in generators that should be monitored during simulation
     * @param filename set filename from calling program instead of input
     *        deck
     * @param buses IDs of buses containing generators
     * @param tags generator IDs for watched generators
     * @param writeFile true if external file is to be written
     */
    void setGeneratorWatch();
    void setGeneratorWatch(const char *filename);
    void setGeneratorWatch(const char *filename,gridpack::utility::Configuration::CursorPtr cursor);

    void setGeneratorWatch(std::vector<int> &buses, std::vector<std::string> &tags,
        bool writeFile = true);
    void setGeneratorWatch(gridpack::utility::Configuration::CursorPtr cursor);


    /**
     * Read in loads that should be monitored during simulation
     */
    void setLoadWatch();

    /**
     * Redirect output from standard out
     * @param filename name of file to write results to
     */
    void open(const char *filename);
    void close();

    /**
     * Print string. This can be used to direct output to the file opened using
     * the open command
     * @param buf string to be printed
     */
    void print(const char *buf);

    /**
     * Check if system is secure
     */
    int isSecure();

    /**
     * Save watch series to an internal data vector
     * @param flag if true, save time series data
     */
    void saveTimeSeries(bool flag);

    /**
     * Return global map of timer series values
     * @return map of time series indices (local to global)
     */
    std::vector<int> getTimeSeriesMap();

    /**
     * Return vector of time series data for watched generators
     * @return vector of time series for generators on this processor
     */
    std::vector<std::vector<double> > getGeneratorTimeSeries();

    /**
     * Return a list of original bus IDs and tags for all monitored
     * generators
     * @param bus_ids list of original bus indices for monitored generators
     * @param gen_ids list of tags for monitored generators
     */
    void getListWatchedGenerators(std::vector<int> &bus_ids,
        std::vector<std::string> &gen_ids);

    /**
     * @return true if no frequency violations occured on monitored generators
     */
    bool frequencyOK();

    /**
     * Scale generator real power. If zone less than 1 then scale all
     * generators in the area.
     * @param scale factor to scale real power generation
     * @param area index of area for scaling generation
     * @param zone index of zone for scaling generation
     */
    void scaleGeneratorRealPower(double scale, int area, int zone);

    /**
     * Scale load power. If zone less than 1 then scale all
     * loads in the area.
     * @param scale factor to scale load real power
     * @param area index of area for scaling load
     * @param zone index of zone for scaling load
     */
    void scaleLoadPower(double scale, int area, int zone);

    /**
     * Return the total real power load for all loads in the zone. If zone
     * less than 1, then return the total load for the area
     * @param area index of area
     * @param zone index of zone
     * @return total load
     */
    double getTotalLoadRealPower(int area, int zone);

    /**
     * Return the current real power generation and the maximum and minimum total
     * power generation for all generators in the zone. If zone is less than 1
     * then return values for all generators in the area
     * @param area index of area
     * @param zone index of zone
     * @param total total real power generation
     * @param pmin minimum allowable real power generation
     * @param pmax maximum available real power generation
     */
    void getGeneratorMargins(int area, int zone, double *total, double *pmin,
        double *pmax);

    /**
     * Reset power of loads and generators to original values
     */
    void resetPower();

    /**
     * Write real time path rating diagnostics
     * @param src_area generation area
     * @param src_zone generation zone
     * @param load_area load area
     * @param load_zone load zone
     * @param gen_scale scale factor for generation
     * @param load_scale scale factor for loads
     * @param file name of file containing diagnostics
     */
    void writeRTPRDiagnostics(int src_area, int src_zone, int load_area,
        int load_zone, double gen_scale, double load_scale, const char *file);

    /**
     * Get a list of buses that had frequency violations
     * @return a list of buses that had frequency failures
     */
    std::vector<int> getFrequencyFailures();

    /**
     * Set parameters for monitoring frequency
     * @param flag true if frequency monitoring is turned on
     * @param maxFreq maximum allowable frequency deviation
     */
    void setFrequencyMonitoring(bool flag, double maxFreq);

    /**
     * Get observations and store them internally
     * @param cursor configuration pointer to observation block
     */
    void setObservations(gridpack::utility::Configuration::CursorPtr cursor);

    /**
     * Get bus and generator IDs for all observations
     * @param genBuses host IDs for all observed generators
     * @param genIDs character identifiers for all observed generators
     * @param loadBuses host IDs for all observed dynamic loads
     * @param loadIDs character identifiers for all observed dynamic loads
     * @param busIDs bus IDs for all observed buses
     */
    void getObservationLists(std::vector<int> &genBuses,
        std::vector<std::string> &genIDs,  std::vector<int> &loadBuses,
        std::vector<std::string> &loadIDs, std::vector<int> &busIDs);
		
	/**
	* Get bus and generator IDs for all observations including bus frequency ob
	* @param genBuses host IDs for all observed generators
	* @param genIDs character identifiers for all observed generators
	* @param loadBuses host IDs for all observed dynamic loads
	* @param loadIDs character identifiers for all observed dynamic loads
	* @param busIDs bus IDs for all observed buses
	* @param busfreqIDs bus IDs for all observed buses for bus frequency
	*/
	void getObservationLists_withBusFreq(
		std::vector<int> &genBuses, std::vector<std::string> &genIDs,
		std::vector<int> &loadBuses, std::vector<std::string> &loadIDs,
		std::vector<int> &busIDs, std::vector<int> &busfreqIDs);

    /**
     * Get current values of observations
     * @param vMag voltage magnitude for observed buses
     * @param vAng voltage angle for observed buses
     * @param rSpd rotor speed on observed generators
     * @param rAng rotor angle on observed generators
	 * @param genP real power on observed generators
     * @param genQ reactive power on observed generators
     * @param fOnline fraction of load shed
     */
    void getObservations(std::vector<double> &vMag, std::vector<double> &vAng,
        std::vector<double> &rSpd, std::vector<double> &rAng,
		std::vector<double> &genP, std::vector<double> &genQ,
        std::vector<double> &fOnline);
		
	/**
	* Get current values of observations including bus frequency ob
	* @param vMag voltage magnitude for observed buses
	* @param vAng voltage angle for observed buses
	* @param rSpd rotor speed on observed generators
	* @param rAng rotor angle on observed generators
	* @param genP real power on observed generators
	* @param genQ reactive power on observed generators
	* @param fOnline fraction of load shed
	* @param busfreq frequency of the buses in ob list
	*/
	void getObservations_withBusFreq(
		std::vector<double> &vMag, std::vector<double> &vAng,
		std::vector<double> &rSpd, std::vector<double> &rAng,
		std::vector<double> &genP, std::vector<double> &genQ,
		std::vector<double> &fOnline, std::vector<double> &busfreq);

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
    void getZoneLoads(std::vector<double> &load_p, std::vector<double> &load_q,
        std::vector<int> &zone_id) const;

    /**
     * Return total active and reactive generator power for each zone
     * @param generator_p active generator power for all zones
     * @param generator_q reactive generator power for all zones
     * @param zone_id label for each zone
     */
    void getZoneGeneratorPower(std::vector<double> &generator_p,
        std::vector<double> &generator_q, std::vector<int> &zone_id) const;

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

    /*
      Get the time-step
    */
    double getTimeStep();

    /*
      Set the time-step
    */
    void setTimeStep(double time_step);

    /*
      Set simulation end time
    */
    void setFinalTime(double final_time);

    /*
      Get simulation end time
    */
    double getFinalTime();

    /*
      Get current time
    */
    double getCurrentTime();

  /**
   * Transfer data from power flow to dynamic simulation
   * @param pf_network power flow network
   * @param ds_network dynamic simulation network
   */
  void transferPFtoDS(
    boost::shared_ptr<gridpack::powerflow::PFNetwork>
    pf_network,
    boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork>
    ds_network);


  // Network analytics queries

    /// Network query: Get the number of buses
  int totalBuses(void) const;

  /// Network query: Get the number of branches
  int totalBranches(void) const;

  /// Network query: Get the branch indexes connected to a bus
  std::vector<int> getConnectedBranches(int oidx) const;

  /// Network query: Get the bus indexes connected to a branch
  void getBranchEndpoints(const int& idx, int *fbus, int *tbus) const;

  /// Network query: Get the number of generators
  int numGenerators(void) const;

  /// Network query: Get the number of generators on a specific bus
  int numGenerators(const int& bus_idx) const;

  /// Network query: Get the number of loads
  int numLoads(void) const;

  /// Network query: Get the number of loads
  int numLoads(const int& bus_idx) const;

  /// Network query: Get the number of lines
  int numLines(void) const;

  /// Network query: Get the number of lines in a specific branch
  int numLines(const int& branch_idx) const;

  /// Network query: Get the number of storage units
  int numStorage(void) const;

  /// Network query: Get the number of storage units on a specific bus
  int numStorage(const int& bus_idx) const;

  /// Network query: Get a value from the bus' data collection
  template <typename T>
  bool
  getBusInfo(const int& bus_idx, const std::string& field,
             T& value, const int& dev_idx = -1)
  {
    bool ok(false);
    if (p_analytics) {
      ok = p_analytics->getBusInfo(bus_idx, field, value, dev_idx);
    }
    return ok;
  }
  

  
  private:

  double p_current_time; /* Current time */
  double p_sim_time;    // Simulation time
  double p_time_step;    /* Time-step */

  /**
     setLineStatus - Sets the line status and updates the associated
     branch and bus objects. 

     @param: from_idx - from bus number
     @param: to_idx - to bus number
     @param: ckt_id - circuit id
     @param: status - new line status

     Note: This method is called by handleEvents method to
           update the branch status and update the bus/branch
	   objects. It sets up values in the bus and branch objects
	   so that incrementMatrix method called on the network Ybus
	   uses these values to remove the branch contributions from
	   the Y-bus matrix
  **/
  void setLineStatus(int from_idx, int to_idx, std::string ckt_id, int status);

  /**
     setGenStatus - Sets the generator status and updates the associated
     bus objects. 

     @param: bus_idx - bus number
     @param: gen_id - generator id
     @param: status - new generator status

     Note: This method is called by handleEvents method to
           update the generator status and update the bus
	   object. It sets up values in the bus objects
	   so that incrementMatrix method called on the network Ybus
	   uses these values to remove the generator contributions from
	   the Y-bus matrix
  **/
  void setGenStatus(int bus_idx, std::string gen_id, int status);

  /**
     Handle any events
  **/
  void handleEvents();
  /**
     run one step of dynamics simulation
  **/
  void runonestep();

  /*
    Update Norton current injected in the network
    predcorrflag = 0 => Predictor stage
    predcorrflag = 1 => Corrector stage
  */
  void getCurrent(int predcorrflag);

  /*
    Solve Network equations
    predcorrflag = 0 => Predictor stage
    predcorrflag = 1 => Corrector stage
    returns true if netwok converged
  */
  bool solveNetwork(int predcorrflag);


  /**
   * Utility function to convert faults that are in event list into
   * internal data structure that can be used by code
   */
  void setFaultEvents();
  
  /**
   * Open file (specified in input deck) to write generator results to.
   * Rotor angle and speeds from generators specified in input deck will be
     * written to this file at specified time intervals
     */
    void openGeneratorWatchFile();

    /**
     * Close file contain generator watch results
     */
    void closeGeneratorWatchFile();

    /**
     * Open file (specified in input deck) to write load results to.
     * Data from loads specified in input deck will be
     * written to this file at specified time intervals
     */
    void openLoadWatchFile();

    /**
     * Close file contain generator watch results
     */
    void closeLoadWatchFile();

    /**
     * Save time series data for watched generators
     */
    void saveTimeStep();

    /**
     * Check to see if frequency variations on monitored generators are okay
     * @param start time at which to start monitoring
     * @param time current value of time
     * @return true if all watched generators are within acceptable bounds
     */
    bool checkFrequency(double start, double time);

    /**
     * Check to see if frequency variations on monitored generators are okay
     * @param limit maximum upper limit on frequency deviation
     * @return true if all watched generators are within acceptable bounds
     */
    bool checkFrequency(double limit);

    /**
     * Get a list of unique zones in the system
     * @param zones a complete list of all zones in the network
     */
    void getZoneList(std::vector<int> &zones) const;

    /**
     * Template function for modifying generator parameters in data collection
     * for specified bus
     * @param bus_id bus ID
     * @param gen_id two character token specifying generator on bus
     * @param gen_par string representing dictionary name of data element
     *                to be modified
     * @param value new value of parameter
     * @return return false if parameter is not found
     */
    template <typename T>
    bool p_modifyDataCollectionGenParam(
        int bus_id, std::string gen_id, std::string genParam, T value)
    {
      std::vector<int> indices = p_network->getLocalBusIndices(bus_id);
      if (indices.size() > 0) {
        int i;
        bool ret = false;
        for (i=0; i<indices.size(); i++) {
          boost::shared_ptr<gridpack::component::DataCollection> data =
            p_network->getBusData(indices[i]);
          int ngen;
          if (data->getValue(GENERATOR_NUMBER, &ngen)) {
            int igen = -1;
            T tval;
            int j;
            std::string gID;
            for (j = 0; j<ngen; j++) {
              // Get index of generator
              if (data->getValue(GENERATOR_ID,&gID,j)) {
                if (gID == gen_id) {
                  igen = j;
                  break;
                }
              }
            }
            if (igen >= 0) {
              if (data->setValue(genParam.c_str(),value,igen)) {
                ret = true;
              }
            }
          }
        }
        return ret;
      }
      return false;
    }

    /**
     * Template function for modifying load parameters in data collection for
     * specified bus
     * @param bus_id bus ID
     * @param load_id two character token specifying load on bus
     * @param load_par string representing dictionary name of data element
     *                to be modified
     * @param value new value of parameter
     * @return return false if parameter is not found
     */
    template <typename T>
    bool p_modifyDataCollectionLoadParam(int bus_id, std::string load_id,
          std::string loadParam, T value)
    {
      std::vector<int> indices = p_network->getLocalBusIndices(bus_id);
      if (indices.size() > 0) {
        int i;
        bool ret = false;
        for (i=0; i<indices.size(); i++) {
          boost::shared_ptr<gridpack::component::DataCollection> data =
            p_network->getBusData(indices[i]);
          int nload;
          if (data->getValue(LOAD_NUMBER, &nload)) {
            int iload = -1;
            T tval;
            int j;
            std::string lID;
            for (j = 0; j<nload; j++) {
              if (data->getValue(LOAD_ID,&lID,j)) {
                if (lID == load_id) {
                  iload = j;
                  break;
                }
              }
            }
            if (iload >= 0) {
              if (data->setValue(loadParam.c_str(),value,iload)) {
                ret = true;
              }
            }
          }
        }
        return ret;
      }
      return false;
    }

    /**
     * Template function for modifying parameters in data collection for
     * specified bus
     * @param bus_id bus ID
     * @param bus_par string representing dictionary name of data element
     *                to be modified
     * @param value new value of parameter
     * @return return false if parameter is not found
     */
    template <typename T>
    bool p_modifyDataCollectionBusParam(int bus_id, std::string busParam, 
        T value)
    {
      std::vector<int> indices = p_network->getLocalBusIndices(bus_id);
      if (indices.size() > 0) {
        int i;
        bool ret = false;
        for (i=0; i<indices.size(); i++) {
          boost::shared_ptr<gridpack::component::DataCollection> data =
            p_network->getBusData(indices[i]);
          T tval;
          if (data->setValue(busParam.c_str(),value)) {
            ret = true;
          }
        }
        return ret;
      }
      return false;
    }

    std::vector<gridpack::dynamic_simulation::Event> p_events;

    // pointer to network
    boost::shared_ptr<DSFullNetwork> p_network;

    // communicator for network
    gridpack::parallel::Communicator p_comm;

    // pointer to factory
    boost::shared_ptr<DSFullFactory> p_factory;
	
	// whether iteratively solve the network interface 
	bool p_biterative_solve_network;
	double ITER_TOL;  // iteratively solve the network interface tolerance, defined in xml file
	int MAX_ITR_NO;   // iteratively solve the network interface max iteration number, defined in xml file
	
	// whether print out debug information for iteratively solve the network interface 
	bool p_iterative_network_debug;
	
	//for the generator observations, output the generator power based on system base or generator base
	bool p_generator_observationpower_systembase;


    // Current step count?
    int p_S_Steps;

    // pointer to bus IO module
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<DSFullNetwork> >
      p_busIO;

    // pointer to branch IO module
    boost::shared_ptr<gridpack::serial_io::SerialBranchIO<DSFullNetwork> >
      p_branchIO;

    // pointer to configuration module
    gridpack::utility::Configuration *p_config;

    // file name for generator watch file
    std::string p_gen_watch_file;

    // suppress watch file output but still monitor watched generators
    bool p_suppress_watch_files;

    // flag indicating whether or not to use application supplied generator
    // watch file name or name from input deck
    bool p_internal_watch_file_name;

    // flag indicating whether or not generators to be monitored have
    // already been read in
    bool p_generators_read_in;

    // Flag indicating that generators are to be monitored
    bool p_generatorWatch;

    // Frequency to write out generator watch results
    int p_generatorWatchFrequency;

    // local bus indices of generators that are being monitored
    std::vector<int> p_gen_buses;

    // IDs of generators that are being monitored
    std::vector<std::string> p_gen_ids;

    // Global indices of buses that own monitored generators
    std::vector<int> p_watch_bus_ids;

    // Tags of generators that are being monitors
    std::vector<std::string> p_watch_gen_ids;
	
	//branches need to be tripped at a specific dynamic simulation step
	std::vector<gridpack::dynamic_simulation::DSFullBranch*> p_vbranches_need_to_trip;
	
	//flag indicating whether there is/are branches need to be tripped at a specific dynamic simulation step
	bool bapplyLineTripAction;
	
	//flag indicating whether there is/are bus load P or Q change at a specific dynamic simulation step
	bool bapplyLoadChangeP;
	bool bapplyLoadChangeQ;
	
	std::vector<gridpack::dynamic_simulation::DSFullBus*> p_vbus_need_to_changeP;
	std::vector<gridpack::dynamic_simulation::DSFullBus*>p_vbus_need_to_changeQ;


    // Monitor generators for instability
    bool p_monitorGenerators;
    double p_maximumFrequency;

    // Frequency deviations for simulation are okay
    bool p_frequencyOK;

    // pointer to bus IO module that is used for generator results
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<DSFullNetwork> >
      p_generatorIO;

    // Flag indicating that loads are to be monitored
    bool p_loadWatch;

    // Frequency to write out load watch results
    int p_loadWatchFrequency;

    // bus indices of loads that are being monitored
    std::vector<int> p_load_buses;

    // pointer to bus IO module that is used for load results
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<DSFullNetwork> >
      p_loadIO;

   // Keep track of whether or not systsem is secure
   int p_insecureAt;

   // Global list of all generators that are being watched
   std::map<std::pair<int,std::string>, int> p_watch_list;

   // Flag to save time series
   bool p_save_time_series;

   // Vector of times series from watched generators
   std::vector<std::vector<double> > p_time_series;

   // Record bus ID where frequency violation occured
   std::vector<int> p_violations;

   // Report observations even if the corresponding network elements do not
   // exist
   bool p_report_dummy_obs;

   // Observation data structures
   std::vector<int> p_obs_genBus;
   std::vector<std::string> p_obs_genIDs;
   std::vector<int> p_obs_loadBus;
   std::vector<std::string> p_obs_loadIDs;
   std::vector<int> p_obs_vBus;
   std::vector<int> p_obs_vBusfreq;

   std::vector<int> p_obs_lGenIdx;
   std::vector<int> p_obs_GenIdx;
   std::vector<int> p_obs_lGenBus;
   std::vector<std::string> p_obs_lGenIDs;
   std::vector<int> p_obs_gActive;
   std::vector<int> p_obs_gUse;

   std::vector<int> p_obs_lLoadIdx;
   std::vector<int> p_obs_LoadIdx;
   std::vector<int> p_obs_lLoadBus;
   std::vector<std::string> p_obs_lLoadIDs;
   std::vector<int> p_obs_lActive;
   std::vector<int> p_obs_lUse;

   std::vector<int> p_obs_lVIdx;
   std::vector<int> p_obs_VIdx;
   std::vector<int> p_obs_lVBus;
   std::vector<int> p_obs_vActive;
   std::vector<int> p_obs_vUse;
   
   std::vector<int> p_obs_lVIdxfreq;
   std::vector<int> p_obs_VIdxfreq;
   std::vector<int> p_obs_lVBusfreq;
   std::vector<int> p_obs_vActivefreq;
   std::vector<int> p_obs_vUsefreq;

   boost::shared_ptr<gridpack::parallel::GlobalVector<double> > p_obs_vMag;
   boost::shared_ptr<gridpack::parallel::GlobalVector<double> > p_obs_vAng;
   boost::shared_ptr<gridpack::parallel::GlobalVector<double> > p_obs_vBusfreqVal;
   boost::shared_ptr<gridpack::parallel::GlobalVector<double> > p_obs_rSpd;
   boost::shared_ptr<gridpack::parallel::GlobalVector<double> > p_obs_rAng;
   boost::shared_ptr<gridpack::parallel::GlobalVector<double> > p_obs_genP;
   boost::shared_ptr<gridpack::parallel::GlobalVector<double> > p_obs_genQ;
   boost::shared_ptr<gridpack::parallel::GlobalVector<double> > p_obs_fOnline;
   
   // below are all variables originally defined the solve function, now define them as class private members
   boost::shared_ptr < gridpack::mapper::FullMatrixMap<DSFullNetwork> > ybusMap_sptr;  
   boost::shared_ptr<gridpack::math::Matrix> orgYbus;  
   boost::shared_ptr<gridpack::math::Matrix> ybusyl;  
   boost::shared_ptr<gridpack::math::Matrix> ybuspg;
   boost::shared_ptr<gridpack::math::Matrix> ybus_jxd;
   boost::shared_ptr<gridpack::math::Matrix> ybus;
   boost::shared_ptr<gridpack::math::Matrix> ybus_fy;
   boost::shared_ptr<gridpack::math::Matrix> ybus_posfy;
   
   boost::shared_ptr < gridpack::mapper::BusVectorMap<DSFullNetwork> > ngenMap_sptr; 
   boost::shared_ptr<gridpack::math::Vector> volt;
   
   boost::shared_ptr < gridpack::mapper::BusVectorMap<DSFullNetwork> > nbusMap_sptr;
   boost::shared_ptr<gridpack::math::Vector> INorton_full;
   boost::shared_ptr<gridpack::math::Vector> INorton_full_chk;
   boost::shared_ptr<gridpack::math::Vector> volt_full;
   double max_INorton_full;
   
   boost::shared_ptr<gridpack::math::LinearSolver> solver_sptr;
   boost::shared_ptr<gridpack::math::LinearSolver> solver_fy_sptr;
   boost::shared_ptr<gridpack::math::LinearSolver> solver_posfy_sptr;

   // analytics module
   boost::shared_ptr<gridpack::analysis::NetworkAnalytics<DSFullNetwork> >
     p_analytics;
   
   int simu_total_steps;
   int S_Steps;
   int last_S_Steps;
   int steps3, steps2, steps1;
   double h_sol1, h_sol2;
   int flagP, flagC;
   int Simu_Current_Step;
   bool p_bDynSimuDone;
   
   const double sysFreq = 60.0;
   double pi = 4.0*atan(1.0);
   const double basrad = 2.0 * pi * sysFreq;
   //gridpack::ComplexType jay(0.0, 1.0);
   
   //timer for record the excution time for each main modules
   //gridpack::utility::CoarseTimer *timer;
   int t_solve;
   int t_execute_steps;
   int t_presolve;
   int t_misc;
   int t_mode;
   int t_ybus;
   int t_init;
   int t_mIf;
   int t_psolve;
   int t_vmap;
   int t_volt;
   int t_predictor;
   int t_csolve;
   int t_corrector;
   int t_secure;
   int t_cmIf;
   
};

} // dynamic simulation
} // gridpack
#endif
