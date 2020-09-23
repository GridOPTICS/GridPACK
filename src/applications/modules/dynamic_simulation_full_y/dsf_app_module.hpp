/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsf_app_module.hpp
 * @author Shuangshuang Jin
 * @date   2020-02-12 12:17:37 d3g096
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
#include "dsf_factory.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/math/math.hpp"
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
     * Set up exchange buffers and other internal parameters and initialize
     * network components using data from data collection
     */
    void initialize();

    /**
     * Reinitialize calculation from data collections
     */
    void reload();

    /**
     * Execute the time integration portion of the application
     */
    void solve(gridpack::dynamic_simulation::Event fault);
	
	/**
     * initialization before the time step integration starts 
     */
    void solvePreInitialize(gridpack::dynamic_simulation::Event fault);
	
	/**
     * Execute only one simulation time step 
     */
    void executeOneSimuStep();
	
	/**
     * execute load shedding	 
     */
    void applyLoadShedding(int bus_number, std::string loadid, double percentage);
	
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
     * Read faults from external file and form a list of faults
     * @param cursor pointer to open file contain fault or faults
     * @return a list of fault events
     */
    std::vector<gridpack::dynamic_simulation::Event>
      getFaults(gridpack::utility::Configuration::CursorPtr cursor);

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
    void setGeneratorWatch(std::vector<int> &buses, std::vector<std::string> &tags,
        bool writeFile = true);

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
     * Get current values of observations
     * @param vMag voltage magnitude for observed buses
     * @param vAng voltage angle for observed buses
     * @param rSpd rotor speed on observed generators
     * @param rAng rotor angle on observed generators
     * @param fOnline fraction of load shed
     */
    void getObservations(std::vector<double> &vMag, std::vector<double> &vAng,
        std::vector<double> &rSpd, std::vector<double> &rAng,
        std::vector<double> &fOnline);

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

  private:
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

    std::vector<gridpack::dynamic_simulation::Event> p_faults;

    // pointer to network
    boost::shared_ptr<DSFullNetwork> p_network;

    // communicator for network
    gridpack::parallel::Communicator p_comm;

    // pointer to factory
    boost::shared_ptr<DSFullFactory> p_factory;

    // Simulation time
    double p_sim_time;

    // Time step
    double p_time_step;

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

   // flag indicating that a line is being tripped somewhere in the system.
   // Needed to help manage matrix modications, which is a collective operation
   bool p_lineTrip;
	

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

   // Observation data structures
   std::vector<int> p_obs_genBus;
   std::vector<std::string> p_obs_genIDs;
   std::vector<int> p_obs_loadBus;
   std::vector<std::string> p_obs_loadIDs;
   std::vector<int> p_obs_vBus;

   std::vector<int> p_obs_lGenIdx;
   std::vector<int> p_obs_GenIdx;
   std::vector<int> p_obs_lGenBus;
   std::vector<std::string> p_obs_lGenIDs;
   std::vector<int> p_obs_gActive;

   std::vector<int> p_obs_lLoadIdx;
   std::vector<int> p_obs_LoadIdx;
   std::vector<int> p_obs_lLoadBus;
   std::vector<std::string> p_obs_lLoadIDs;
   std::vector<int> p_obs_lActive;

   std::vector<int> p_obs_lVIdx;
   std::vector<int> p_obs_VIdx;
   std::vector<int> p_obs_lVBus;
   std::vector<int> p_obs_vActive;

   boost::shared_ptr<gridpack::parallel::GlobalVector<double> > p_obs_vMag;
   boost::shared_ptr<gridpack::parallel::GlobalVector<double> > p_obs_vAng;
   boost::shared_ptr<gridpack::parallel::GlobalVector<double> > p_obs_rSpd;
   boost::shared_ptr<gridpack::parallel::GlobalVector<double> > p_obs_rAng;
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
