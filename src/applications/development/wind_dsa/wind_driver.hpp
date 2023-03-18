/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wind_driver.hpp
 * @author Bruce Palmer
 * @date   February 10, 2014
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#ifndef _wind_driver_h_
#define _wind_driver_h_

#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "gridpack/applications/modules/dynamic_simulation_full_y/dsf_app_module.hpp"

namespace gridpack {
namespace contingency_analysis {

class QuantileAnalysis
{
  public:
  /**
   * Basic constructor
   * @param comm communicator used for analysis
   * @param nwatch number of generators being watched
   * @param nconf number of scenarios
   * @param nsteps number of timesteps being stored
   */
  QuantileAnalysis(gridpack::parallel::Communicator comm, int nwatch, int nconf, int nsteps);

  /**
   * Basic destructor
   */
  ~QuantileAnalysis();

  /**
   * Save data for a single time step for a single generator
   * @param cfg_idx scenario index for time series
   * @param gen_idx generator index for time series
   * @param vals vector of time series values for a generator
   */
  void saveData(int cfg_idx, int gen_idx, std::vector<double> &vals);

  /**
   * Save variable names
   * @param name vector of variable names
   */
  void saveVarNames(std::vector<std::string> &names);

  /**
   * Stream data in storage array
   */
  void writeData();

  /**
   * Calculate quantiles and write them to a file
   * @param quantiles values describing quantiles to be calculated.
   *                  These values should be between 0 and 1.
   * @param dt magnitude time step (in seconds)
   */
  void exportQuantiles(std::vector<double> quantiles, double dt);

  private:

    gridpack::parallel::Communicator p_comm;
    int p_nwatch;
    std::vector<std::string> p_var_names;
    int p_nconf;
    int p_nsteps;
    int p_GA;
};

typedef gridpack::dynamic_simulation::DSFullNetwork DSFullNetwork;

enum ContingencyType{Generator, Branch};

// Calling program for contingency analysis application
class WindDriver
{
  public:
    /**
     * Basic constructor
     */
    WindDriver(void);

    /**
     * Basic destructor
     */
    ~WindDriver(void);

    /**
     * Transfer data from power flow to dynamic simulation
     * @param pf_network power flow network
     * @param ds_network dynamic simulation network
     */
     void transferPFtoDS(boost::shared_ptr<gridpack::powerflow::PFNetwork>
        pf_network,
        boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork>
        ds_network);

    /**
     * Read faults from external file and form a list of faults
     * @param cursor pointer to open file contain fault or faults
     * @return a list of fault events
     */
    std::vector<gridpack::dynamic_simulation::Event>
      getEvents(gridpack::utility::Configuration::CursorPtr cursor);

    /**
     * Modify real power parameters in DS network based on values in
     * wind and load files
     * @param genIDs list of bus IDs containing generators to be modified
     * @param genTags list of tags for generators to be modified
     * @param windVals list of new real power values for generators
     * @param loadIDs list of bus IDs containing loads to be modified
     * @param loadTags list of tags for loads to be modified
     * @param loadVals list of new real power values for loads
     * @param pf_network pointer to PF network
     * @param ds_network pointer to DS network
     */
    void resetData(std::vector<int> &genIDs,
        std::vector<std::string> &genTags,
        std::vector<double> &windVals,
        std::vector<int> &loadIDs,
        std::vector<std::string> &loadTags,
        std::vector<double> &loadVals,
        boost::shared_ptr<gridpack::powerflow::PFNetwork> &pf_network,
        boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork> &ds_network);

    /**
     * Read list of branches to monitor in power flow calculation and store
     * contents in internal variables
     * @param cursor pointer to open file containing list of branchs
     */
    void getWatchLines(gridpack::utility::Configuration::CursorPtr cursor);

    /**
     * Get power flow values from watched branches
     * @param p real power values
     * @param q reactive power values
     * @param network pointer to PF network
     */
    void getWatchedBranches(std::vector<double> &p, std::vector<double> &q,
        boost::shared_ptr<gridpack::powerflow::PFNetwork> &network);

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);

    private:

    boost::shared_ptr<DSFullNetwork> p_network;

    /**
     * parameters for getting power flow data
     */
    bool p_watch_lines;
    std::vector<int> p_from;
    std::vector<int> p_to;
    std::vector<std::string> p_tags;
};

} // contingency analysis 
} // gridpack
#endif
