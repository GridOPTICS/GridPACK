/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsf_app_module.hpp
 * @author Shuangshuang Jin
 * @date   Feb 04, 2015
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#ifndef _dsf_app_module_h_
#define _dsf_app_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "dsf_factory_module.hpp"


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
     */
    void readGenerators(void);

    /**
     * Set up exchange buffers and other internal parameters and initialize
     * network components using data from data collection
     */
    void initialize();

    /**
     * Execute the time integration portion of the application
     */
    void solve(gridpack::dynamic_simulation::DSFullBranch::Event fault);

    /**
     * Write out final results of dynamic simulation calculation to standard output
     */
    void write(const char *signal);

    /**
     * Read faults from external file and form a list of faults
     * @param cursor pointer to open file contain fault or faults
     * @return a list of fault events
     */
    std::vector<gridpack::dynamic_simulation::DSFullBranch::Event>
      getFaults(gridpack::utility::Configuration::CursorPtr cursor);

    /**
     * Read in generators that should be monitored during simulation
     */
    void setGeneratorWatch();

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

    std::vector<gridpack::dynamic_simulation::DSFullBranch::Event> p_faults;

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


    // Flag indicating that generators are to be monitored
    bool p_generatorWatch;

    // Frequency to write out generator watch results
    int p_watchFrequency;

    // bus indices of generators that are being monitored
    std::vector<int> p_gen_buses;

    // tags of generators that are being monitored
    std::vector<std::string> p_tags;

    // pointer to bus IO module that is used for generator results
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<DSFullNetwork> >
      p_generatorIO;

   // Keep track of whether or not systsem is secure
   int p_insecureAt;

};

} // dynamic simulation
} // gridpack
#endif
