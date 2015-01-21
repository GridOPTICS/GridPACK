/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_app_module.hpp
 * @author Shuangshuang Jin
 * @date   September 19, 2013
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#ifndef _ds_app_module_h_
#define _ds_app_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/include/gridpack.hpp"
#include "ds_factory_module.hpp"

namespace gridpack {
namespace dynamic_simulation {

    // Calling program for dynamic simulation application

class DSAppModule
{
  public:
    /**
     * Basic constructor
     */
    DSAppModule(void);

    /**
     * Basic destructor
     */
    ~DSAppModule(void);

    /**
     * Read in and partition the powerflow network. The input file is read
     * directly from the Dynamic_simulation block in the configuration file so no
     * external file names or parameters need to be passed to this routine
     * @param network pointer to a DSNetwork object. This should not have any
     * buses or branches defined on it.
     * @param config pointer to open configuration file
     */
    void readNetwork(boost::shared_ptr<DSNetwork> &network,
        gridpack::utility::Configuration *config);

    /**
     * Assume that DSNetwork already exists and just cache an internal pointer
     * to it. This routine does not call the partition function. Also read in
     * simulation parameters from configuration file
     * @param network pointer to a complete DSNetwork object.
     * @param config pointer to open configuration file
     */
    void setNetwork(boost::shared_ptr<DSNetwork> &network,
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
    void solve();

    /**
     * Write out final results of dynamic simulation calculation to standard output
     */
    void write();

  private:

    /**
     * Utility function to convert faults that are in event list into
     * internal data structure that can be used by code
     */
    void setFaultEvents();

    std::vector<gridpack::dynamic_simulation::DSBranch::Event> p_faults;

    // pointer to network
    boost::shared_ptr<DSNetwork> p_network;

    // communicator for network
    gridpack::parallel::Communicator p_comm;

    // pointer to factory
    boost::shared_ptr<DSFactoryModule> p_factory;

    // Simulation time
    double p_sim_time;

    // Time step
    double p_time_step;

    // Current step count?
    int p_S_Steps;

    // pointer to bus IO module
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<DSNetwork> >
      p_busIO;

    // pointer to branch IO module
    boost::shared_ptr<gridpack::serial_io::SerialBranchIO<DSNetwork> >
      p_branchIO;

    // pointer to configuration module
    gridpack::utility::Configuration *p_config;
};

} // dynamic simulation
} // gridpack
#endif
