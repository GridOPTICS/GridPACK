/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_app_module.hpp
 * @author Yousu Chen, Bruce Palmer
 * @date   1/23/2015
 * @Last modified 1/23/2015
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#ifndef _se_app_module_h_
#define _se_app_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "se_factory_module.hpp"

namespace gridpack {
namespace state_estimation {

// Calling program for state estimation application
// Calling program for state estimation application
class SEAppModule
{
  public:
    /**
     * Basic constructor
     */
    //SEAppModule(gridpack::parallel::Communicator comm);
    SEAppModule(void);

    /**
     * Basic destructor
     */
    ~SEAppModule(void);

    /**
     * Get list of measurements from external file
     * @param cursor pointer to contingencies in input deck
     * @return vector of measurements
     */
    std::vector<gridpack::state_estimation::Measurement> getMeasurements(
        gridpack::utility::Configuration::ChildCursors measurements);

    /**
     * Read in and partition the network. The input file is read
     * directly from the state_estimation block in the configuration file so no
     * external file names or parameters need to be passed to this routine
     * @param network pointer to a SENetwork object. This should not have any
     * buses or branches defined on it.
     * @param config pointer to open configuration file
     */
    void readNetwork(boost::shared_ptr<SENetwork> &network,
        gridpack::utility::Configuration *config);

    /**
     * Assume that SENetwork already exists and just cache an internal pointer
     * to it. This routine does not call the partition function. Also read in
     * simulation parameters from configuration file
     * @param network pointer to a complete SENetwork object.
     * @param config pointer to open configuration file
     */
    void setNetwork(boost::shared_ptr<SENetwork> &network,
        gridpack::utility::Configuration *config);

    /**
     * Read branch and bus measurements. These will come from a separate file.
     * The name of this file comes from the input configuration file.
     */
    void readMeasurements(void);

    /**
     * Set up exchange buffers and other internal parameters and initialize
     * network components using data from data collection
     */
    void initialize();

    /**
     * Solve the state estimation problem
     */
    void solve();

    /**
     * Write final results of state estimation calculation to standard output
     */
    void write();

    /**
     * Save results of state estimation calculation to data collection objects
     */
    void saveData();

    private:

    // pointer to network
    boost::shared_ptr<SENetwork> p_network;

    // communicator for network
    gridpack::parallel::Communicator p_comm;

    // pointer to factory
    boost::shared_ptr<SEFactoryModule> p_factory;

    // pointer to bus IO module
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<SENetwork> >
      p_busIO;

    // pointer to branch IO module
    boost::shared_ptr<gridpack::serial_io::SerialBranchIO<SENetwork> >
      p_branchIO;

    // pointer to configuration module
    gridpack::utility::Configuration *p_config;

    // maximum number of iterations
    int p_max_iteration;

    // convergence tolerance
    double p_tolerance;
};

} // state estimation
} // gridpack
#endif
