/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   kds_app_module.hpp
 * @author Da Meng and Yousu Chen 
 * @date   1/06/2015
 *
 * @brief
 *
 * @ Modified by Xinya Li 6/30/2015
 */
// -------------------------------------------------------------

#ifndef _kalmer_app_h_
#define _kalmer_app_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "kds_factory_module.hpp"

namespace gridpack {
namespace kalman_filter {

// Calling program for state estimation application
// Calling program for state estimation application
class KalmanApp
{
  public:
    /**
     * Basic constructor
     */
    //KalmanApp(gridpack::parallel::Communicator comm);
    KalmanApp(void);

    /**
     * Basic destructor
     */
    ~KalmanApp(void);

    /**
     * Read in and partition the network. The input file is read
     * directly from the state_estimation block in the configuration file
     * so no external file names or parameters need to be passed to this
     * routine
     * @param network pointer to a KalmanNetwork object. This should not
     * have any buses or branches defined on it.
     * @param config pointer to open configuration file
     */
    void readNetwork(boost::shared_ptr<KalmanNetwork> &network,
        gridpack::utility::Configuration *config);

    /**
     * Assume that KalmanNetwork already exists and just cache an internal
     * pointer to it. This routine does not call the partition function.
     * Also read in simulation parameters from configuration file
     * @param network pointer to a complete KalmanNetwork object.
     * @param config pointer to open configuration file
     */
    void setNetwork(boost::shared_ptr<KalmanNetwork> &network,
        gridpack::utility::Configuration *config);

    /**
     * Set up exchange buffers and other internal parameters and initialize
     * network components using data from data collection
     */
    void initialize();

    /**
     * Perform the Kalman Filter simulation
     */
    void solve();

    private:

    /**
     * Utility function to convert faults that are in event list into
     * internal data structure that can be used by code
     * @param cursors list of cursors pointing to individual events in input
     * deck
     * @return list of event data structures
     */
    std::vector<gridpack::kalman_filter::KalmanBranch::Event>
      setFaultEvents(std::vector<gridpack::utility::Configuration::CursorPtr >
          cursors);


    /**
     * Get the time series measurements from a file
     * @param filename name of external file with time series data
     * @param series a vector of pointers to arrays of time series data
     * @param nsteps the number of timesteps in each of the time series
     * @param keys list of bus indices that own the time series data
     */
    void getTimeSeries(std::string filename, std::vector<double*> &series,
        double &delta_t, int &nsteps, std::vector<int> &keys);

    /**
     * Set up time series data for all buses
     * @param network pointer to KalmanNetwork object
     * @param cursor pointer to data in input deck
     */
    void setTimeData(boost::shared_ptr<KalmanNetwork> &network,
        gridpack::utility::Configuration::CursorPtr cursor);

    // Pointer to network
    boost::shared_ptr<KalmanNetwork> p_network;

    // Pointer to factory
    boost::shared_ptr<KalmanFactory> p_factory;

    // Communicator for network
    gridpack::parallel::Communicator p_comm;

    // Pointer to configuration module
    gridpack::utility::Configuration *p_config;

    // Number of simulation steps
    int p_nsteps;

    // Time step increment
    double p_delta_t;

    // Time step
    double p_time_step;

    // Simulation time
    double p_sim_time;

    // Known fault
    int p_KnownFault;

    // Time offset
    int p_TimeOffset;

    // Check equation flag
    int p_CheckEqn;

    // Fault list
    std::vector<gridpack::kalman_filter::KalmanBranch::Event> p_faults;

    // Matrix parameters
    double p_Rm1;
    double p_Rm1n;
    double p_N_inv;

    // Serial IO modules
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<
      gridpack::kalman_filter::KalmanNetwork> > p_busIO;
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<
      gridpack::kalman_filter::KalmanNetwork> > p_deltaIO;
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<
      gridpack::kalman_filter::KalmanNetwork> > p_omegaIO;
    boost::shared_ptr<gridpack::serial_io::SerialBranchIO<
      gridpack::kalman_filter::KalmanNetwork> > p_branchIO;
};

} // state estimation
} // gridpack
#endif
