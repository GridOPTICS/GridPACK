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
#include "gridpack/include/gridpack.hpp"
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

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);

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

    private:
      //boost::shared_ptr<KalmanNetwork> p_network;
      //boost::shared_ptr<KalmanFactory> p_factory;
    gridpack::parallel::Communicator p_comm;
    int p_nsteps;
    double p_delta_t;
};

} // state estimation
} // gridpack
#endif
