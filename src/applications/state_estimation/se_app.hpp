/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_app.hpp
 * @author Yousu Chen 
 * @date   2/24/2014 
 * @Last modified 8/5/2014
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#ifndef _se_app_h_
#define _se_app_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/include/gridpack.hpp"
#include "se_factory.hpp"

namespace gridpack {
namespace state_estimation {

// Calling program for state estimation application
// Calling program for state estimation application
class SEApp
{
  public:
    /**
     * Basic constructor
     */
    //SEApp(gridpack::parallel::Communicator comm);
    SEApp(void);

    /**
     * Basic destructor
     */
    ~SEApp(void);

    /**
     * Get list of measurements from external file
     * @param cursor pointer to contingencies in input deck
     * @return vector of measurements
     */
    std::vector<gridpack::state_estimation::Measurement> getMeasurements(
        gridpack::utility::Configuration::ChildCursors measurements);

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    //void execute(gridpack::state_estimation::Measurement measurement);
    void execute(int argc, char** argv);

    private:
      //boost::shared_ptr<SENetwork> p_network;
      //boost::shared_ptr<SEFactory> p_factory;
    gridpack::parallel::Communicator p_comm;
};

} // state estimation
} // gridpack
#endif
