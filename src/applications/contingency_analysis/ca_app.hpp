/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ca_app.hpp
 * @author Yousu Chen 
 * @date   January 20, 2014
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#ifndef _ca_app_h_
#define _ca_app_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parallel/communicator.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/applications/contingency_analysis/ca_factory.hpp"
#include "gridpack/applications/contingency_analysis/ca_driver.hpp"

namespace gridpack {
namespace contingency_analysis {

// Calling program for contingency analysis application
class CAApp
{
  public:
    /**
     * Basic constructor
     * @param comm communicator that hosts contingency application
     */
    CAApp(gridpack::parallel::Communicator comm);

    /**
     * Basic destructor
     */
    ~CAApp(void);

    /**
     * Initialize application by reading in grid network, partioning it and
     * setting up buffers and indices
     * @param argc number of arguments
     * @param argv list of character strings
     */ 
    void init(int argc, char** argv);

    /**
     * Execute application for a particular contingency
     * @param contingency data structure describing the contingency
     */
    void execute(gridpack::contingency_analysis::Contingency contingency);

  private:
    boost::shared_ptr<CANetwork> p_network;
    boost::shared_ptr<CAFactory> p_factory;
};

} // contingency analysis 
} // gridpack
#endif
