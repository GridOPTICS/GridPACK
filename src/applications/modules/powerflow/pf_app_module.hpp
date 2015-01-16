/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_app_module.hpp
 * @author Bruce Palmer
 * @date   2014-01-28 11:30:49 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _pf_app_module_h_
#define _pf_app_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "pf_factory_module.hpp"

namespace gridpack {
namespace powerflow {

// Calling program for powerflow application

class PFAppModule
{
  public:
    /**
     * Basic constructor
     */
    PFAppModule(void);

    /**
     * Basic destructor
     */
    ~PFAppModule(void);

    /**
     * Read in and partition the powerflow network. The input file is read
     * directly from the Powerflow block in the configuration file so no
     * external file names or parameters need to be passed to this routine
     * @param network pointer to a PFNetwork object. This should not have any
     * buses or branches defined on it.
     * @param config point to open configuration file
     */
    void readNetwork(boost::shared_ptr<PFNetwork> &network,
                     gridpack::utility::Configuration *config);

    /**
     * Set up exchange buffers and other internal parameters and initialize
     * network components using data from data collection
     */
    void initialize();

    /**
     * Execute the iterative solve portion of the application
     */
    void solve();

    /**
     * Write out results of powerflow calculation to standard output
     */
    void write();

  private:

    // pointer to network
    boost::shared_ptr<PFNetwork> p_network;

    // communicator for network
    gridpack::parallel::Communicator p_comm;

    // pointer to factory
    boost::shared_ptr<PFFactoryModule> p_factory;

    // maximum number of iterations
    int p_max_iteration;

    // convergence tolerance
    double p_tolerance;

    // pointer to bus IO module
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<PFNetwork> > p_busIO;

    // pointer to branch IO module
    boost::shared_ptr<gridpack::serial_io::SerialBranchIO<PFNetwork> > p_branchIO;

    // pointer to configuration module
    gridpack::utility::Configuration *p_config;
};

} // powerflow
} // gridpack
#endif
