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

#ifndef _vs_app_module_h_
#define _vs_app_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "vs_factory_module.hpp"

namespace gridpack {
namespace voltage_stability {

// Calling program for state estimation application
// Calling program for state estimation application
class VSAppModule
{
  public:
    /**
     * Basic constructor
     */
    //VSAppModule(gridpack::parallel::Communicator comm);
    VSAppModule(void);

    /**
     * Basic destructor
     */
    ~VSAppModule(void);

    void transferPFtoVS(boost::shared_ptr<gridpack::powerflow::PFNetwork> pf_network,boost::shared_ptr<VSNetwork> VS_network);
    void solvepowerflow(void);
    
    void write_FVSI();
  
    private:

    // pointer to network
    //boost::shared_ptr<VSNetwork> p_network;
    
    // Network pointer
    boost::shared_ptr<VSNetwork> vs_network;
    // Power flow application pointer
    gridpack::powerflow::PFAppModule *p_pfapp;
  
    // communicator for network
    gridpack::parallel::Communicator p_comm;

    // pointer to factory
    boost::shared_ptr<VSFactoryModule> p_factory;

    // pointer to bus IO module
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<VSNetwork> >
      p_busIO;

    // pointer to branch IO module
    boost::shared_ptr<gridpack::serial_io::SerialBranchIO<VSNetwork> >
      p_branchIO;

    // pointer to configuration module
    gridpack::utility::Configuration *p_config;
    gridpack::utility::Configuration::CursorPtr p_configcursor;
    char p_configfile[256];

    // maximum number of iterations
    int p_max_iteration;

    // convergence tolerance
    double p_tolerance;
};

} // state estimation
} // gridpack
#endif
