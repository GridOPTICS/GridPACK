/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   uc_app.hpp
 * @author 
 * @date   
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _uc_app_h_
#define _uc_app_h_

#include "uc_factory.hpp"

namespace gridpack {
namespace unit_commitment {

// Calling program for unit commitment application

class UCApp
//  : public gridpack::optimization::VariableVisitor 
{
  public:
    /**
     * Basic constructor
     */
    UCApp(void);

    /**
     * Basic destructor
     */
    ~UCApp(void);

    /**
     * Get time series data for loads and reserves and move them to individual
     * buses
     * @param filename name of file containing load and reserves time series
     * data
     */
    void getLoadsAndReserves(const char* filename);

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);


  private:

    boost::shared_ptr<UCNetwork> p_network;
};

} // unit_commitment
} // gridpack
#endif
