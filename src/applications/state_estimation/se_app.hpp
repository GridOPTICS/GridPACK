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
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#ifndef _se_app_h_
#define _se_app_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/applications/state_estimation/se_factory.hpp"

namespace gridpack {
namespace state_estimation {

    // Calling program for state estimation application

class SEApp
{
  public:
    /**
     * Basic constructor
     */
    SEApp(void);

    /**
     * Basic destructor
     */
    ~SEApp(void);

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);

    private:
};

} // state estimation
} // gridpack
#endif
