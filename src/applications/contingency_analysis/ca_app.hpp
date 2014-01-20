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
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/applications/contingency_analysis/ca_factory.hpp"

namespace gridpack {
namespace contingency_analysis {

// Calling program for contingency analysis application
class CAApp
{
  public:
    /**
     * Basic constructor
     */
    CAApp(void);

    /**
     * Basic destructor
     */
    ~CAApp(void);

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);

    private:
};

} // contingency analysis 
} // gridpack
#endif
