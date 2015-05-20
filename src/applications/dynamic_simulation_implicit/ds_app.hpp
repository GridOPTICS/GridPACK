/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_app.hpp
 * @author Shrirang Abhyankar
 * @date   2015-01-23 d3g293
 * 
 * @brief  
 * Application class definition
 * 
 */
// -------------------------------------------------------------

#ifndef _ds_app_h_
#define _ds_app_h_

#include "boost/smart_ptr/shared_ptr.hpp"

namespace gridpack {
namespace dsimplicit {

// Calling program for implicit dynamics simulation application

class DSApp
{
  public:
    /**
     * Basic constructor
     */
    DSApp(void);

    /**
     * Basic destructor
     */
    ~DSApp(void);

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);

  private:
};

} // dsimplicit
} // gridpack
#endif
