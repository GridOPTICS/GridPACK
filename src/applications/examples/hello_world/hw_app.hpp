/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hw_app.hpp
 * @author Bruce Palmer
 * @date   2013-08-08 10:18:32 d3g293
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _hw_app_h_
#define _hw_app_h_

#include "boost/smart_ptr/shared_ptr.hpp"

namespace gridpack {
namespace hello_world {

// Calling program for hello world application

class HWApp
{
  public:
    /**
     * Basic constructor
     */
    HWApp(void);

    /**
     * Basic destructor
     */
    ~HWApp(void);

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);

  private:
};

} // hello_world
} // gridpack
#endif
