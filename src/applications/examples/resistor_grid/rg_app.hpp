/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   rg_app.hpp
 * @author Bruce Palmer
 * @date   2014-01-28 d3g293
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _rg_app_h_
#define _rg_app_h_

#include "boost/smart_ptr/shared_ptr.hpp"

namespace gridpack {
namespace resistor_grid {

// Calling program for resistor application

class RGApp
{
  public:
    /**
     * Basic constructor
     */
    RGApp(void);

    /**
     * Basic destructor
     */
    ~RGApp(void);

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);

  private:
};

} // resistor_grid
} // gridpack
#endif
