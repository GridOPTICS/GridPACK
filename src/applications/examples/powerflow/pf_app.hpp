/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_app.hpp
 * @author Bruce Palmer
 * @date   2014-01-28 11:30:49 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _pf_app_h_
#define _pf_app_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "pf_factory.hpp"

namespace gridpack {
namespace powerflow {

// Calling program for powerflow application

class PFApp
{
  public:
    /**
     * Basic constructor
     */
    PFApp(void);

    /**
     * Basic destructor
     */
    ~PFApp(void);

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);

  private:
};

} // powerflow
} // gridpack
#endif
