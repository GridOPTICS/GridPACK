/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hadrec_app_module.hpp
 * @author Bruce Palmer
 * @date   2014-01-28 11:30:49 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _hadrec_app_module_h_
#define _hadrec_app_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/configuration/configuration.hpp"

namespace gridpack {
namespace hadrec {

// Calling program for powerflow application

class HADRECAppModule
{
  public:
    /**
     * Basic constructor
     */
    HADRECAppModule(void);

    /**
     * Basic destructor
     */
    ~HADRECAppModule(void);

  private:

};

} // hadrec
} // gridpack
#endif
