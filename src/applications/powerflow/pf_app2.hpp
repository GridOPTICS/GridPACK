/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_app2.hpp
 * @author William A. Perkins
 * @date   2014-01-28 11:30:13 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _pf_app_h_
#define _pf_app_h_

#include "boost/smart_ptr/scoped_ptr.hpp"
#include "pf_factory.hpp"
#include "gridpack/timer/coarse_timer.hpp"

namespace gridpack {
namespace powerflow {

// Calling program for powerflow application

class PFApp2
{
public:
  /**
   * Basic constructor
   */
  PFApp2(void);
  
  /**
   * Basic destructor
   */
  ~PFApp2(void);
  
  /**
   * Execute application
   */
  void execute(int argc, char** argv);
  
private:
};

} // powerflow
} // gridpack
#endif
