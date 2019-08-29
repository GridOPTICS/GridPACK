/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   GainBlockClass.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 11, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _gainblockclass_h_
#define _gainblockclass_h_

#include "boost/smart_ptr/shared_ptr.hpp"

namespace gridpack {
namespace dynamic_simulation {
class GainBlockClass 
{
  public:
    /**
     * Basic constructor
     */
    GainBlockClass();

    /**
     * Basic destructor
     */
    ~GainBlockClass();

    /**
     * @param theX 
     */
    double XtoY(double theX);

    /**
     * @param theX 
     */
    double YtoX(double theY);

  private:

    double X[5];
    double Y[5];
    int Count;
    bool ExtrapolatePastEnd;
   
};
}  // dynamic_simulation
}  // gridpack
#endif
