/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   GainBlockClass.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   10/02/2019 
 * 
 * @brief  
 * 
 * 
 */

#ifndef _gainblockclass_h_
#define _gainblockclass_h_

#include "boost/smart_ptr/shared_ptr.hpp"

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
#endif
