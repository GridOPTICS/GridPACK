/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   BackLashClass.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 10, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _backlashclass_h_
#define _backlashclass_h_

#include "boost/smart_ptr/shared_ptr.hpp"

namespace gridpack {
namespace dynamic_simulation {
class BackLashClass 
{
  public:
    /**
     * Basic constructor
     */
    BackLashClass();

    /**
     * Basic destructor
     */
    ~BackLashClass();

    /**
     * @param theDb2 
     * @param InitOutput 
     */
    void Initialize(double theDb2, double InitOutput);

    /**
     * @param theInput
     */ 
    double Output(double theInput);

  private:

    double Db2;
    double LastOutput;
   
};
}  // dynamic_simulation
}  // gridpack
#endif
