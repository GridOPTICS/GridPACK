/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   DBIntClass.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 10, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _dbintclass_h_
#define _dbintclass_h_

#include "boost/smart_ptr/shared_ptr.hpp"

namespace gridpack {
namespace dynamic_simulation {
class DBIntClass 
{
  public:
    /**
     * Basic constructor
     */
    DBIntClass();

    /**
     * Basic destructor
     */
    ~DBIntClass();

    /**
     * @param theDb1 
     * @param theEps 
     * @param theInit
     */
    void Initialize(double theDb1, double theEps, double theInit);

    /**
     * @param theInput
     */ 
    double Output(double anInput);

  private:

    double Db1;
    double Eps;
    int State;
    double InitialValue;
   
};
}  // dynamic_simulation
}  // gridpack
#endif
