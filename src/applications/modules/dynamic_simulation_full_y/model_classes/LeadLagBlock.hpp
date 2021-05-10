/*
 *     Copyright (c) 2021 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   LeadLagBlock.hpp
 * @author Renke Huang renke.huang@pnnl.gov
 * @Last modified:   May 7, 2021
 * 
 * @brief  
 *   
 *    In -> (1 + s*T2) / (1 + s*T1) ->Out
 * 
 */

#ifndef _leadlagblock_h_
#define _leadlagblock_h_

//#include "boost/smart_ptr/shared_ptr.hpp"

namespace gridpack {
namespace dynamic_simulation {
class LeadLagBlock 
{
  public:
    /**
     * Basic constructor
     */
    LeadLagBlock();

    /**
     * Basic destructor
     */
    ~LeadLagBlock();

    /**
     * @init, initialize the block with parameters and the states  x0, x1, as well as dx0, dx1 
     */
    double init(double dOut, double T1, double T2, double t_inc);

    double predictor(double In, double t_inc, bool flag);
	
	double corrector(double In, double t_inc, bool flag);

	// varaibles, make it public for easy access
    double T1, T2;
    double x0, x1, dx0, dx1;
   
};
}  // dynamic_simulation
}  // gridpack
#endif
