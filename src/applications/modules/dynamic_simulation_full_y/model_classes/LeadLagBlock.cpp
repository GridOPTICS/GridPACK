/*
 *     Copyright (c) 2021 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   LeadLagBlock.cpp
 * @author Renke Huang renke.huang@pnnl.gov
 * @Last modified:   May 7, 2021
 * 
 * @brief  
 * 
 *    In -> (1 + s*T2) / (1 + s*T1) ->Out
 * 
 */

//#include "boost/smart_ptr/shared_ptr.hpp"
#include "LeadLagBlock.hpp"
#include <cstdio>

/**
 * Basic constructor
 */
gridpack::dynamic_simulation::LeadLagBlock::LeadLagBlock(void)
{
	
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::LeadLagBlock::~LeadLagBlock(void)
{
}

double gridpack::dynamic_simulation::LeadLagBlock::init(double dOut, double T11, double T21, double t_inc1)
{
	T1 = T11;
	T2 = T21;
	double dIn;
	
	dIn = dOut;
	
	if (T1<4.0*t_inc1){
		x0 = 0.0;
	}else{
		x0 = dOut*(1.0 - T2/T1);
	}
	
	x1 = x0;
	dx0 = 0.0;
	dx1 = 0.0;
	
	return dIn;
}

double gridpack::dynamic_simulation::LeadLagBlock::predictor(double In, double t_inc, bool flag)
{
    if (!flag) {
		x0 = x1;
	}  
	
	//printf("--- debug LeadLagBlock predictor test 1: In = %f, x = %f, T2 = %f, T1 = %f \n", In, x0, T2, T1);
	
	if (T1<4.0*t_inc){
		dx0 = 0.0;
	}else{
		dx0 = ( In * (1.0 - T2/T1) - x0 ) / T1;

	} // finished dx compuatation
	
	//compute output
	double dOut;
	if (T1<4.0*t_inc){
		dOut = In;
	}else{
		dOut = In*T2/T1 + x0;
	}
	
	x1 = x0 + dx0 * t_inc;
	
	return dOut;
}

double gridpack::dynamic_simulation::LeadLagBlock::corrector(double In, double t_inc, bool flag)
{	
	if (T1<4.0*t_inc){
		dx1 = 0.0;
	}else{
		dx1 = ( In * (1.0 - T2/T1) - x1 ) / T1;

	} // finished dx compuatation
	
	//compute output
	double dOut;
	if (T1<4.0*t_inc){
		dOut = In;
	}else{
		dOut = In*T2/T1 + x1;
	}
	
	x1 = x0 + (dx0 + dx1) / 2.0 * t_inc;
	
	return dOut;
}