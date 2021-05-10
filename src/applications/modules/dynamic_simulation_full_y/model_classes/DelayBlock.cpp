/*
 *     Copyright (c) 2021 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   DelayBlock.cpp
 * @author Renke Huang renke.huang@pnnl.gov
 * @Last modified:   May 7, 2021
 * 
 * @brief  
 * 
 * 
 */

//#include "boost/smart_ptr/shared_ptr.hpp"
#include "DelayBlock.hpp"

/**
 * Basic constructor
 */
gridpack::dynamic_simulation::DelayBlock::DelayBlock(void)
{
	
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::DelayBlock::~DelayBlock(void)
{
}

double gridpack::dynamic_simulation::DelayBlock::init(double dOut, double Ts)
{
	Ts = Ts;
	double dIn;
	dIn = dOut;
	x0 = dOut;
	
	x1 = x0;
	dx0 = 0.0;
	dx1 = 0.0;
	
	return dIn;
}

double gridpack::dynamic_simulation::DelayBlock::predictor(double In, double t_inc, bool flag)
{
    if (!flag) {
		x0 = x1;
	}  
	
	if (Ts < 4.0*t_inc){
		dx0 = 0.0;
		x0 = In;
	}else{
		dx0 = (In - x0)/Ts;

	} // finished dx compuatation
	
	//compute output
	double dOut = x0;
	
	x1 = x0 + dx0 * t_inc;
	
	return dOut;
}

double gridpack::dynamic_simulation::DelayBlock::corrector(double In, double t_inc, bool flag)
{
	
	if (Ts < 4.0*t_inc){
		dx1 = 0.0;
		x1 = In;
	}else{
		dx1 = (In - x1)/Ts;

	} // finished dx compuatation
	
	//compute output
	double dOut = x1;
	
	x1 = x0 + (dx0 + dx1) / 2.0 * t_inc;
	
	return dOut;
}