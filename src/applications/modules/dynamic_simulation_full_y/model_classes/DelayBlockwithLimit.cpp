/*
 *     Copyright (c) 2021 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   DelayBlockwithLimit.cpp
 * @author Renke Huang renke.huang@pnnl.gov
 * @Last modified:   May 7, 2021
 * 
 * @brief  
 * 
 * 
 */

//#include "boost/smart_ptr/shared_ptr.hpp"
#include "DelayBlockwithLimit.hpp"

/**
 * Basic constructor
 */
gridpack::dynamic_simulation::DelayBlockwithLimit::DelayBlockwithLimit(void)
{
	
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::DelayBlockwithLimit::~DelayBlockwithLimit(void)
{
}

double gridpack::dynamic_simulation::DelayBlockwithLimit::init(double dOut, double Ts, double Max, double Min)
{
	Ts = Ts;
	Max = Max;
	Min = Min;
	double dIn;
	dIn = dOut;
	x0 = dOut;
	
	//adjust the Max min values to make sure flat start
	if (dOut>Max) Max = dOut + 0.1;
	if (dOut<Min) Min = dOut - 0.1;
	
	x1 = x0;
	dx0 = 0.0;
	dx1 = 0.0;
	
	return dIn;
}

double gridpack::dynamic_simulation::DelayBlockwithLimit::predictor(double In, double t_inc, bool flag)
{
    if (!flag) {
		x0 = x1;
	}  
	
	if (Ts < 4.0*t_inc){
		x0 = In;
	}
	if (x0>Max) x0 = Max;
	if (x0<Min) x0 = Min;
	
	if (Ts < 4.0*t_inc){
		dx0 = 0.0;
	}else{
		dx0 = (In - x0)/Ts;

	}
	
	if ( dx0 > 0.0 && x0 >= Max) dx0 = 0.0;
	if ( dx0 < 0.0 && x0 >= Min) dx0 = 0.0;
	// finished dx compuatation
	
	//compute output
	double dOut = x0;
	if (dOut > Max) dOut = Max;
	if (dOut < Min) dOut = Min;
	
	x1 = x0 + dx0 * t_inc;
	
	return dOut;
}

double gridpack::dynamic_simulation::DelayBlockwithLimit::corrector(double In, double t_inc, bool flag)
{
	
	if (Ts < 4.0*t_inc){
		x1 = In;
	}
	if (x1>Max) x1 = Max;
	if (x1<Min) x1 = Min;
	
	if (Ts < 4.0*t_inc){
		dx1 = 0.0;
	}else{
		dx1 = (In - x1)/Ts;

	}
	
	if ( dx1 > 0.0 && x1 >= Max) dx1 = 0.0;
	if ( dx1 < 0.0 && x1 >= Min) dx1 = 0.0;
	// finished dx compuatation
	
	//compute output
	double dOut = x1;
	if (dOut > Max) dOut = Max;
	if (dOut < Min) dOut = Min;
	
	x1 = x0 + (dx0 + dx1) / 2.0 * t_inc;
	
	return dOut;
}