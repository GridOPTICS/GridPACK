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
#include <cstdio>

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

double gridpack::dynamic_simulation::DelayBlockwithLimit::init(double dOut, double Ts1, double Max1, double Min1)
{
	Ts = Ts1;
	Max = Max1;
	Min = Min1;
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
	//printf ("----delayblock test 1: x0: %15.11f, x1: %15.11f,  Ts: %15.11f, t_inc: %15.11f \n", x0, x1, Ts, t_inc);
    if (!flag) {
		x0 = x1;
	}  
	
	//printf ("----delayblock test 2: x0: %15.11f \n", x0);
	
	if (Ts < 4.0*t_inc){
		x0 = In;
	}
	
	//printf ("----delayblock test 3: x0: %15.11f \n", x0);
	if (x0>Max) x0 = Max;
	if (x0<Min) x0 = Min;
	
	//printf ("----delayblock test 4: x0: %15.11f, max: %f, min: %f \n", x0, Max, Min);
	
	if (Ts < 4.0*t_inc){
		dx0 = 0.0;
	}else{
		dx0 = (In - x0)/Ts;
	}
	
	//printf ("----delayblock test 5: dx0: %15.11f  \n", dx0);
	
	if ( dx0 > 0.0 && x0 >= Max) dx0 = 0.0;
	if ( dx0 < 0.0 && x0 <= Min) dx0 = 0.0;
	
	//printf ("----delayblock test 6: dx0: %15.11f  \n", dx0);
	// finished dx compuatation
	
	//compute output
	double dOut = x0;
	
	//printf ("----delayblock test 7: dOut: %15.11f  \n", dOut);
	if (dOut > Max) dOut = Max;
	if (dOut < Min) dOut = Min;
	
	//printf ("----delayblock test 8: dOut: %15.11f  \n", dOut);
	
	x1 = x0 + dx0 * t_inc;
	
	//printf ("----delayblock test 9: x1: %15.11f  \n", x1);
	
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
	if ( dx1 < 0.0 && x1 <= Min) dx1 = 0.0;
	// finished dx compuatation
	
	//compute output
	double dOut = x1;
	if (dOut > Max) dOut = Max;
	if (dOut < Min) dOut = Min;
	
	x1 = x0 + (dx0 + dx1) / 2.0 * t_inc;
	
	return dOut;
}