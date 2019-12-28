/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   psssim.hpp
 * @author Renke Huang
 * @Last modified:   Aug. 20, 2019
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_pss_model.hpp"
#include "psssim.hpp"

#define TS_THRESHOLD 4

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::PsssimModel::PsssimModel(void)
{
  dx1pss = 0.0;
  dx2pss = 0.0;
  dx3pss = 0.0;
  dx1pss_1 = 0.0;
  dx2pss_1 = 0.0;
  dx3pss_1 = 0.0;
  wideareafreq = 0.0;
  kp = 7.0/60.0;
  
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::PsssimModel::~PsssimModel(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to
 * PsssimModel
 */
void gridpack::dynamic_simulation::PsssimModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(PSSSIM_INPUTTYPE, &inputtype, idx)) inputtype = -1; // TBD: UEL
  if (!data->getValue(PSSSIM_BUS1, &bus1, idx)) bus1 = -1; // TBD: VOS
  if (!data->getValue(PSSSIM_BUS2, &bus2, idx)) bus2 = -1; // Tr
  if (!data->getValue(PSSSIM_BUS3, &bus3, idx)) bus3 = -1; // TBD: Vimax
  if (!data->getValue(PSSSIM_BUS4, &bus4, idx)) bus4 = -1; // TBD: Vimin
  if (!data->getValue(PSSSIM_BUS5, &bus5, idx)) bus5 = -1; // Tc
  if (!data->getValue(PSSSIM_BUS6, &bus6, idx)) bus6 = -1; // Tb
  if (!data->getValue(PSSSIM_GAINK, &gaink, idx)) gaink = 100.0; // TBD: Tc1
  if (!data->getValue(PSSSIM_TW, &tw, idx)) tw = 10.0; // TBD: Tb1
  if (!data->getValue(PSSSIM_T1, &t1, idx)) t1 = 0.05; // Ka
  if (!data->getValue(PSSSIM_T2, &t2, idx)) t2 = 0.015; // Ta
  if (!data->getValue(PSSSIM_T3, &t3, idx)) t3 = 0.08; // TBD: Vamax
  if (!data->getValue(PSSSIM_T4, &t4, idx)) t4 = 0.01; // TBD: Vamin
  if (!data->getValue(PSSSIM_MAXOUT, &maxout, idx)) maxout = 0.2; // Vrmax
  if (!data->getValue(PSSSIM_MINOUT, &minout, idx)) minout = -0.05; // Vrmin
  
  printf("----!!renke debug:  psssim load:  %d, %d, %d, %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f \n", inputtype, bus1, bus2, gaink, tw, t1, t2, t3, t4, maxout, minout); 

  //if (!data->getValue(EXCITER_TA1, &Ta1, idx)) Ta1 = 0.0; // Ta1
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::PsssimModel::init(double mag, double ang, double ts)
{
  x1pss = genspd;
  x2pss = 0.0;
  x3pss = 0.0;
  
  x1pss_1 = genspd;
  x2pss_1 = 0.0;
  x3pss_1 = 0.0;
  
  pssout_vstab = 0.0;
  psscon1 = t1/t2;
  psscon2 = t3/t4; 
  
//  if (t4!=0.0) {
//	 psscon2 = t3/t4; 
//  }else{
//	 psscon2 = 1.0; 
//  }

  printf("----renke debug: psssim init:  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f \n", x1pss_1, x2pss_1, x3pss_1, pssout_vstab, psscon1, psscon2); 
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::PsssimModel::predictor(double t_inc, bool flag)
{
	if (!flag) {
    x1pss = x1pss_1;
    x2pss = x2pss_1;
    x3pss = x3pss_1;
    }
  
	double var1, var2, var3;
	
	//double kp = 7.0/60;
	double addwidearea;
	
	if (inputtype == 1){
		addwidearea = kp*wideareafreq;
		printf ("-------------!renke debug: PsssimModel::predictor wide area freq: %12.6f, addwidearea: %12.6f \n", wideareafreq, addwidearea);
	}else{
		addwidearea = 0.0;
	}
	
	var1 = (genspd + addwidearea - x1pss)/tw;
	dx1pss = var1;
	
	var2 = psscon1*gaink*var1 + x2pss;
	dx2pss = ( (1.0-psscon1)*gaink*var1 - x2pss )/t2;
	
	var3 = psscon2*var2 + x3pss;
    dx3pss = ( (1.0-psscon2)*var2 - x3pss )/t4;
	
	x1pss_1 = x1pss + dx1pss * t_inc;
	x2pss_1 = x2pss + dx2pss * t_inc;
	x3pss_1 = x3pss + dx3pss * t_inc;

	//pssout_vstab = min(maxout,  max(var3,-minout));
	double tmp;
	if (var3 > (-maxout)){
		tmp = var3;
	}else {
		tmp = -maxout;
	}
	
	if (tmp < maxout){
		pssout_vstab = tmp;
	}else {
		pssout_vstab = maxout;
	} 

  //printf("----renke debug: psssim predictor:  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f \n", x1pss_1, x2pss_1, x3pss_1, dx1pss, dx2pss, dx3pss, pssout_vstab, var1, var2, var3); 
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::PsssimModel::corrector(double t_inc, bool flag)
{
	double var1, var2, var3;
	
	//double kp = -7.0/60;
	double addwidearea;
	
	if (inputtype == 1){
		addwidearea = kp*wideareafreq;
		printf ("-------------!renke debug: PsssimModel::corrector: wide area freq: %12.6f, addwidearea: %12.6f \n", wideareafreq, addwidearea);
	}else{
		addwidearea = 0.0;
	}

	var1 = (genspd + addwidearea - x1pss_1)/tw;
	dx1pss_1 = var1;
	
	var2 = psscon1*gaink*var1 + x2pss_1;
	dx2pss_1 = ( (1.0-psscon1)*gaink*var1 - x2pss_1 )/t2;
	
	var3 = psscon2*var2 + x3pss_1;
    dx3pss_1 = ( (1.0-psscon2)*var2 - x3pss_1 )/t4;
	
	x1pss_1 = x1pss + (dx1pss + dx1pss_1) / 2.0 * t_inc;
	x2pss_1 = x2pss + (dx2pss + dx2pss_1) / 2.0 * t_inc;
	x3pss_1 = x3pss + (dx3pss + dx3pss_1) / 2.0 * t_inc;
	
	//pssout_vstab = min(maxout,  max(var3,-minout));
	double tmp;
	if (var3 > (-maxout)){
		tmp = var3;
	}else {
		tmp = -maxout;
	}
	
	if (tmp < maxout){
		pssout_vstab = tmp;
	}else {
		pssout_vstab = maxout;
	} 

  //printf("psssim pssout_vstab: %f\n", pssout_vstab);
}

double gridpack::dynamic_simulation::PsssimModel::getVstab( )
{
	if (inputtype == -1){
		return 0.0;
	}else{
		return pssout_vstab;
	}

}

void gridpack::dynamic_simulation::PsssimModel::setOmega(double omega)
{
   genspd = omega + 1.0;
}

double gridpack::dynamic_simulation::PsssimModel::getBusFreq(int busnum)
{
	return 0.0;
}

void gridpack::dynamic_simulation::PsssimModel::setWideAreaFreqforPSS(double freq)
{
	wideareafreq = freq;
}

