/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   reeca1.cpp
 * @author Renke Huang
 * @Last modified:   May. 16, 2021
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_exciter_model.hpp"
#include "reeca1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Reeca1Model::Reeca1Model(void)
{

	Vterm = 1.0;
	
	OptionToModifyLimitsForInitialStateLimitViolation = true;
	bmodel_debug = false;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Reeca1Model::~Reeca1Model(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::Reeca1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  p_sbase = 100.0;

  data->getValue(BUS_NUMBER,&p_bus_id);
  data->getValue(GENERATOR_ID,&p_ckt,idx);

  if (!data->getValue(GENERATOR_MBASE, &MVABase, idx)) MVABase = 100.0; // MVABase
  if (!data->getValue(GENERATOR_REECA_TRV , &trv, idx))     trv = 0.02; 
  if (!data->getValue(GENERATOR_REECA_DBD1 ,    &dbd1, idx))    dbd1 = -0.05; 
  if (!data->getValue(GENERATOR_REECA_DBD2 , &dbd2, idx))    dbd2 = 0.05; 
  if (!data->getValue(GENERATOR_REECA_VDIP , &vdip, idx))    vdip = -99.0; 
  if (!data->getValue(GENERATOR_REECA_VUP , &vup, idx))     vup = 99.0; 
  if (!data->getValue(GENERATOR_REECA_KQV , &kqv, idx))     kqv = 0.0; 
  if (!data->getValue(GENERATOR_REECA_IMAX , &imax, idx))    imax = 2.0; 
  if (!data->getValue(GENERATOR_REECA_TIQ , &tiq, idx))     tiq = 0.02; 
  if (!data->getValue(GENERATOR_REECA_DPMAX , &dpmax, idx))   dpmax = 999.0; 
  if (!data->getValue(GENERATOR_REECA_DPMIN , &dpmin, idx))   dpmin = -999.0; 
  if (!data->getValue(GENERATOR_REECA_PMAX , &pmax, idx))    pmax = 1.5; 
  if (!data->getValue(GENERATOR_REECA_PMIN ,   &pmin, idx))    pmin = 0.04; 
  if (!data->getValue(GENERATOR_REECA_TPORD , &tpord, idx))  tpord = 0.02; 
  ipmax = imax;
  ipmin = -imax;
  iqmax = imax;
  iqmin = -imax;
  
  if (bmodel_debug){
	printf("\n--------Reeca1Model parameters: MVABase = %12.6f, trv = %12.6f, dbd1 = %12.6f, dbd2 = %12.6f, vdip = %12.6f, vup = %12.6f, kqv = %12.6f, imax = %12.6f, tiq = %12.6f, dpmax = %12.6f, dpmin = %12.6f, pmax = %12.6f, pmin = %12.6f, tpord = %12.6f  \n", 
	MVABase, trv, dbd1, dbd2, vdip, vup, kqv, imax, tiq, dpmax, dpmin, pmax, pmin, tpord);
  }
  
}

/**
 * Initialize generator model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 */
void gridpack::dynamic_simulation::Reeca1Model::init(double mag,
    double ang, double ts)
{
  //printf("Step0 gen%d mag = %f\n", p_bus_id, mag);
    Vterm = mag;
    //presentMag = mag;
    //Theta = ang;
    //presentAng = ang;
    
    //assume already get ipcmd and iqcmd from generator model by set
    vfilt = delayblock_x1vfilt.init(Vterm, trv);
    
    double iqtmp = delayblock_x3qv.init( iqcmd, tiq);  
    qext = iqtmp * vfilt;
    
    double pordtmp = ipcmd*vfilt;
	//if (bmodel_debug) printf ("Reeca1Model::init, test 1 tpord: %f, pmax: %f, pmin: %f \n", tpord, pmax, pmin);
    pref = delayblocklimit_x2pord.init(pordtmp, tpord, pmax, pmin);
	
	//if (bmodel_debug) printf ("Reeca1Model::init, test 2 tpord: %f, pmax: %f, pmin: %f \n", delayblocklimit_x2pord.Ts, delayblocklimit_x2pord.Max, delayblocklimit_x2pord.Min);
  
    x1vfilt_0 = delayblock_x1vfilt.x0;
	x2pord_0 = delayblocklimit_x2pord.x0;
	x3qv_0 = delayblock_x3qv.x0;
	
    x1vfilt_1 = x1vfilt_0;
	x2pord_1 = x2pord_0;
	x3qv_1 = x3qv_0;
  
  if (bmodel_debug){
    printf("Reeca1Model init, bus: %d, Vterm = %f, Ipcmd = %f, Iqcmd = %f, pref= %f, qext= %f, \n", p_bus_id, Vterm, ipcmd, iqcmd, pref, qext);
	printf("Reeca1Model init, bus: %d, x1vfilt_0 = %f, x2pord_0 = %f, x3qv_0 = %f, \n", p_bus_id, x1vfilt_0, x2pord_0, x3qv_0);
  }

}


/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Reeca1Model::predictor(
    double t_inc, bool flag)
{
	
	//v filter
	vfilt =  delayblock_x1vfilt.predictor(Vterm, t_inc, flag);
	
	///--------Q path----------------------------
	double vfilt_clip = vfilt;
	if (vfilt <= 0.01) vfilt_clip = 0.01;
	
	double tmpin = qext/vfilt_clip;
	double tmpoutq = delayblock_x3qv.predictor(tmpin, t_inc, flag);
	
	double iqinj = 0.0;
	//-----------later add iqinj logic
	
	iqcmd = tmpoutq+iqinj;
	if (iqcmd > iqmax) iqcmd = iqmax;
	if (iqcmd < iqmin) iqcmd = iqmin;
	
	//-----------P path--------------------------
	
	//-----------later add the P ramping up and down limit
	//printf ("----!!!!repc debug Reeca1Model predictor test 0: delayblocklimit_x2pord.x0:  %15.11f,  \n", delayblocklimit_x2pord.x0);
	
	double tmppord = delayblocklimit_x2pord.predictor(pref, t_inc, flag);
	
	//printf ("----!!!!repc debug Reeca1Model predictor test 1: delayblocklimit_x2pord.x0:  %15.11f,  \n", delayblocklimit_x2pord.x0),
	
	ipcmd = tmppord/vfilt_clip;
	
	//printf ("----!!!!repc debug Reeca1Model predictor test 2: pref:  %15.11f, tmppord:  %15.11f, vfilt_clip: %15.11f, ipcmd, %15.11f, ipmax: %15.11f,, ipmin: %15.11f, \n", 
	//pref, tmppord, vfilt_clip, ipcmd, ipmax, ipmin);
	
	if (ipcmd > ipmax) ipcmd = ipmax;
	if (ipcmd < ipmin) ipcmd = ipmin;
	
	if (bmodel_debug){	
		printf("---Reeca1Model predictor dx, pbusid: %d, %12.6f, %12.6f, %12.6f\n", p_bus_id, delayblock_x1vfilt.dx0, delayblocklimit_x2pord.dx0, delayblock_x3qv.dx0);
	}
	
	if (bmodel_debug){
		printf("------renke debug in Reeca1Model::predictor other, pref = %12.6f, qext = %12.6f, vfilt_clip = %12.6f, ipcmd = %12.6f, iqcmd = %12.6f, \n", 
		pref, qext, vfilt_clip, ipcmd, iqcmd);
		printf("------Reeca1Model::predictor states: bus, %8d, %12.6f, %12.6f, %12.6f  \n",
          p_bus_id,  delayblock_x1vfilt.x0, delayblocklimit_x2pord.x0, delayblock_x3qv.x0);
	}
	
}


/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Reeca1Model::corrector(
    double t_inc, bool flag)
{
  	//v filter
	vfilt =  delayblock_x1vfilt.corrector(Vterm, t_inc, flag);
	
	///--------Q path----------------------------
	double vfilt_clip = vfilt;
	if (vfilt <= 0.01) vfilt_clip = 0.01;
	
	double tmpin = qext/vfilt_clip;
	double tmpoutq = delayblock_x3qv.corrector(tmpin, t_inc, flag);
	
	double iqinj = 0.0;
	//-----------later add iqinj logic
	
	iqcmd = tmpoutq+iqinj;
	if (iqcmd > iqmax) iqcmd = iqmax;
	if (iqcmd < iqmin) iqcmd = iqmin;
	
	//-----------P path--------------------------
	
	//-----------later add the P ramping up and down limit
	
	double tmppord = delayblocklimit_x2pord.corrector(pref, t_inc, flag);
	ipcmd = tmppord/vfilt_clip;
	
	if (ipcmd > ipmax) ipcmd = ipmax;
	if (ipcmd < ipmin) ipcmd = ipmin;
	
	if (bmodel_debug){	
		printf("---Reeca1Model corrector dx: pbusid: %d, %12.6f, %12.6f, %12.6f\n", p_bus_id, delayblock_x1vfilt.dx0, delayblocklimit_x2pord.dx0, delayblock_x3qv.dx0);
	}
	
	if (bmodel_debug){
		printf("------renke debug in Reeca1Model::corrector other, pref = %12.6f, qext = %12.6f, vfilt_clip = %12.6f, ipcmd = %12.6f, iqcmd = %12.6f, \n", 
		pref, qext, vfilt_clip, ipcmd, iqcmd);
		printf("------Reeca1Model::corrector states: bus, %8d, %12.6f, %12.6f, %12.6f  \n",
          p_bus_id,  delayblock_x1vfilt.x0, delayblocklimit_x2pord.x0, delayblock_x3qv.x0);
	}
  
}

/** 
 * Set the value of the Vterminal
 * @return value of field current
 */
void gridpack::dynamic_simulation::Reeca1Model::setVterminal(double mag)
{
  Vterm = mag;
}

void gridpack::dynamic_simulation::Reeca1Model::setIpcmdIqcmd(double ipcmd1, double iqcmd1)
{
	ipcmd = ipcmd1;
	iqcmd = iqcmd1;
}

void gridpack::dynamic_simulation::Reeca1Model::setPrefQext(double pref1, double qext1) 
{
	pref = pref1;
	qext = qext1;
}

double gridpack::dynamic_simulation::Reeca1Model::getPref()
{
	return pref;
}

double gridpack::dynamic_simulation::Reeca1Model::getQext()
{
	return qext;
}

double gridpack::dynamic_simulation::Reeca1Model::getIpcmd()
{
	return ipcmd;
}

double gridpack::dynamic_simulation::Reeca1Model::getIqcmd()
{
	return iqcmd;
}


// Yuan added below 2020-6-23
/** 
 * Set the exciter bus number
 * @return value of exciter bus number
 */
void gridpack::dynamic_simulation::Reeca1Model::setExtBusNum(int ExtBusNum)
{
	p_bus_id = ExtBusNum;
}	

/** 
 * Set the exciter generator id
 * @return value of generator id
 */
void gridpack::dynamic_simulation::Reeca1Model::setExtGenId(std::string ExtGenId)
{
	p_ckt = ExtGenId;
}	
// Yuan added above 2020-6-23

