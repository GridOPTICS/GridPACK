/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   repca1.cpp
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
#include "base_plant_model.hpp"
#include "repca1.hpp"

#define TS_THRESHOLD 4

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Repca1Model::Repca1Model(void)
{
  bmodel_debug = false;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Repca1Model::~Repca1Model(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to
 * Repca1Model
 */
void gridpack::dynamic_simulation::Repca1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{

  if (!data->getValue(GENERATOR_REPCA_KC,    &kc     ,idx))   kc =    0.02;
  if (!data->getValue(GENERATOR_REPCA_TFLTR, &tfltr  ,idx))   tfltr = 0.02;
  if (!data->getValue(GENERATOR_REPCA_DBD1,  &dbd1   , idx))  dbd1 =  -0.0;
  if (!data->getValue(GENERATOR_REPCA_DBD2,  &dbd2   , idx))  dbd2  = 0.0;
  if (!data->getValue(GENERATOR_REPCA_EMAX,  &emax   , idx))  emax  = 0.3;
  if (!data->getValue(GENERATOR_REPCA_EMIN,  &emin   , idx))  emin   = -0.3;
  if (!data->getValue(GENERATOR_REPCA_QMAX,  &qmax   , idx))  qmax  = 0.56;
  if (!data->getValue(GENERATOR_REPCA_QMIN,  &qmin   , idx))  qmin  = -0.56;
  if (!data->getValue(GENERATOR_REPCA_KP,    &kp     , idx))  kp   =  18.0;
  if (!data->getValue(GENERATOR_REPCA_KI,    &ki     , idx))  ki    = 5.0;
  if (!data->getValue(GENERATOR_REPCA_TFT,   &tft    , idx))  tft    = 0.0;
  if (!data->getValue(GENERATOR_REPCA_TFV,   &tfv    , idx))  tfv    = 0.05;
  if (!data->getValue(GENERATOR_REPCA_FDBD1, &fdbd1   , idx)) fdbd1  = 0.0;
  if (!data->getValue(GENERATOR_REPCA_FDBD2, &fdbd2   , idx)) fdbd2  = 0.0;
												               
  if (!data->getValue(GENERATOR_REPCA_DDN,   &ddn    , idx))  ddn    = 20.0;
  if (!data->getValue(GENERATOR_REPCA_DUP,   &dup    , idx))  dup    =  -10.0;
  if (!data->getValue(GENERATOR_REPCA_TP,    &tp     , idx))  tp     = 0.05;
  if (!data->getValue(GENERATOR_REPCA_FEMAX, &femax   , idx)) femax  = 999.0;
  if (!data->getValue(GENERATOR_REPCA_FEMIN, &femin   , idx)) femin  = -999.0;
  if (!data->getValue(GENERATOR_REPCA_KPG,   &kpg    , idx))  kpg    = 0.1; 
  if (!data->getValue(GENERATOR_REPCA_KIG,   &kig    , idx))  kig    = 0.05;
  if (!data->getValue(GENERATOR_REPCA_PMAX,  &pmax   , idx))  pmax   = 1.5;
  if (!data->getValue(GENERATOR_REPCA_PMIN,  &pmin   , idx))  pmin   =  -1.5;
  if (!data->getValue(GENERATOR_REPCA_TG,    &tg     , idx))  tg     = 0.1;
  
  if (bmodel_debug){
	printf("----!!renke debug:  repca1 model load() function: bus, %d,  kc=%12.6f, tfltr=%12.6f, dbd1=%12.6f, dbd2=%12.6f, emax=%12.6f, emin=%12.6f, qmax=%12.6f, qmin=%12.6f, kp=%12.6f, ki=%12.6f, tft=%12.6f, tfv =%12.6f \n", p_bus_id, kc, tfltr, dbd1, dbd2, emax, emin, qmax, qmin, kp, ki, tft, tfv);
	printf("----!!renke debug:  repca1 model load() function: fdbd1=%12.6f, fdbd2=%12.6f, ddn=%12.6f, dup=%12.6f, tp=%12.6f, femax=%12.6f, femin=%12.6f \n", fdbd1, fdbd2, ddn, dup, tp, femax, femin);   
	printf("----!!renke debug:  repca1 model load() function: kpg=%12.6f, kig=%12.6f, pmax=%12.6f, pmin=%12.6f, tg=%12.6f, \n " , kpg, kig, pmax, pmin, tg); 
  }
  
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Repca1Model::init(double mag, double ang, double ts)
{
  //
  double tmpout;
  tmpout = delayblock_x7pref.init(pref, tg);
  //printf ("----!!!!repc debug init test 1: tmpout: %f, kpg: %f, kig: %f, pmax: %f, pmin: %f\n", tmpout, kpg, kig, pmax, pmin);
  tmpout = piblock_x6ppi.init(tmpout, kpg, kig, pmax, pmin);
  //printf ("----!!!!repc debug init test 2: tmpout: %f, piblock_x6ppi.x0: %f\n", tmpout, piblock_x6ppi.x0);
  
  double pdelaytmp = delayblock_x5pmeas.init(genP, tp);
  plant_pref   = tmpout+ pdelaytmp;
  
  freqref = 1.0;
  busfreq = 1.0;
  
  tmpout = leadlagblock_x4qext.init(qext, tfv, tft, ts);
  tmpout = piblock_x3qpi.init(tmpout, kp, ki, qmax, qmin);
  
  double qdelaytmp = delayblock_x1vmeas.init(genQ*kc+Vterm, tfltr);
  
  vref = tmpout + qdelaytmp;
  
  if (bmodel_debug){
	printf("----renke debug: Repca1Model init states:  bus, %d, x1vmeas: %12.6f,  piblock_x3qpi: %12.6f,  leadlagblock_x4qext: %12.6f,  delayblock_x5pmeas: %12.6f,  piblock_x6ppi: %12.6f,  delayblock_x7pref: %12.6f \n", 
	p_bus_id, delayblock_x1vmeas.x0, piblock_x3qpi.x0, leadlagblock_x4qext.x0, 
	      delayblock_x5pmeas.x0, piblock_x6ppi.x0, delayblock_x7pref.x0); 
		  
    printf("----renke debug: Repca1Model init states:  bus, %d, pref, %12.6f, qext, %12.6f, vref, %12.6f, freqref, %12.6f, plant_pref, %12.6f, genP, %12.6f, genQ, %2.6f, Vterm, %12.6f \n", 
	p_bus_id, pref, qext, vref, freqref, plant_pref, genP, genQ, Vterm); 
  }
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Repca1Model::predictor(double t_inc, bool flag)
{
	/*
	if (!flag) {
    x1pss = x1pss_1;
    x2pss = x2pss_1;
    x3pss = x3pss_1;
    }
	*/
	
	//-----Qexit flow path-------------------
	
	double tmpin = genQ*kc + Vterm;
	//printf ("----!!!!repc debug Repca1Model predictor test 1: genQ: %f, kc: %15.11f, Vterm: %15.11f tmpin: %15.11f\n", genQ, kc, Vterm, tmpin);
	
	double tmpout;
	tmpout = delayblock_x1vmeas.predictor(tmpin, t_inc, flag);
	//printf ("----!!!!repc debug Repca1Model predictor test 2: tmpout:  %15.11f \n", tmpout);
	
	tmpout = vref - tmpout;
	//printf ("----!!!!repc debug Repca1Model predictor test 3: vref: %15.11f, tmpout:  %15.11f \n", vref, tmpout);
	
	//----------pass the deadband ---------add later
	
	if (tmpout>emax) tmpout = emax;
	if (tmpout<emin) tmpout = emin;
	
	//printf ("----!!!!repc debug Repca1Model predictor test 31: tmpout:  %15.11f, emax:  %15.11f, emin:  %15.11f \n", tmpout, emax, emin);
	
	tmpout = piblock_x3qpi.predictor(tmpout, t_inc, flag);
	
	//printf ("----!!!!repc debug Repca1Model predictor test 4: tmpout:  %15.11f \n",  tmpout);
	
	qext = leadlagblock_x4qext.predictor(tmpout, t_inc, flag);
	
	//printf ("----!!!!repc debug Repca1Model predictor test 5: tmpout:  %15.11f \n",  qext);
	
	//----Pref flow path---------------
	tmpin = freqref - busfreq;
	
	//----------pass the deadband ---------add later
		
	if (tmpin >= 0.0){
		tmpin = dup*tmpin;
	}else
	{
		tmpin = ddn*tmpin;
	}
	
	double pouttmp = delayblock_x5pmeas.predictor(genP, t_inc, flag);
	tmpin = tmpin + plant_pref - pouttmp;
	
	if (tmpin>femax) tmpin = femax;
	if (tmpin<femin) tmpin = femin;
	
	tmpout = piblock_x6ppi.predictor(tmpin, t_inc, flag);
	pref = delayblock_x7pref.predictor(tmpout, t_inc, flag);
	
	if (bmodel_debug){
		printf("----renke debug: repca1 predictor dx:  bus, %d, delayblock_x1vmeas %12.6f,  piblock_x3qpi %12.6f,  leadlagblock_x4qext %12.6f,  delayblock_x5pmeas %12.6f,  piblock_x6ppi %12.6f,  delayblock_x7pref %12.6f \n", 
		p_bus_id, delayblock_x1vmeas.dx0, piblock_x3qpi.dx0, leadlagblock_x4qext.dx0, 
	      delayblock_x5pmeas.dx0, piblock_x6ppi.dx0, delayblock_x7pref.dx0); 
		  
		printf("----renke debug test 1: repca1 predictor x:  bus, %d,  %12.6f,   %12.6f,   %12.6f,   %12.6f,   %12.6f,   %12.6f \n", 
		p_bus_id, delayblock_x1vmeas.x0, piblock_x3qpi.x0, leadlagblock_x4qext.x0, 
	      delayblock_x5pmeas.x0, piblock_x6ppi.x0, delayblock_x7pref.x0); 
		 
		printf("----renke debug: repca1 predictor other:  bus, %d, pref = %12.6f, qext = %12.6f, genP = %12.6f, genQ = %12.6f, busfreq = %12.6f, Vterm = %12.6f, \n", 
		p_bus_id, pref, qext, genP, genQ, busfreq, Vterm);
		
		printf("----renke debug test 2: repca1 predictor other:  bus, %d,  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f,  %12.6f, \n", 
		p_bus_id, pref, qext, genP, genQ, busfreq, Vterm);
	}
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Repca1Model::corrector(double t_inc, bool flag)
{
    //-----Qexit flow path-------------------
	double tmpin = genQ*kc + Vterm;
	if (bmodel_debug) printf ("----!!!!repc debug Repca1Model predictor test 1: genQ: %f, kc: %15.11f, Vterm: %15.11f tmpin: %15.11f\n", genQ, kc, Vterm, tmpin);
		
	double tmpout;
	tmpout = delayblock_x1vmeas.corrector(tmpin, t_inc, flag);
	if (bmodel_debug) printf ("----!!!!repc debug Repca1Model predictor test 2: tmpout:  %15.11f \n", tmpout);
		
	tmpout = vref - tmpout;
	if (bmodel_debug) printf ("----!!!!repc debug Repca1Model predictor test 3: vref: %15.11f, tmpout:  %15.11f \n", vref, tmpout);
	
	//----------pass the deadband ---------add later
	
	if (tmpout>emax) tmpout = emax;
	if (tmpout<emin) tmpout = emin;
	if (bmodel_debug) printf ("----!!!!repc debug Repca1Model predictor test 31: tmpout:  %15.11f, emax:  %15.11f, emin:  %15.11f \n", tmpout, emax, emin);
		
	tmpout = piblock_x3qpi.corrector(tmpout, t_inc, flag);
	if (bmodel_debug) printf ("----!!!!repc debug Repca1Model predictor test 4: tmpout:  %15.11f \n",  tmpout);
	
	qext = leadlagblock_x4qext.corrector(tmpout, t_inc, flag);
	if (bmodel_debug) printf ("----!!!!repc debug Repca1Model predictor test 5: tmpout:  %15.11f \n",  qext);
	
	//----Pref flow path---------------
	tmpin = freqref - busfreq;
	
	//----------pass the deadband ---------add later
		
	if (tmpin >= 0.0){
		tmpin = ddn*tmpin;
	}else
	{
		tmpin = dup*tmpin;
	}
	
	double pouttmp = delayblock_x5pmeas.corrector(genP, t_inc, flag);
	tmpin = tmpin + plant_pref - pouttmp;
	
	if (tmpin>femax) tmpin = femax;
	if (tmpin<femin) tmpin = femin;
	
	tmpout = piblock_x6ppi.corrector(tmpin, t_inc, flag);
	pref = delayblock_x7pref.corrector(tmpout, t_inc, flag);
	
	if (bmodel_debug){
		printf("----renke debug: repca1 corrector dx:  bus, %d, delayblock_x1vmeas %12.6f,  piblock_x3qpi %12.6f,  leadlagblock_x4qext %12.6f,  delayblock_x5pmeas %12.6f,  piblock_x6ppi %12.6f,  delayblock_x7pref %12.6f \n", 
		p_bus_id, delayblock_x1vmeas.dx1, piblock_x3qpi.dx1, leadlagblock_x4qext.dx1, 
	    delayblock_x5pmeas.dx1, piblock_x6ppi.dx1, delayblock_x7pref.dx1); 
		
		printf("----renke debug test 1: repca1 corrector x:  bus, %d,  %12.6f,   %12.6f,   %12.6f,   %12.6f,   %12.6f,   %12.6f \n", 
		p_bus_id, delayblock_x1vmeas.x1, piblock_x3qpi.x1, leadlagblock_x4qext.x1, 
	      delayblock_x5pmeas.x1, piblock_x6ppi.x1, delayblock_x7pref.x1);
		 
		printf("----renke debug: repca1 corrector other:  bus, %d, pref = %12.6f, qext = %12.6f, genP = %12.6f, genQ = %12.6f, busfreq = %12.6f, Vterm = %12.6f, \n", 
		p_bus_id, pref, qext, genP, genQ, busfreq, Vterm);
	}
}

void gridpack::dynamic_simulation::Repca1Model::setGenPQV(double P, double Q, double Vt)
{
	genP = P;
	genQ = Q;
	Vterm = Vt;
}

void gridpack::dynamic_simulation::Repca1Model::setBusFreq(double freq)
{
	busfreq  = freq;
}

void gridpack::dynamic_simulation::Repca1Model::setPrefQext(double Pref, double Qext)
{
	pref = Pref;
	qext = Qext;
}

double gridpack::dynamic_simulation::Repca1Model::getPref( )
{
	return pref;
}

double gridpack::dynamic_simulation::Repca1Model::getQext( )
{
	return qext;
}

void gridpack::dynamic_simulation::Repca1Model::setExtBusNum(int ExtBusNum)
{
	p_bus_id = ExtBusNum;
}	

void gridpack::dynamic_simulation::Repca1Model::setExtGenId(std::string ExtGenId)
{
	p_ckt = ExtGenId;
}


