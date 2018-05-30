/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   ieel.cpp
 * @author Shuangshuang Jin 
 * @Last modified:  Oct 25, 2016 
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
#include "base_load_model.hpp"
#include "ieel.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::IeelLoad::IeelLoad(void)
{
  a1 = 0.0;
  a2 = 0.0;
  a3 = 0.0;
  a4 = 0.0;
  a5 = 0.0;
  a6 = 0.0;
  a7 = 0.0;
  a8 = 0.0;
       
  n1 = 0.0;
  n2 = 0.0;
  n3 = 0.0;
  n4 = 0.0;
  n5 = 0.0;
  n6 = 0.0;

  P0 = 0.0; // intial load P value
  Q0 = 0.0; // intial load Q value
  P = 0.0;  // actual load P at each time step
  Q = 0.0;  // actual load Q at each time step
       
  p_INorton = gridpack::ComplexType(0.0, 0.0); 
  nortonY = gridpack::ComplexType(0.0, 0.0); 
    
  vt_init = 1.0;
     
  Vsmin_pu = 0.7;  // constant P will be scaled down if voltage is less than Vsmin_pu
  Vimin_pu = 0.7;  // constant I will be scaled down if voltage is less than Vimin_pu
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::IeelLoad::~IeelLoad(void)
{
}

/**
 * Load parameters from DataCollection object into load model
 * @param data collection of load parameters from input files
 * @param index of load on bus
 */
void gridpack::dynamic_simulation::IeelLoad::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx, double dloadP, double dloadQ, int ibCMPL)
{
  data->getValue(BUS_NUMBER,&p_bus_id);
  data->getValue(LOAD_ID,&p_loadid,idx);
  
  // set the information of dynamic load P, Q, id at the base class level too.
  p_pl = dloadP;
  p_ql = dloadQ;
  setDynLoadP(p_pl);
  setDynLoadQ(p_ql);
  setDynLoadID(p_loadid);
  
  printf("IeelLoad::load(): p_pl: %12.6f, p_ql: %12.6f \n", p_pl, p_ql);
  
  if ( ibCMPL==1 ){  // if load from the composite load model
	if (!data->getValue(LOAD_P1C, &a1))  		a1 = 0.0;
	if (!data->getValue(LOAD_P2C, &a2))  		a2 = 0.0;
	
	a3 = 1.0 - a1-a2;
	//if (!data->getValue(LOAD_A3, &a3))  		a3 = 0.0;
	if (!data->getValue(LOAD_Q1C, &a4))  		a4 = 0.0;
	if (!data->getValue(LOAD_Q2C, &a5))  		a5 = 0.0;
	a6 = 1.0 - a4-a5;
	//if (!data->getValue(LOAD_A6, &a6))  		a6 = 0.0;
	if (!data->getValue(LOAD_PFREQ, &a7))  		a7 = 0.0;
	if (!data->getValue(LOAD_QFREQ, &a8))  		a8 = 0.0;
	
	if (!data->getValue(LOAD_P1E, &n1))  		n1 = 0.0;
	if (!data->getValue(LOAD_P2E, &n2))  		n2 = 0.0;
	n3 = 0.0;
	//if (!data->getValue(LOAD_N3, &n3))  		n3 = 0.0;
	if (!data->getValue(LOAD_Q1E, &n4))  		n4 = 0.0;
	if (!data->getValue(LOAD_Q2E, &n5))  		n5 = 0.0;
	n6 = 0.0;
	//if (!data->getValue(LOAD_N6, &n6))  		n6 = 0.0;
	
  }else{
	if (!data->getValue(LOAD_A1, &a1, idx))  		a1 = 0.0;
	if (!data->getValue(LOAD_A2, &a2, idx))  		a2 = 0.0;
	if (!data->getValue(LOAD_A3, &a3, idx))  		a3 = 0.0;
	if (!data->getValue(LOAD_A4, &a4, idx))  		a4 = 0.0;
	if (!data->getValue(LOAD_A5, &a5, idx))  		a5 = 0.0;
	if (!data->getValue(LOAD_A6, &a6, idx))  		a6 = 0.0;
	if (!data->getValue(LOAD_A7, &a7, idx))  		a7 = 0.0;
	if (!data->getValue(LOAD_A8, &a8, idx))  		a8 = 0.0;
	
	if (!data->getValue(LOAD_N1, &n1, idx))  		n1 = 0.0;
	if (!data->getValue(LOAD_N2, &n2, idx))  		n2 = 0.0;
	if (!data->getValue(LOAD_N3, &n3, idx))  		n3 = 0.0;
	if (!data->getValue(LOAD_N4, &n4, idx))  		n4 = 0.0;
	if (!data->getValue(LOAD_N5, &n5, idx))  		n5 = 0.0;
	if (!data->getValue(LOAD_N6, &n6, idx))  		n6 = 0.0;
  }
    
  printf("IeelLoad::load(): a1: %f, a2: %f, a3: %f, a4: %f, a5: %f, a6: %f, a7: %f, a8: %f, \n", a1, a2, a3, a4, a5, a6, a7, a8);
  printf("IeelLoad::load(): n1: %f, n2: %f, n3: %f, n4: %f, n5: %f, n6: %f,  \n", n1, n2, n3, n4, n5, n6);
  
}

/**
 * Initialize load model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 */
void gridpack::dynamic_simulation::IeelLoad::init(double mag,
    double ang, double ts)
{
  double vt_mag = mag;
  //SJin: What are parameters Pinit_pu, Qinit_pu, and vt? From where to get them?
  //Fake declaration; 
  double Pinit_pu, Qinit_pu, sysMVA;
  sysMVA = 100.0;
  Pinit_pu = p_pl/sysMVA;
  Qinit_pu = p_ql/sysMVA;

  P0 = Pinit_pu;
  Q0 = Qinit_pu;
  P = Pinit_pu;
  Q = Qinit_pu;

  vt_init = vt_mag;
  gridpack::ComplexType tmp (Pinit_pu, -Qinit_pu); 
  nortonY = tmp / vt_mag / vt_mag;

  // check the data to make sure a1+a2+a3 = 1 
  if (a1 + a2 + a3  > 0.0) {
    if (abs(a1 + a2 + a3 -1.0) >1.0E-6) {
      a1 = a1/ (a1 + a2 + a3);
      a2 = a2/ (a1 + a2 + a3);
      a3 = a3/ (a1 + a2 + a3);
    } else {
      a1 = 1.0;
      a2 = 0.0;
      a3 = 0.0;
    }
  }
          
  // check the data to make sure a4+a5+a6 = 1 
  if (a4 + a5 + a6  > 0.0) {
    if (abs(a4 + a5 + a6 - 1.0) >1.0E-6) {
      a4 = a4/ (a4 + a5 + a6);
      a5 = a5/ (a4 + a5 + a6);
      a6 = a6/ (a4 + a5 + a6);
    } else {
      a4 = 1.0;
      a5 = 0.0;
      a6 = 0.0;
    }
  }
}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType gridpack::dynamic_simulation::IeelLoad::INorton()
{
  //SJIN: Matlab getINorton?
  return p_INorton;
}


/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType gridpack::dynamic_simulation::IeelLoad::NortonImpedence()
{
  //SJIN: Matlab getNortonImpedance?
  return nortonY; // refer to gensal.cpp?
}

/**
 * Predict part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::IeelLoad::predictor_currentInjection(bool flag)
{
  //SJin: What are parameter vt_complex and freq? From where to get them?
  //Fake declaration; 
  double freq = presentFreq;
  
  double vt = abs(vt_complex);
  double v = vt/vt_init;
  gridpack::ComplexType PQ_equivY = vt*vt*conj(nortonY);

  // power scaling factor
  double kp = 1.0;
  double ki = 1.0;

  double pi = 4.0*atan(1.0); 
  if (vt < Vsmin_pu) 
    kp = 0.5*(1.0-cos(pi*vt/Vsmin_pu));
   
  if (vt < Vimin_pu)
    ki = sin(pi*vt/2.0/Vimin_pu);

  P = P0*(a1*pow(v,n1) + ki*a2*pow(v,n2) + kp*a3*pow(v,n3))*(1.0 + a7*(freq -1.0));
  Q = Q0*(a4*pow(v,n4) + ki*a5*pow(v,n5) + kp*a6*pow(v,n6))*(1.0 + a8*(freq -1.0));
   
  gridpack::ComplexType tmp(P, Q);
  gridpack::ComplexType dPQ = PQ_equivY - tmp;            
               
  // Inorton is to compensate the power difference
  p_INorton = conj(dPQ/vt_complex);

} 

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::IeelLoad::predictor(
    double t_inc, bool flag)
{
}

/**
 * Corrector part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::IeelLoad::corrector_currentInjection(bool flag)
{
  //SJin: What are parameter vt_complex and freq? From where to get them?
  //Fake declaration; 
  double freq = presentFreq;
  
  double vt = abs(vt_complex);
  double v = vt/vt_init;
  gridpack::ComplexType PQ_equivY = vt*vt*conj(nortonY);

  // power scaling factor
  double kp = 1.0;
  double ki = 1.0;

  double pi = 4.0*atan(1.0); 
  if (vt < Vsmin_pu) 
    kp = 0.5*(1.0-cos(pi*vt/Vsmin_pu));
   
  if (vt < Vimin_pu)
    ki = sin(pi*vt/2.0/Vimin_pu);

  P = P0*(a1*pow(v,n1) + ki*a2*pow(v,n2) + kp*a3*pow(v,n3))*(1.0 + a7*(freq -1.0));
  Q = Q0*(a4*pow(v,n4) + ki*a5*pow(v,n5) + kp*a6*pow(v,n6))*(1.0 + a8*(freq -1.0));
   
  gridpack::ComplexType tmp(P, Q);
  gridpack::ComplexType dPQ = PQ_equivY - tmp;            
               
  // Inorton is to compensate the power difference
  p_INorton = conj(dPQ/vt_complex);

}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::IeelLoad::corrector(
    double t_inc, bool flag)
{
}

/**
 * Set voltage on each load
 */
void gridpack::dynamic_simulation::IeelLoad::setVoltage(
    gridpack::ComplexType voltage)
{
  printf("IeelLoad::setVoltage, %12.6f + j %12.6f \n", real(voltage), imag(voltage));	
  presentMag = abs(voltage);
  presentAng = atan2(imag(voltage), real(voltage));  
  vt_complex = voltage;
  
}

/**
 * Set terminal voltage frequency on each load
 */
void gridpack::dynamic_simulation::IeelLoad::setFreq(double dFreq)
{
  presentFreq = dFreq;
}

/**
 * Write out load state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
bool gridpack::dynamic_simulation::IeelLoad::serialWrite(
    char* string, const int bufsize, const char *signal)
{
  return false;
}
