/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   tgov1.cpp
 * @author Renke Huang
 * @Last modified:   June 17, 2020
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
#include "base_governor_model.hpp"
#include "tgov1.hpp"

#define TS_THRESHOLD 4

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Tgov1Model::Tgov1Model(void)
{
  bdebugmodel = false;
  
  x1pow = 1.0;
  x2val = 1.0;
  x1pow_1 = 1.0;
  x2val_1 = 1.0;
  
  dx1pow = 0.0;
  dx2val = 0.0;
  dx1pow_1 = 0.0;
  dx2val_1 = 0.0;
  
  Pmech = 1.0;
  Pref = 1.0;
  delta_w = 0.0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Tgov1Model::~Tgov1Model(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to
 * Tgov1Model
 */
void gridpack::dynamic_simulation::Tgov1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
   //bdebugmodel = true;
   
  if (bdebugmodel) printf ("----------------!!! DEBUG for TGOV1, parameters are:   \n");
  
  if (!data->getValue(GOVERNOR_R, &R, idx)) R = 0.05; 
  if (bdebugmodel) printf ("R = %8.4f \n", R);
  
  if (!data->getValue(GOVERNOR_T1, &T1, idx)) T1 = 0.5; 
  if (bdebugmodel) printf ("T1 = %8.4f \n", T1);
  
  if (!data->getValue(GOVERNOR_VMAX, &Vmax, idx)) Vmax = 1.0; 
  if (bdebugmodel) printf ("VMAX = %8.4f \n", Vmax);
  
  if (!data->getValue(GOVERNOR_VMIN, &Vmin, idx)) Vmin = 0.0; 
  if (bdebugmodel) printf ("VMIN = %8.4f \n", Vmin);
  
  if (!data->getValue(GOVERNOR_T2, &T2, idx)) T2 = 3.0; 
  if (bdebugmodel) printf ("T2 = %8.4f \n", T2);
  
  if (!data->getValue(GOVERNOR_T3, &T3, idx)) T3 = 10.0; 
  if (bdebugmodel) printf ("T3= %8.4f \n", T3);
  
  if (!data->getValue(GOVERNOR_DT, &Dt, idx)) Dt = 0.0; 
  if (bdebugmodel) printf ("Dt = %8.4f \n", Dt);
  
  //bdebugmodel = false;
  
  //I do not think we need the MVA BASE of the corresponding generator, all the values computed in this TGOV1 model is per-unit
  //if (!data->getValue(GENERATOR_MBASE, &GenMVABase, idx)) GenMVABase = 0.0;
  //printf ("Mbase = %8.4f \n", GenMVABase);

}

/**
 * Initialize governor model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Tgov1Model::init(double mag, double ang, double ts)
{
  if (bdebugmodel)  printf("----------init, Tgov1: Pmech = %f\n", Pmech);
  
  // State 1
  if (T3 < TS_THRESHOLD * ts){
		x1pow = 0.0;
  }else{
		x1pow = Pmech*(1.0 - T2/T3);
  }
  
  // State 2
  x2val = Pmech;
  
  if (x2val>=Vmax) printf("-------------warning!!!!!!!!!, TGOV1 model init error: x2val with value %12.6f is bigger than Vmax %12.6f \n", x2val, Vmax);
  if (x2val<=Vmin) printf("-------------warning!!!!!!!!!, TGOV1 model init error: x2val with value %12.6f is smaller than Vmin %12.6f \n", x2val, Vmin);
  
  Pref = Pmech*R;
  
  x1pow_1 = x1pow;
  x2val_1 = x2val;
  
  if (bdebugmodel)  printf("-----------init, Tgov1 init, x1pow: %12.6f, x2val: %12.6f, Pref: %12.6f \n", x1pow, x2val, Pref);
  
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Tgov1Model::predictor(double t_inc, bool flag)
{
  if (!flag) {
	x1pow = x1pow_1;
	x2val = x2val_1;
  }

  if (bdebugmodel)  printf("--------------debug Tgov1Model::predictor, Tgov1 Pref: %12.6f, delta_w: %12.6f \n", Pref, delta_w);
  double TempIn = (Pref - delta_w)/R;
  double TempOut;
  
  // State 2
  
  if (T1 < TS_THRESHOLD * t_inc) x2val = TempIn;
  if (x2val>Vmax) x2val = Vmax;
  if (x2val<Vmin) x2val = Vmin;

  if (T1 < TS_THRESHOLD * t_inc){
	  dx2val = 0.0;
  }else{
	  dx2val = (TempIn-x2val)/T1;
  }
  
  if ( (dx2val>0.0) && (x2val >= Vmax) ) dx2val = 0.0;
  if ( (dx2val<0.0) && (x2val <= Vmin) ) dx2val = 0.0;
  
  TempIn = x2val;
  if (TempIn > Vmax) TempIn = Vmax;
  if (TempIn < Vmin) TempIn = Vmin;
  
  // State 1
  
  if (T3 < TS_THRESHOLD * t_inc){
	  dx1pow = 0.0;
  }else{
	  dx1pow = ( TempIn*(1.0 - T2/T3) - x1pow )/T3;
  }
  
  if (T3 < TS_THRESHOLD * t_inc){
	  TempOut = TempIn;
  }else{
	  TempOut = TempIn*T2/T3 + x1pow;
  }
  
  Pmech = TempOut - delta_w * Dt;

  x1pow_1 = x1pow + dx1pow * t_inc;
  x2val_1 = x2val + dx2val * t_inc;

  if (bdebugmodel)  printf("--------------debug Tgov1Model::predictor, Tgov1 dx: dx1pow: %12.6f, dx2val: %12.6f \n", dx1pow, dx2val);
  if (bdebugmodel)  printf("--------------debug Tgov1Model::predictor, Tgov1 x: x1pow_1: %12.6f, x2val_1: %12.6f\n", x1pow_1, x2val_1);
  
  if (bdebugmodel)  printf("--------------debug Tgov1Model::predictor, Tgov1 Pmech = %12.6f\n", Pmech);
  
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Tgov1Model::corrector(double t_inc, bool flag)
{
  double TempIn = (Pref - delta_w)/R;
  double TempOut;
  
  // State 2
  
  if (T1 < TS_THRESHOLD * t_inc) x2val_1 = TempIn;
  if (x2val_1>Vmax) x2val_1 = Vmax;
  if (x2val_1<Vmin) x2val_1 = Vmin;

  if (T1 < TS_THRESHOLD * t_inc){
	  dx2val_1 = 0.0;
  }else{
	  dx2val_1 = (TempIn-x2val_1)/T1;
  }
  
  if ( (dx2val_1>0.0) && (x2val_1 >= Vmax) ) dx2val_1 = 0.0;
  if ( (dx2val_1<0.0) && (x2val_1 <= Vmin) ) dx2val_1 = 0.0;
  
  TempIn = x2val_1;
  if (TempIn > Vmax) TempIn = Vmax;
  if (TempIn < Vmin) TempIn = Vmin;
  
  // State 1
  
  if (T3 < TS_THRESHOLD * t_inc){
	  dx1pow_1 = 0.0;
  }else{
	  dx1pow_1 = ( TempIn*(1.0 - T2/T3) - x1pow_1 )/T3;
  }
  
  if (T3 < TS_THRESHOLD * t_inc){
	  TempOut = TempIn;
  }else{
	  TempOut = TempIn*T2/T3 + x1pow_1;
  }
  
  Pmech = TempOut - delta_w * Dt;

  x1pow_1 = x1pow + (dx1pow + dx1pow_1) / 2.0 * t_inc;
  x2val_1 = x2val + (dx2val + dx2val_1) / 2.0 * t_inc;
 

  if (bdebugmodel)  printf("--------------debug Tgov1Model::corrector, Tgov1 dx: dx1pow_1: %12.6f, dx2val_1: %12.6f \n", dx1pow_1, dx2val_1);
  if (bdebugmodel)  printf("--------------debug Tgov1Model::corrector, Tgov1 x: x1pow_1: %12.6f, x2val_1: %12.6f\n", x1pow_1, x2val_1);
  
  if (bdebugmodel)  printf("--------------debug Tgov1Model::corrector, Tgov1 Pmech = %12.6f\n", Pmech);
 
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power
 */
void gridpack::dynamic_simulation::Tgov1Model::setMechanicalPower(double pmech)
{
  Pmech = pmech; 
}

/**
 * Set the rotor speed deviation inside the governor
 * @param delta_o value of the rotor speed deviation
 */
void gridpack::dynamic_simulation::Tgov1Model::setRotorSpeedDeviation(double delta_w)
{
  delta_w = delta_w;
}

/** 
 * Get the value of the mechanical power
 * @return value of mechanical power
 */
double gridpack::dynamic_simulation::Tgov1Model::getMechanicalPower()
{
  return Pmech; 
}

/** 
 * Get the value of the rotor speed deviation
 * 
 */
/*double gridpack::dynamic_simulation::Tgov1Model::getRotorSpeedDeviation()
{
  return w;
}*/

/** 
 * Set the governor generator bus number
 */
 /**
void gridpack::dynamic_simulation::Tgov1Model::setExtBusNum(int ExtBusNum)
{
	p_bus_id = ExtBusNum;
}
*/	

/** 
 * Set the governor generator id
 */
 /**
void gridpack::dynamic_simulation::Tgov1Model::setExtGenId(std::string ExtGenId)
{
	p_ckt = ExtGenId;
}
*/
