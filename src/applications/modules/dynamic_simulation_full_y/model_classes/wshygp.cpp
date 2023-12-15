/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   wshygp.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 11, 2015
 * @Latested modification with control blocks: Aug 28, 2023
 * @Last modified by Yuan: November 21, 2023
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_governor_model.hpp"
#include "wshygp.hpp"

#define TS_THRESHOLD 4
#define KI_THRESHOLD 1e-6

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::WshygpModel::WshygpModel(void)
{
  w = 0.0;
  GenMVABase = 100.0; 
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::WshygpModel::~WshygpModel(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to
 * WshygpModel
 */
void gridpack::dynamic_simulation::WshygpModel::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(GOVERNOR_DB1, &Db1, idx)) Db1 = 0.0; //printf ("Db1 = %8.4f \n", Db1);
  if (!data->getValue(GOVERNOR_ERR, &Err, idx)) Err = 0.0; //printf ("Err = %8.4f \n", Err);
  if (!data->getValue(GOVERNOR_TD, &Td, idx)) Td = 0.0; //printf ("TD = %8.4f \n", TD);
  if (!data->getValue(GOVERNOR_KI, &KI, idx)) KI = 0.0; //printf ("KI = %8.4f \n", KI);
  if (!data->getValue(GOVERNOR_TF, &Tf, idx)) Tf = 0.0; //printf ("TF = %8.4f \n", TF);
  if (!data->getValue(GOVERNOR_KD, &KD, idx)) KD = 0.0; //printf ("KD = %8.4f \n", KD);
  if (!data->getValue(GOVERNOR_KP, &KP, idx)) KP = 0.0; //printf ("Kp = %8.4f \n", Kp);
  if (!data->getValue(GOVERNOR_R, &R, idx)) R = 0.0;    //printf ("R = %8.4f \n", R);
  if (!data->getValue(GOVERNOR_TT, &Tt, idx)) Tt = 0.0; //printf ("Tt = %8.4f \n", Tt);
  if (!data->getValue(GOVERNOR_KG, &KG, idx)) KG = 0.0; //printf ("KG = %8.4f \n", KG);
  if (!data->getValue(GOVERNOR_TP, &Tp, idx)) Tp = 0.0; //printf ("TP = %8.4f \n", TP);
  if (!data->getValue(GOVERNOR_VELOPEN, &VELopen, idx)) VELopen = 0.0; //printf ("VELopen = %8.4f \n", VELopen);
  if (!data->getValue(GOVERNOR_VELCLOSE, &VELclose, idx)) VELclose = 0.0; //printf ("VELclose = %8.4f \n", VELclose);
  if (!data->getValue(GOVERNOR_PMAX, &Pmax, idx)) Pmax = 0.0; //printf ("Pmax = %8.4f \n", Pmax);
  if (!data->getValue(GOVERNOR_PMIN, &Pmin, idx)) Pmin = 0.0; //printf ("Pmin = %8.4f \n", Pmin);
  if (!data->getValue(GOVERNOR_DB2, &Db2, idx)) Db2 = 0.0; //printf ("Db2 = %8.4f \n", Db2);
  if (!data->getValue(GOVERNOR_GV1, &Gv1, idx)) Gv1 = 0.0; // Gv1
  if (!data->getValue(GOVERNOR_PGV1, &PGv1, idx)) PGv1 = 0.0; // PGv1
  if (!data->getValue(GOVERNOR_GV2, &Gv2, idx)) Gv2 = 0.0; // Gv2
  if (!data->getValue(GOVERNOR_PGV2, &PGv2, idx)) PGv2 = 0.0; // PGv2
  if (!data->getValue(GOVERNOR_GV3, &Gv3, idx)) Gv3 = 0.0; // Gv3
  if (!data->getValue(GOVERNOR_PGV3, &PGv3, idx)) PGv3 = 0.0; // PGv3
  if (!data->getValue(GOVERNOR_GV4, &Gv4, idx)) Gv4 = 0.0; // Gv4
  if (!data->getValue(GOVERNOR_PGV4, &PGv4, idx)) PGv4 = 0.0; // PGv4
  if (!data->getValue(GOVERNOR_GV5, &Gv5, idx)) Gv5 = 0.0; // Gv5
  if (!data->getValue(GOVERNOR_PGV5, &PGv5, idx)) PGv5 = 0.0; // PGv5
  if (!data->getValue(GOVERNOR_ATURB, &Aturb, idx)) Aturb = 0.0; //printf ("Aturb = %8.4f \n", Aturb);
  if (!data->getValue(GOVERNOR_BTURB, &Bturb, idx)) Bturb = 0.0; //printf ("Bturb = %8.4f \n", Bturb);
  if (!data->getValue(GOVERNOR_TTURB, &Tturb, idx)) Tturb = 0.0; //printf ("Tturb = %8.4f \n", Tturb);
  if (!data->getValue(GOVERNOR_TRATE, &Trate, idx)) Trate = 0.0; //printf ("Trate = %8.4f \n", Trate);
  

  // printf("Db1=%f,Err=%f,Td=%f,KI=%f,Tf=%f,KD=%f,KP=%f,R=%f,Tt=%f,KG=%f,Tp=%f,VELopen=%f,VELclose=%f,Pmax=%f,Pmin=%f,Db2=%f\n",
  // Db1,Err,Td,KI,Tf,KD,KP,R,Tt,KG,Tp,VELopen,VELclose,Pmax,Pmin,Db2);
  // printf("Gv1=%f,PGv1=%f,Gv2=%f,PGv2=%f,Gv3=%f,PGv3=%f,Gv4=%f,PGv4=%f,Gv5=%f,PGv5=%f,Aturb=%f,Bturb=%f,Tturb=%f,Trate=%f\n",
  // Gv1,PGv1,Gv2,PGv2,Gv3,PGv3,Gv4,PGv4,Gv5,PGv5,Aturb,Bturb,Tturb,Trate);

  
  zero_TD = false;
  zero_KI = false;
  zero_TF = false;
  zero_TT = false;
  zero_TP = false;
  zero_TTURB_BTURB = false;
  Tt = 0.0; // hard coded because there is no setPelec method in this and base_governor class
  Pelec = 0.0; // the same reason
  
  OptionToModifyLimitsForInitialStateLimitViolation = true;
}

/**
 * Initialize governor model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::WshygpModel::init(double mag, double ang, double ts)
{
  // parameter sanity checks
  if (Td < TS_THRESHOLD * ts) zero_TD = true;
  if (KI < KI_THRESHOLD) zero_KI = true;
  if (Tf < TS_THRESHOLD * ts) zero_TF = true;
  if (Tt < TS_THRESHOLD * ts) zero_TT = true;
  if (Tp < TS_THRESHOLD * ts) zero_TP = true;
  if (abs(Tturb * Bturb) < TS_THRESHOLD * ts) zero_TTURB_BTURB = true;
  
  if (zero_TD) printf("Warning: Td=%f is less than %d times timestep=%f.\n", Td, TS_THRESHOLD, ts);
  if (zero_KI) printf("Warning: KI=%f is less than %f.\n", KI, KI_THRESHOLD);
  if (zero_TF) printf("Warning: Tf=%f is less than %d times timestep=%f.\n", Tf, TS_THRESHOLD, ts);
  if (zero_TT) printf("Warning: Tt=%f is less than %d times timestep=%f.\n", Tt, TS_THRESHOLD, ts);
  if (zero_TP) printf("Warning: Tp=%f is less than %d times timestep=%f.\n", Tp, TS_THRESHOLD, ts);
  if (zero_TTURB_BTURB) printf("Warning: Tturb * Bturb=%f is less than %d times timestep=%f.\n", Tturb * Bturb, TS_THRESHOLD, ts);
  
  if (KP == 0.0 && KI < KI_THRESHOLD) {
	  KP = 1.0;
	  printf("KP and KI cannot be both zeros; reset KP=1.0.\n");
  }
  
  if (Tf < TS_THRESHOLD * ts) {
	  Tf = TS_THRESHOLD * ts;
	  zero_TF = false;
	  printf("Tf is reset to be equal to %d times timestep=%f.\n", TS_THRESHOLD, ts);
  }
  if (KD < TS_THRESHOLD * ts) {
	  KD = TS_THRESHOLD * ts;
	  printf("KD is reset to be equal to %d times timestep=%f.\n", TS_THRESHOLD, ts);
  }
  
  if (zero_TTURB_BTURB) printf("Since Tturb * Bturb is too small, this block will be treated as a gain=1.0 block.\n");

  // set model parameters
  Db1_blk.setparams(Db1, Err);
  if (!zero_TD) Filter_blk_d.setparams(1.0, Td);  // else {a gain=1 block}
  if (!zero_KI) PIControl_blk.setparams(KP, KI); // else {a gain=KP block}

  double a[2], b[2];
  a[0] = Tf; a[1] = 1.0;
  b[0] = KD; b[1] = 0.0;
  Feedback_blk_f.setcoeffs(a, b);

  if (!zero_TT) Filter_blk_t.setparams(1.0, Tt); // else {a gain=1 block}
  if (!zero_TP) Filter_blk_p.setparams(KG, Tp, VELclose, VELopen, -1000.0, 1000.0); // else {a gain=KG block}

  Integrator_blk.setparams(1.0);
  Integrator_blk.setylimits(Pmin, Pmax);
  
  Db2_blk.setparams(Db2, 2);  

  double uin[5], yin[5];
  uin[0] = Gv1; yin[0] = PGv1;
  uin[1] = Gv2; yin[0] = PGv2;
  uin[2] = Gv3; yin[0] = PGv3;
  uin[3] = Gv4; yin[0] = PGv4;
  uin[4] = Gv5; yin[0] = PGv5;
  NGV_blk.setparams(5, uin, yin);

  if (!zero_TTURB_BTURB) Leadlag_blk.setparams(Aturb*Tturb, Bturb*Tturb); // else {Aturb*Tturb=0, a gain=1 block}
  
  // block initialization
  double u1, u2, u3, u4, u5, u6, u7, u8, u9, u10;

  u10 = Pmech / (Trate / GenMVABase);

  u9 = Leadlag_blk.init_given_y(u10);

  // u8 = NGV_blk.init_given_y(u9); 
  u8 = u9; // treat NGV_blk as gain=1
  
  GV = u8;

  u7 = u8; // db2
  Db2_blk.init_given_y(GV);

  u6 = Integrator_blk.init_given_y(u7);
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
	  if (u7 > Pmax) Pmax = u7 + 0.1;
	  if (u7 < Pmin) Pmin = u7 - 0.1;
  }

  if (!zero_TP) u5 = Filter_blk_p.init_given_y(u6);
  else u5 = u6/KG; // just a gain=KG block
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
	  if (!zero_TP) {
		  if (u6 > VELopen) VELopen = u6 + 0.1;
		  if (u6 < VELclose) VELclose = u6 - 0.1;
	  }
  }
  
  CV = u5 + GV;
  
  // in this model, the input and output of Feedback_blk_f should be both zero
  // we can first initialize PIControl_blk using output, then the input u will be zero
  // we then use input u to initialize Feedback_blk_f
  u3 = PIControl_blk.init_given_y(CV);  // u5 should be zero
  double temp_y;
  temp_y = Feedback_blk_f.init_given_u(u3); // temp_y should also be zero
  
  if (!zero_TD) u2 = Filter_blk_d.init_given_y(u3);
  else u2 = u3; // just a gain=1 block

  if (Tt > 0.0) u4 = Filter_blk_t.init_given_u(Pelec);
  else u4 = 0.0;

  u1 = w;
  Db1_blk.init_given_u(w);

  if (Tt == 0.0) Pref = u2  + u1 + CV * R;
  else if (Tt > 0) Pref = u2 + u1 + u4 * R;
  
}

/**
 * computeModel - Updates the model states and gets the output
 */
void gridpack::dynamic_simulation::WshygpModel::computeModel(double t_inc,IntegrationStage int_flag)
{
  double u0, u1, u2, u3, u4, u5, u6, u7, u8, u9, u10;

  u0 = Db1_blk.getoutput(w);
  
  if (Tt == 0) {
    if (!zero_TD) u1 = Filter_blk_d.getoutput(Pref - u0 - CV * R, t_inc, int_flag, true);
    else u1 = Pref - u0 - CV * R; // a gain=1 block 
  } else if (Tt > 0) {
    if (!zero_TT) u4 = Filter_blk_t.getoutput(Pelec);
	else u4 = Pelec; // a gain=1 block
	
    if (!zero_TD) u1 = Filter_blk_d.getoutput(Pref - u0 - u4 * R, t_inc, int_flag, true);
	else u1 = Pref - u0 - u4 * R; // a gain=1 block 
  }
  
  if (!zero_KI) u2 = PIControl_blk.getoutput(u1, t_inc, int_flag, true);
  else u2 = u1 * KP;  // a gain=KP block

  u3 = Feedback_blk_f.getoutput(u1, t_inc, int_flag, true);

  u5 = u2 + u3;
  CV = u5;

  u5 = CV - GV;

  if (!zero_TP) u6 = Filter_blk_p.getoutput(u5, t_inc, int_flag, true);
  else u6 = u5 * KG;  // a gain=KG block

  u7 = Integrator_blk.getoutput(u6, t_inc, int_flag, true);

  u8 = Db2_blk.getoutput(u7);

  // u9 = NGV_blk.getoutput(u8);
  u9 = u8; // currently treat the NGV block as a gain=1 block

  if (!zero_TTURB_BTURB) u10 = Leadlag_blk.getoutput(u9, t_inc, int_flag, true);
  else u10 = u9;  // a gain=1 block

  Pmech = u10 * Trate / GenMVABase;

}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::WshygpModel::predictor(double t_inc, bool flag)
{
  computeModel(t_inc,PREDICTOR);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::WshygpModel::corrector(double t_inc, bool flag)
{
  computeModel(t_inc,CORRECTOR);
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power
 */
void gridpack::dynamic_simulation::WshygpModel::setMechanicalPower(double pmech)
{
  Pmech = pmech; 
}

// void gridpack::dynamic_simulation::WshygpModel::setPelec(double pelec)
// {
  // Pelec = pelec; 
// }

// void gridpack::dynamic_simulation::WshygpModel::setPref(double pref)
// {
  // Pref = pref; 
// }

/**
 * Set the rotor speed deviation inside the governor
 * @param delta_o value of the rotor speed deviation
 */
void gridpack::dynamic_simulation::WshygpModel::setRotorSpeedDeviation(double delta_o)
{
  w = delta_o;
}

/** 
 * Get the value of the mechanical power
 * @return value of mechanical power
 */
double gridpack::dynamic_simulation::WshygpModel::getMechanicalPower()
{
  return Pmech; 
}

/**
 * Set internal state parameter in governor
 * @param name character string corresponding to state variable
 * @param value new value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::WshygpModel::setState(std::string name,
    double value)
{
  return false;
}

/**
 * Get internal state parameter in governor
 * @param name character string corresponding to state variable
 * @param value current value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::WshygpModel::getState(std::string name,
    double *value)
{
  return false;
}

/** 
 * Get the value of the rotor speed deviation
 * 
 */
/*double gridpack::dynamic_simulation::WshygpModel::getRotorSpeedDeviation()
{
  return w;
}*/
