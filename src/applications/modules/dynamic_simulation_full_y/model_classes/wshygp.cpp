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
 * @Latested modification with control blocks: Aug 23, 2023
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
#include "wshygp.hpp"

#define TS_THRESHOLD 4

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::WshygpModel::WshygpModel(void)
{
  //SecondGenExists = false;
  //OptionToModifyLimitsForInitialStateLimitViolation = true;
  w = 0.0;
  /*dx1Pmech = 0;
  dx2Td = 0;
  dx3Int = 0;
  dx4Der = 0;
  dx5Pelec = 0;
  dx6Valve = 0;
  dx7Gate = 0;
  dx1Pmech_1 = 0;
  dx2Td_1 = 0;
  dx3Int_1 = 0;
  dx4Der_1 = 0;
  dx5Pelec_1 = 0;
  dx6Valve_1 = 0;
  dx7Gate_1 = 0;*/
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

  if (!data->getValue(GENERATOR_MBASE, &GenMVABase, idx)) GenMVABase = 0.0; //printf ("Mbase = %8.4f \n", GenMVABase);

  Db1_blk.setparams(Db1, Err); 
  Filter_blk_d.setparams(1.0, Td);
  PIControl_blk.setparams(KP, KI);

  double a[2], b[2];
  a[0] = Tf; a[1] = 1.0;
  b[0] = KD; b[1] = 0.0;
  Feedback_blk_f.setcoeffs(a, b);

  Filter_blk_t.setparams(1.0, Tt);

  Filter_blk_p.setparams(KG, Tp, VELclose, VELopen, 0.0, 0.0);

  Integrator_blk.setparams(1.0);
  Db2_blk.setparams(Db2, Db2);  

  double uin[5], yin[5];
  uin[0] = Gv1; yin[0] = PGv1;
  uin[1] = Gv2; yin[0] = PGv2;
  uin[2] = Gv3; yin[0] = PGv3;
  uin[3] = Gv4; yin[0] = PGv4;
  uin[4] = Gv5; yin[0] = PGv5;
  NGV_blk.setparams(5, uin, yin);

  Leadlag_blk.setparams(Aturb*Tturb, Bturb*Tturb);

}

/**
 * Initialize governor model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::WshygpModel::init(double mag, double ang, double ts)
{
  double u1, u2, u3, u4, u5, u6, u7, u8, u9, u10;

  u10 = Pmech / (Trate / GenMVABase);

  u9 = Leadlag_blk.init_given_y(u10);

  u8 = NGV_blk.init_given_y(u9); 
  GV = u8;

  u7 = u8; // db2

  u6 = Integrator_blk.init_given_y(u7);

  u5 = Filter_blk_p.init_given_y(u6);

  u4 = Filter_blk_t.init_given_y(Pelec);

  u3 = u5 - PIControl_blk.init_given_y(u4) - Feedback_blk_f.init_given_y(u4);

  u2 = Filter_blk_d.init_given_y(u3);

  u1 = Pref - u2 - u4 / R; // or: u1 = Pref - u2 - u5 / R

  lastValue = u4 / R; // or: lastValue = u5 / R;

  w = u1;



  /*printf("wshygp: Pmech = %f\n", Pmech);
  // State 1
  double PGV = Pmech * GenMVABase / Trate;
  if (Bturb * Tturb < TS_THRESHOLD * ts) x1Pmech = 0;
  else x1Pmech = PGV * (1 - Aturb / Bturb);
  // State 7
  double GV = GainBlock.YtoX(PGV);
  // Initialize the BackLash
  BackLash.Initialize(Db2, GV);
  x7Gate = GV;
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (GV > Pmax) Pmax = GV;
    if (GV < Pmin) Pmin = GV;
  }
  // State 6
  x6Valve = 0;
  // State 5
  x5Pelec = GenPelec * GenMVABase / Trate; // TBD: where GenPelec is initilalized?
  // State 4
  x4Der = 0;
  // State 3
  double TempIn;
  if (KI == 0) {
    TempIn = GV / Kp;
    x3Int = 0;
  } else {
    TempIn = 0;
    x3Int = GV;
  }
  // State 2
  x2Td = TempIn;
  // Pref
  if (Tt == 0) {
    Pref = GV * R + TempIn;
  } else {
    Pref = x5Pelec + TempIn;
  }
  // Initialize the Intentional Deadband
  DBInt.Initialize(Db1, Err, w);
  printf("wshygp init: %f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1Pmech, x2Td, x3Int, x4Der, x5Pelec, x6Valve, x7Gate);*/
}

/**
 * computeModel - Updates the model states and gets the output
 */
void gridpack::dynamic_simulation::WshygpModel::computeModel(double t_inc,IntegrationStage int_flag)
{
  double u0, u1, u2, u3, u4, u5, u6, u7, u8, u9, u10;

  u0 = Db1_blk.getoutput(w);
  
  u1 = Filter_blk_d.getoutput(u0 + Pref - lastValue, t_inc, int_flag, true);

  u2 = PIControl_blk.getoutput(u1, t_inc, int_flag, true);

  u3 = Feedback_blk_f.getoutput(u1, t_inc, int_flag, true);

  u5 = u2 + u3; // or u5 = u2 + u3 + Filter_blk_t.getoutput(Pelec, t_inc, int_flag, true); 

  u5 = u5 - GV;

  u6 = Filter_blk_p.getoutput(u5, t_inc, int_flag, true);

  u7 = Integrator_blk.getoutput(u6, t_inc, int_flag, true);

  u8 = Db2_blk.getoutput(u7);

  u9 = NGV_blk.getoutput(u8);

  u10 = Leadlag_blk.getoutput(u9, t_inc, int_flag, true);

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
  /*if (!flag) {
    x1Pmech = x1Pmech_1;
    x2Td = x2Td_1;
    x3Int = x3Int_1;
    x4Der = x4Der_1;
    x5Pelec = x5Pelec_1;
    x6Valve = x6Valve_1;
    x7Gate = x7Gate_1;
  }

  // State 5
  if (Tt > TS_THRESHOLD * t_inc) {
    dx5Pelec = (GenPelec * GenMVABase /Trate - x5Pelec) / Tt;
  } else {
    x5Pelec = GenPelec * GenMVABase /Trate;
    dx5Pelec = 0;
  }
  // State 2
  double CV = x3Int + (x4Der + x2Td * KD / TF) + x2Td * Kp;
  double TempIn;
  if (Tt == 0) {
    x5Pelec = GenPelec * GenMVABase /Trate;
    TempIn = Pref - CV * R - DBInt.Output(w);
  } else {
    TempIn = Pref - x5Pelec  - DBInt.Output(w);
  }
  if (TD > TS_THRESHOLD * t_inc) {
    dx2Td = (TempIn - x2Td) / TD;
  } else {
    x2Td = TempIn;
    dx2Td = 0;
  }
  // State 3
  dx3Int = KI * TempIn;
  // State 4
  dx4Der = (-TempIn * KI / TF - x4Der) / TF;
  // State 6
  // Must recalculate CV because x2Td may have changed if TD = 0
  CV = x3Int + (x4Der + x2Td * KD / TF) + x2Td * Kp;
  if (x7Gate > Pmax) x7Gate = Pmax; // nonwindup
  if (x7Gate < Pmin) x7Gate = Pmin; // nonwindup
  double GV = BackLash.Output(x7Gate);
  //printf ("wshygp predictor GV = %8.4f, x1Pmech = %8.4f, x7Gate=%8.4f \n",GV, x1Pmech, x7Gate);
  TempIn = CV - GV;
  if (TP < TS_THRESHOLD * t_inc) x6Valve = TempIn * KG;
  if (x6Valve > VELopen) x6Valve = VELopen;
  if (x6Valve < VELclose) x6Valve = VELclose;
  if (TP < TS_THRESHOLD * t_inc) dx6Valve = 0;
  else dx6Valve = (TempIn * KG - x6Valve) / TP;
  if (dx6Valve > 0 && x6Valve >= VELopen) dx6Valve = 0; 
  if (dx6Valve < 0 && x6Valve <= VELclose) dx6Valve = 0;
  // State 7
  // Note: non windup stuff handled above
  //printf ("wshygp predictor x6Valve = %8.4f, \n",x6Valve);
  dx7Gate = x6Valve;
  if (dx7Gate > 0 && x7Gate >= Pmax) dx6Valve = 0;
  if (dx7Gate < 0 && x7Gate <= Pmin) dx6Valve = 0;
  // State 1
  double PGV = GainBlock.XtoY(GV);
  //printf ("wshygp predictor PGV = %8.4f, \n",PGV);
  if (Bturb * Tturb < TS_THRESHOLD * t_inc) dx1Pmech = 0;
  else dx1Pmech = (PGV * (1 - Aturb / Bturb) - x1Pmech) / (Bturb * Tturb);

  x1Pmech_1 = x1Pmech + dx1Pmech * t_inc;
  x2Td_1 = x2Td + dx2Td * t_inc;
  x3Int_1 = x3Int + dx3Int * t_inc;
  x4Der_1 = x4Der + dx4Der * t_inc;
  x5Pelec_1 = x5Pelec + dx5Pelec * t_inc;
  x6Valve_1 = x6Valve + dx6Valve * t_inc;
  x7Gate_1 = x7Gate + dx7Gate * t_inc;

  //printf("wshygp dx: %f\t%f\t%f\t%f\t%f\t%f\t%f\n", dx1Pmech, dx2Td, dx3Int, dx4Der, dx5Pelec, dx6Valve, dx7Gate);
  //printf("wshygp x: %f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1Pmech_1, x2Td_1, x3Int_1, x4Der_1, x5Pelec_1, x6Valve_1, x7Gate_1);

  // Note: you may want to cash the value PGV and keep it around
  GV = BackLash.Output(x7Gate);
  PGV = GainBlock.XtoY(GV);
  if (Bturb * Tturb < TS_THRESHOLD * t_inc) {
    Pmech = PGV * Trate / GenMVABase;
  } else {
    Pmech = (PGV * Aturb / Bturb + x1Pmech) * Trate / GenMVABase;
  } 
  
 // printf("wshygp Pmech = %f\n", Pmech);*/
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::WshygpModel::corrector(double t_inc, bool flag)
{
  computeModel(t_inc,CORRECTOR);
  /*// State 5
  if (Tt > TS_THRESHOLD * t_inc) {
    dx5Pelec_1 = (GenPelec * GenMVABase /Trate - x5Pelec_1) / Tt;
  } else {
    x5Pelec_1 = GenPelec * GenMVABase /Trate;
    dx5Pelec_1 = 0;
  }
  // State 2
  double CV = x3Int_1 + (x4Der_1 + x2Td_1 * KD / TF) + x2Td_1 * Kp;
  double TempIn;
  if (Tt == 0) {
    x5Pelec_1 = GenPelec * GenMVABase /Trate;
    TempIn = Pref - CV * R - DBInt.Output(w);
  } else {
    TempIn = Pref - x5Pelec_1  - DBInt.Output(w);
  }
  if (TD > TS_THRESHOLD * t_inc) {
    dx2Td_1 = (TempIn - x2Td_1) / TD;
  } else {
    x2Td_1 = TempIn;
    dx2Td_1 = 0;
  }
  // State 3
  dx3Int_1 = KI * TempIn;
  // State 4
  dx4Der_1 = (-TempIn * KI / TF - x4Der_1) / TF;
  // State 6
  // Must recalculate CV because x2Td may have changed if TD = 0
  CV = x3Int_1 + (x4Der_1 + x2Td_1 * KD / TF) + x2Td_1 * Kp;
  if (x7Gate > Pmax) x7Gate = Pmax; // nonwindup
  if (x7Gate < Pmin) x7Gate = Pmin; // nonwindup
  double GV = BackLash.Output(x7Gate);
  TempIn = CV - GV;
  if (TP < TS_THRESHOLD * t_inc) x6Valve = TempIn * KG;
  if (x6Valve > VELopen) x6Valve = VELopen;
  if (x6Valve < VELclose) x6Valve = VELclose;
  if (TP < TS_THRESHOLD * t_inc) dx6Valve = 0;
  else dx6Valve = (TempIn * KG - x6Valve) / TP;
  if (dx6Valve > 0 && x6Valve >= VELopen) dx6Valve = 0; 
  if (dx6Valve < 0 && x6Valve <= VELclose) dx6Valve = 0;
  // State 7
  // Note: non windup stuff handled above
  dx7Gate = x6Valve;
  if (dx7Gate > 0 && x7Gate >= Pmax) dx6Valve = 0;
  if (dx7Gate < 0 && x7Gate <= Pmin) dx6Valve = 0;
  // State 1
  double PGV = GainBlock.XtoY(GV);
  if (Bturb * Tturb < TS_THRESHOLD * t_inc) dx1Pmech_1 = 0;
  else dx1Pmech_1 = (PGV * (1 - Aturb / Bturb) - x1Pmech_1) / (Bturb * Tturb);

  x1Pmech_1 = x1Pmech + (dx1Pmech + dx1Pmech_1) / 2.0 * t_inc;
  x2Td_1 = x2Td + (dx2Td + dx2Td_1) / 2.0 * t_inc;
  x3Int_1 = x3Int + (dx3Int + dx3Int_1) / 2.0 * t_inc;
  x4Der_1 = x4Der + (dx4Der + dx4Der_1) / 2.0 * t_inc;
  x5Pelec_1 = x5Pelec + (dx5Pelec + dx5Pelec_1) / 2.0 * t_inc;
  x6Valve_1 = x6Valve + (dx6Valve + dx6Valve_1) / 2.0 * t_inc;
  x7Gate_1 = x7Gate + (dx7Gate + dx7Gate_1) / 2.0 * t_inc;
 
 // printf("wshygp dx: %f\t%f\t%f\t%f\t%f\t%f\t%f\n", dx1Pmech, dx2Td, dx3Int, dx4Der, dx5Pelec, dx6Valve, dx7Gate);
  //printf("wshygp x: %f\t%f\t%f\t%f\t%f\t%f\t%f\n", x1Pmech_1, x2Td_1, x3Int_1, x4Der_1, x5Pelec_1, x6Valve_1, x7Gate_1);

  // Note: you may want to cash the value PGV and keep it around
  GV = BackLash.Output(x7Gate_1);
  PGV = GainBlock.XtoY(GV);
  if (Bturb * Tturb < TS_THRESHOLD * t_inc) {
    Pmech = PGV * Trate / GenMVABase;
  } else {
    Pmech = (PGV * Aturb / Bturb + x1Pmech) * Trate / GenMVABase;
  } 
  
 // printf("wshygp Pmech = %f\n", Pmech);*/
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power
 */
void gridpack::dynamic_simulation::WshygpModel::setMechanicalPower(double pmech)
{
  Pmech = pmech; 
}

void gridpack::dynamic_simulation::WshygpModel::setPelec(double pelec)
{
  Pelec = pelec; 
}

void gridpack::dynamic_simulation::WshygpModel::setPref(double pref)
{
  Pref = pref; 
}

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
