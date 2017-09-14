/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   wsieg1.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 11, 2015
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_governor_model.hpp"
#include "wsieg1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Wsieg1Model::Wsieg1Model(void)
{
  SecondGenExists = false;
  OptionToModifyLimitsForInitialStateLimitViolation = false;
  w = 0.0;
  dx1LL = 0;
  dx2GovOut = 0;
  dx3Turb1 = 0;
  dx4Turb2 = 0;
  dx5Turb3 = 0;
  dx6Turb4 = 0;
  dx1LL_1 = 0;
  dx2GovOut_1 = 0;
  dx3Turb1_1 = 0;
  dx4Turb2_1 = 0;
  dx5Turb3_1 = 0;
  dx6Turb4_1 = 0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Wsieg1Model::~Wsieg1Model(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to
 * Wsieg1Model
 */
void gridpack::dynamic_simulation::Wsieg1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  //if (!data->getValue(GOVERNOR_JBUS, &JBUS, idx)) JBUS = 0.0; // JBUS 
  //if (!data->getValue(GOVERNOR_M, &M, idx)) M = 0.0; // M
  if (!data->getValue(GOVERNOR_K, &K, idx)) K = 0.0; // K
  if (!data->getValue(GOVERNOR_T1, &T1, idx)) T1 = 0.0; // T1
  if (!data->getValue(GOVERNOR_T2, &T2, idx)) T2 = 0.0; // T2
  if (!data->getValue(GOVERNOR_T3, &T3, idx)) T3 = 0.0; // T3
  if (!data->getValue(GOVERNOR_UO, &Uo, idx)) Uo = 0.0; // Uo
  if (!data->getValue(GOVERNOR_UC, &Uc, idx)) Uc = 0.0; // Uc
  if (!data->getValue(GOVERNOR_PMAX, &Pmax, idx)) Pmax = 0.0; // Pmax
  if (!data->getValue(GOVERNOR_PMIN, &Pmin, idx)) Pmin = 0.0; // Pmin
  if (!data->getValue(GOVERNOR_T4, &T4, idx)) T4 = 0.0; // T4
  if (!data->getValue(GOVERNOR_K1, &K1, idx)) K1 = 0.0; // K1
  if (!data->getValue(GOVERNOR_K2, &K2, idx)) K2 = 0.0; // K2
  if (!data->getValue(GOVERNOR_T5, &T5, idx)) T5 = 0.0; // T5
  if (!data->getValue(GOVERNOR_K3, &K3, idx)) K3 = 0.0; // K3
  if (!data->getValue(GOVERNOR_K4, &K4, idx)) K4 = 0.0; // K4
  if (!data->getValue(GOVERNOR_T6, &T6, idx)) T6 = 0.0; // T6
  if (!data->getValue(GOVERNOR_K5, &K5, idx)) K5 = 0.0; // K5
  if (!data->getValue(GOVERNOR_K6, &K6, idx)) K6 = 0.0; // K6
  if (!data->getValue(GOVERNOR_T7, &T7, idx)) T7 = 0.0; // T7
  if (!data->getValue(GOVERNOR_K7, &K7, idx)) K7 = 0.0; // K7
  if (!data->getValue(GOVERNOR_K8, &K8, idx)) K8 = 0.0; // K8
  if (!data->getValue(GOVERNOR_DB1, &Db1, idx)) Db1 = 0.0; // Db1
  if (!data->getValue(GOVERNOR_ERR, &Err, idx)) Err = 0.0; // Err
  if (!data->getValue(GOVERNOR_DB2, &Db2, idx)) Db2 = 0.0; // Db2
  /*if (!data->getValue(GOVERNOR_GV1, &Gv1, idx)) Gv1 = 0.0; // Gv1
  if (!data->getValue(GOVERNOR_PGV1, &PGv1, idx)) PGv1 = 0.0; // PGv1
  if (!data->getValue(GOVERNOR_GV2, &Gv2, idx)) Gv2 = 0.0; // Gv2
  if (!data->getValue(GOVERNOR_PGV2, &PGv2, idx)) PGv2 = 0.0; // PGv2
  if (!data->getValue(GOVERNOR_GV3, &Gv3, idx)) Gv3 = 0.0; // Gv3
  if (!data->getValue(GOVERNOR_PGV3, &PGv3, idx)) PGv3 = 0.0; // PGv3
  if (!data->getValue(GOVERNOR_GV4, &Gv4, idx)) Gv4 = 0.0; // Gv4
  if (!data->getValue(GOVERNOR_PGV4, &PGv4, idx)) PGv4 = 0.0; // PGv4
  if (!data->getValue(GOVERNOR_GV5, &Gv5, idx)) Gv5 = 0.0; // Gv5
  if (!data->getValue(GOVERNOR_PGV5, &PGv5, idx)) PGv5 = 0.0; // PGv5*/
  if (!data->getValue(GOVERNOR_IBLOCK, &Iblock, idx)) Iblock = 0.0; // Iblock
}

/**
 * Initialize governor model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Wsieg1Model::init(double mag, double ang, double ts)
{
  ///printf("wsieg1: Pmech1 = %f, Pmech2 = %f\n", Pmech1, Pmech2);
  double PGV;
  if (K1 + K3 + K5 + K7 > 0) 
    PGV = Pmech1 / (K1 + K3 + K5 + K7);
  else if (K2 + K4 + K6 + K8 > 0) 
    PGV = Pmech2 / (K2 + K4 + K6 + K8);
  else 
    PGV = 0;
  if (SecondGenExists && (Pmech2 != 0) && (K2 + K4 + K6 + K8 > 0) && (PGV != 0)) {
    double temp = Pmech2 / PGV * (K2 + K4 + K6 + K8);
    K2 = temp * K2;
    K4 = temp * K4;
    K6 = temp * K6;
    K8 = temp * K8;
  }
  x6Turb4 = PGV;
  x5Turb3 = PGV;
  x4Turb2 = PGV;
  x3Turb1 = PGV;
  double GV = GainBlock.YtoX(PGV); // TBD: check GainBlock?
  //printf("GV = %f\n", GV);
  x2GovOut = GV;
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (GV > Pmax) Pmax = GV;
    if (GV < Pmin) Pmin = GV;
  }
  Pref = GV;
  // Initialize the Backlash
  BackLash.Initialize(Db2, GV);
  // Initialize the Intentional Deadband
  DBInt.Initialize(Db1, Err, w); // TBD: has w been set at gensal init step? yes
  // Note: (GV > Pmax) or (GV < Pmin) is an initial state violation
  if (Iblock == 1 && Pmin == 0) Pmin = GV;
  if (Iblock == 2 && Pmax == 0) Pmax = GV;
  if (Iblock == 3 && Pmin == 0) Pmin = GV;
  if (Iblock == 3 && Pmax == 0) Pmax = GV;
  if (T1 > 4 * ts) x1LL = GV * (1 - T2 / T1);
  else x1LL = GV;
  //printf("T1 = %f, T2 = %f, ts = %f\n", T1, T2, ts);
  ///printf("wsieg1 init: %f\t%f\t%f\t%f\t%f\t%f\n", x1LL, x2GovOut, x3Turb1, x4Turb2, x5Turb3, x6Turb4);
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wsieg1Model::predictor(double t_inc, bool flag)
{
  if (!flag) {
    x1LL = x1LL_1;
    x2GovOut = x2GovOut_1;
    x3Turb1 = x3Turb1_1;
    x4Turb2 = x4Turb2_1;
    x5Turb3 = x5Turb3_1;
    x6Turb4 = x6Turb4_1;
  }
  // State 1
  //printf("w = %f\n", w);
  double TempIn1 = K * w;//DBInt.Output(w);
  double TempOut;
  if (T1 > 4 * t_inc) {
    dx1LL = (TempIn1 * ( 1 - T2 / T1) - x1LL) / T1;
    TempOut = TempIn1 * (T2 / T1) + x1LL;
  } else 
    TempOut = TempIn1;
  //printf("T1 = %f, T2 = %f, x1LL = %f, K = %f, TempIn1 = %f\n", T1, T2, x1LL, K, TempIn1);
  // State 2
  // enforce non-windup limits
  double TempIn2;
  if (x2GovOut > Pmax) x2GovOut = Pmax;
  else if (x2GovOut < Pmin) x2GovOut = Pmin;
  double GV = BackLash.Output(x2GovOut);
  if (T3 < 4 * t_inc) TempIn2 = (+ Pref - TempOut - GV) / (4 * t_inc);
  else TempIn2  = (+ Pref - TempOut - GV) / T3;
  if (TempIn2 > Uo) TempIn2 = Uo;
  else if (TempIn2 < Uc) TempIn2 = Uc;
  dx2GovOut = TempIn2;
  //printf("TempIn1 = %f, TempOut = %f, w = %f, TempIn2 = %f\n", TempIn1, TempOut, w, TempIn2);
  // enforce non-windup limits
  if (dx2GovOut > 0 && x2GovOut >= Pmax) dx2GovOut = 0;
  else if (dx2GovOut <0 && x2GovOut <= Pmin) dx2GovOut = 0;
  // State 3
  double PGV = GainBlock.XtoY(GV);
  if (T4 < 4 * t_inc) {
    x3Turb1 = PGV;
    dx3Turb1 = 0;
  } else
    dx3Turb1 = (PGV - x3Turb1) / T4;
  // State 4
  if (T5 < 4 * t_inc) {
    x4Turb2 = x3Turb1;
    dx4Turb2 = 0;
  } else
    dx4Turb2 = (x3Turb1 - x4Turb2) / T5;
  // State 5
  if (T6 < 4 * t_inc) {
    x5Turb3 = x4Turb2;
    dx5Turb3 = 0;
  } else
    dx5Turb3 = (x4Turb2 - x5Turb3) / T6;
  // State 6
  if (T7 < 4 * t_inc) {
    x6Turb4 = x5Turb3;
    dx6Turb4 = 0;
  } else
    dx6Turb4 = (x5Turb3 - x6Turb4) / T7;

  x1LL_1 = x1LL + dx1LL * t_inc;
  x2GovOut_1 = x2GovOut + dx2GovOut * t_inc;
  x3Turb1_1 = x3Turb1 + dx3Turb1 * t_inc;
  x4Turb2_1 = x4Turb2 + dx4Turb2 * t_inc;
  x5Turb3_1 = x5Turb3 + dx5Turb3 * t_inc;
  x6Turb4_1 = x6Turb4 + dx6Turb4 * t_inc;

  ///printf("wsieg1 dx: %f\t%f\t%f\t%f\t%f\t%f\n", dx1LL, dx2GovOut, dx3Turb1, dx4Turb2, dx5Turb3, dx6Turb4);
  ///printf("wsieg1 x: %f\t%f\t%f\t%f\t%f\t%f\n", x1LL_1, x2GovOut_1, x3Turb1_1, x4Turb2_1, x5Turb3_1, x6Turb4_1);

  Pmech1 = x3Turb1_1 * K1 + x4Turb2_1 * K3 + x5Turb3_1 * K5 + x6Turb4_1 * K7;
  Pmech2 = x3Turb1_1 * K2 + x4Turb2_1 * K4 + x5Turb3_1 * K6 + x6Turb4_1 * K8;
  
  ///printf("wsieg1 Pmech1 = %f, Pmech2 = %f\n", Pmech1, Pmech2);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wsieg1Model::corrector(double t_inc, bool flag)
{
  // State 1
  double TempIn1 = K * w;//DBInt.Output(w);
  double TempOut;
  if (T1 > 4 * t_inc) {
    dx1LL_1 = (TempIn1 * ( 1 - T2 / T1) - x1LL_1) / T1;
    TempOut = TempIn1 * (T2 / T1) + x1LL_1;
  } else 
    TempOut = TempIn1;
  // State 2
  // enforce non-windup limits
  double TempIn2;
  if (x2GovOut_1 > Pmax) x2GovOut_1 = Pmax;
  else if (x2GovOut_1 < Pmin) x2GovOut_1 = Pmin;
  double GV = BackLash.Output(x2GovOut_1);
  if (T3 < 4 * t_inc) TempIn2 = (+ Pref - TempOut - GV) / (4 * t_inc);
  else TempIn2  = (+ Pref - TempOut - GV) / T3;
  if (TempIn2 > Uo) TempIn2 = Uo;
  else if (TempIn2 < Uc) TempIn2 = Uc;
  dx2GovOut_1 = TempIn2;
  // enforce non-windup limits
  if (dx2GovOut_1 > 0 && x2GovOut_1 >= Pmax) dx2GovOut_1 = 0;
  else if (dx2GovOut_1 <0 && x2GovOut_1 <= Pmin) dx2GovOut_1 = 0;
  // State 3
  double PGV = GainBlock.XtoY(GV);
  if (T4 < 4 * t_inc) {
    x3Turb1_1 = PGV;
    dx3Turb1_1 = 0;
  } else
    dx3Turb1_1 = (PGV - x3Turb1_1) / T4;
  // State 4
  if (T5 < 4 * t_inc) {
    x4Turb2_1 = x3Turb1_1;
    dx4Turb2_1 = 0;
  } else
    dx4Turb2_1 = (x3Turb1_1 - x4Turb2_1) / T5;
  // State 5
  if (T6 < 4 * t_inc) {
    x5Turb3_1 = x4Turb2_1;
    dx5Turb3_1 = 0;
  } else
    dx5Turb3_1 = (x4Turb2_1 - x5Turb3_1) / T6;
  // State 6
  if (T7 < 4 * t_inc) {
    x6Turb4_1 = x5Turb3_1;
    dx6Turb4_1 = 0;
  } else
    dx6Turb4_1 = (x5Turb3_1 - x6Turb4_1) / T7;

  x1LL_1 = x1LL + (dx1LL + dx1LL_1) / 2.0 * t_inc;
  x2GovOut_1 = x2GovOut + (dx2GovOut + dx2GovOut_1) / 2.0 * t_inc;
  x3Turb1_1 = x3Turb1 + (dx3Turb1 + dx3Turb1_1) / 2.0 * t_inc;
  x4Turb2_1 = x4Turb2 + (dx4Turb2 + dx4Turb2_1) / 2.0 * t_inc;
  x5Turb3_1 = x5Turb3 + (dx5Turb3 + dx5Turb3_1) / 2.0 * t_inc;
  x6Turb4_1 = x6Turb4 + (dx6Turb4 + dx6Turb4_1) / 2.0 * t_inc;
 
  ///printf("wsieg1 dx: %f\t%f\t%f\t%f\t%f\t%f\n", dx1LL_1, dx2GovOut_1, dx3Turb1_1, dx4Turb2_1, dx5Turb3_1, dx6Turb4_1);
  ///printf("wsieg1 x: %f\t%f\t%f\t%f\t%f\t%f\n", x1LL_1, x2GovOut_1, x3Turb1_1, x4Turb2_1, x5Turb3_1, x6Turb4_1);

  Pmech1 = x3Turb1_1 * K1 + x4Turb2_1 * K3 + x5Turb3_1 * K5 + x6Turb4_1 * K7;
  Pmech2 = x3Turb1_1 * K2 + x4Turb2_1 * K4 + x5Turb3_1 * K6 + x6Turb4_1 * K8;

  ///printf("wsieg1 Pmech1 = %f, Pmech2 = %f\n", Pmech1, Pmech2);
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power
 */
void gridpack::dynamic_simulation::Wsieg1Model::setMechanicalPower(double pmech)
{
  Pmech1 = pmech; 
  Pmech2 = pmech; 
}

/**
 * Set the rotor speed deviation inside the governor
 * @param delta_o value of the rotor speed deviation
 */
void gridpack::dynamic_simulation::Wsieg1Model::setRotorSpeedDeviation(double delta_o)
{
  w = delta_o;
}

/** 
 * Get the value of the mechanical power
 * @return value of mechanical power
 */
double gridpack::dynamic_simulation::Wsieg1Model::getMechanicalPower()
{
  return Pmech1; 
}

/** 
 * Get the value of the rotor speed deviation
 * 
 */
/*double gridpack::dynamic_simulation::Wsieg1Model::getRotorSpeedDeviation()
{
  return w;
}*/
