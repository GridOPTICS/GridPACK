/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wsieg1.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   10/02/19
 *  
 * @brief  
 *
 *
 */

#include <wsieg1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

#define TS_THRESHOLD 1

Wsieg1Gov::Wsieg1Gov(void)
{
  x1LL = 0.0; 
  x2GovOut = 0.0; 
  x3Turb1 = 0.0; 
  x4Turb2 = 0.0; 
  x5Turb3 = 0.0; 
  x6Turb4 = 0.0;
  dx1LL = 0.0;
  dx2GovOut = 0.0;
  dx3Turb1 = 0.0;
  dx4Turb2 = 0.0;
  dx5Turb3 = 0.0;
  dx6Turb4 = 0.0;
  K = 0.0;
  T1 = 0.0;
  T2 = 0.0;
  T3 = 0.0;
  Uo = 0.0;
  Uc = 0.0;
  Pmax = 0.0;
  Pmin = 0.0;
  T4 = 0.0;
  K1 = 0.0;
  K2 = 0.0;
  T5 = 0.0;
  K3 = 0.0;
  K4 = 0.0;
  T6 = 0.0;
  K5 = 0.0;
  K6 = 0.0;
  T7 = 0.0;
  K7 = 0.0;
  K8 = 0.0;
  SecondGenExists = false;
  OptionToModifyLimitsForInitialStateLimitViolation = false;
}

Wsieg1Gov::~Wsieg1Gov(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to BaseGovernorModel
 */
void Wsieg1Gov::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseGovModel::load(data,idx); // load parameters in base governor model
  
  // load parameters for the model type
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
  if (!data->getValue(GOVERNOR_IBLOCK, &Iblock, idx)) Iblock = 0.0; // Iblock

  //printf("K=%f,T1=%f,T2=%f,T3=%f,Uo=%f,Uc=%f,Pmax=%f,Pmin=%f,T4=%f,K1=%f,K2=%f,T5=%f,K3=%f,K4=%f,T6=%f,K5=%f,K6=%f,T7=%f,K7=%f,K8=%f,Db1=%f,Err=%f,Db2=%f,Iblock=%f\n", K, T1, T2, T3, Uo, Uc, Pmax, Pmin, T4, K1, K2, T5, K3, K4, T6, K5, K6, T7, K7, K8, Db1, Err, Db2, Iblock);

  // Convert governor parameters from machine base to MVA base
  /*p_H *= mbase/sbase;
  p_D *= mbase/sbase;
  p_Xdp /= mbase/sbase;*/

}

/**
 * Initialize governor model before calculation
 * @param [output] values - array where initialized governor variables should be set
 */
void Wsieg1Gov::init(gridpack::ComplexType* values) 
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
  printf("wsieg1 init: %f\t%f\t%f\t%f\t%f\t%f\n", x1LL, x2GovOut, x3Turb1, x4Turb2, x5Turb3, x6Turb4);
  values[0] = x1LL;
  values[1] = x2GovOut;
  values[2] = x3Turb1;
  values[3] = x4Turb2;
  values[4] = x5Turb3;
  values[5] = x6Turb4;
}

/**
 * Write output from governors to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Wsieg1Gov::serialWrite(char *string, const int bufsize,const char *signal)
{
}

double Wsieg1Gov::getAngle(void)
{
}

/**
 * Write out governor state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Wsieg1Gov::write(const char* signal, char* string)
{
}

/**
 *  Set the number of variables for this governor model
 *  @param [output] number of variables for this model
 */
bool Wsieg1Gov::vectorSize(int *nvar) const
{
  *nvar = 6;
  return true;
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto governors
 * @param values array containing governor state variables
*/
void Wsieg1Gov::setValues(gridpack::ComplexType *values)
{
  if(p_mode == XVECTOBUS) {
    x1LL = real(values[0]);
    x2GovOut = real(values[1]);
    x3Turb1 = real(values[2]);
    x4Turb2 = real(values[3]);
    x5Turb3 = real(values[4]);
    x6Turb4 = real(values[5]);
  } else if(p_mode == XDOTVECTOBUS) {
    dx1LL = real(values[0]);
    dx2GovOut = real(values[1]);
    dx3Turb1 = real(values[2]);
    dx4Turb2 = real(values[3]);
    dx5Turb3 = real(values[4]);
    dx6Turb4 = real(values[5]);
  }
}

/**
 * Return the values of the governor vector block
 * @param values: pointer to vector values
 * @return: false if governor does not contribute
 *        vector element
 */
bool Wsieg1Gov::vectorValues(gridpack::ComplexType *values)
{
  int x1_idx = 0;
  int x2_idx = 1;
  int x3_idx = 2;
  int x4_idx = 3;
  int x5_idx = 4;
  int x6_idx = 5;
  // On fault (p_mode == FAULT_EVAL flag), the governor variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    //values[delta_idx] = values[dw_idx] = 0.0;
    values[x1_idx] = 0.0;
    values[x2_idx] = 0.0;
    values[x3_idx] = 0.0;
    values[x4_idx] = 0.0;
    values[x5_idx] = 0.0;
    values[x6_idx] = 0.0;
  } else if(p_mode == RESIDUAL_EVAL) {
    // Governor equations
    //printf("...........%f\t%f\t%f\t%f\t%f\n", x1LL, x2GovOut, x3Turb1, x4Turb2, x5Turb3, x6Turb4);

    // State 1
    //printf("w = %f\n", w);
    double TempIn1 = K * w;//DBInt.Output(w);
    double TempOut;
    if (T1 > 4 * t_inc) {
        values[x1_idx] = (TempIn1 * ( 1 - T2 / T1) - x1LL) / T1 - dx1LL;
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
    values[x2_idx] = TempIn2 - dx2GovOut;
    //printf("TempIn1 = %f, TempOut = %f, w = %f, TempIn2 = %f\n", TempIn1, TempOut, w, TempIn2);
    // enforce non-windup limits
    if (dx2GovOut > 0 && x2GovOut >= Pmax) values[x2_idx] = - dx2GovOut;
    else if (dx2GovOut <0 && x2GovOut <= Pmin) values[x2_idx] = - dx2GovOut;
    // State 3
    double PGV = GainBlock.XtoY(GV);
    if (T4 < 4 * t_inc) {
        x3Turb1 = PGV;
        values[x3_idx] = - dx3Turb1;
    } else
        values[x3_idx] = (PGV - x3Turb1) / T4 - dx3Turb1;
    // State 4
    if (T5 < 4 * t_inc) {
        x4Turb2 = x3Turb1;
        values[x4_idx] = - dx4Turb2;
    } else
        values[x4_idx] = (x3Turb1 - x4Turb2) / T5 - dx4Turb2;
    // State 5
    if (T6 < 4 * t_inc) {
        x5Turb3 = x4Turb2;
        values[x5_idx] = - dx5Turb3;
    } else
        values[x5_idx] = (x4Turb2 - x5Turb3) / T6 - dx5Turb3;
    // State 6
    if (T7 < 4 * t_inc) {
        x6Turb4 = x5Turb3;
        values[x6_idx] = - dx6Turb4;
    } else
        values[x6_idx] = (x5Turb3 - x6Turb4) / T7 - dx6Turb4;

    ///printf("wsieg1 dx: %f\t%f\t%f\t%f\t%f\t%f\n", dx1LL, dx2GovOut, dx3Turb1, dx4Turb2, dx5Turb3, dx6Turb4);
    ///printf("wsieg1 x: %f\t%f\t%f\t%f\t%f\t%f\n", x1LL_1, x2GovOut_1, x3Turb1_1, x4Turb2_1, x5Turb3_1, x6Turb4_1);
  
    Pmech1 = x3Turb1 * K1 + x4Turb2 * K3 + x5Turb3 * K5 + x6Turb4 * K7;
    Pmech2 = x3Turb1 * K2 + x4Turb2 * K4 + x5Turb3 * K6 + x6Turb4 * K8;
    
    ///printf("wsieg1 Pmech1 = %f, Pmech2 = %f\n", Pmech1, Pmech2);      
  }
  
  return true;
}

/**
 * Return the governor current injection (in rectangular form) 
 * @param [output] IGD - real part of the governor current // SJin: match to Ir
 * @param [output] IGQ - imaginary part of the governor current // SJin: match to Ii 
*/
void Wsieg1Gov::getCurrent(double *IGD, double *IGQ)
{
}

/**
 * Return the matrix entries
 * @param [output] nval - number of values set
 * @param [output] row - row indics for matrix entries
 * @param [output] col - col indices for matrix entries
 * @param [output] values - matrix entries
 * return true when matrix entries set
 */
bool Wsieg1Gov::matrixDiagEntries(int *nval,int *row, int *col, gridpack::ComplexType *values)
{
  int idx = 0;
  if(p_mode == FAULT_EVAL) { // SJin: put values 1 along diagonals, 0 along off diagonals
  // On fault (p_mode == FAULT_EVAL flag), the governor variables are held constant. This is done by setting the diagonal matrix entries to 1.0 and all other entries to 0. The residual function values are already set to 0.0 in the vector values function. This results in the equation 1*dx = 0.0 such that dx = 0.0 and hence x does not get changed.
    row[idx] = 0; col[idx] = 0;
    values[idx] = 1.0;
    idx++;
    row[idx] = 1; col[idx] = 1;
    values[idx] = 1.0;
    idx++;
    row[idx] = 2; col[idx] = 2;
    values[idx] = 1.0;
    idx++;
    row[idx] = 3; col[idx] = 3;
    values[idx] = 1.0;
    idx++;
    row[idx] = 4; col[idx] = 4;
    values[idx] = 1.0;
    idx++;
    row[idx] = 5; col[idx] = 5;
    values[idx] = 1.0;
    idx++;
    row[idx] = 6; col[idx] = 6;
    values[idx] = 1.0;
    idx++;
    *nval = idx;
  } /*else if(p_mode == DIG_DV) { // SJin: Jacobian matrix block Jgy
    // These are the partial derivatives of the governor currents (see getCurrent function) w.r.t to the voltage variables VD and VQ

    *nval = idx;
  } else if(p_mode == DFG_DV) {  // SJin: Jacobian matrix block Jfyi
    // These are the partial derivatives of the governor equations w.r.t variables VD and VQ  

    *nval = idx;
  } else if(p_mode == DIG_DX) { // SJin: Jacobian matrix block Jgx
    // These are the partial derivatives of the governor currents (see getCurrent) w.r.t governor variables

    *nval = idx;
  } else { // SJin: Jacobin matrix block Jfxi
    // Partials of governor equations w.r.t governor variables
 
    *nval = idx;
  }*/
  return true;
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power 
 */
void Wsieg1Gov::setMechanicalPower(double pmech)
{
  Pmech1 = pmech;
  Pmech2 = pmech;
  printf("Pmech1 in WSIEG1 = %f\n", pmech);
}

/**
 * Set the rotor speed deviation parameter inside the governor
 * @param delta_o value of the rotor speed deviation 
 */
void Wsieg1Gov::setRotorSpeedDeviation(double delta_o)
{
  w = delta_o;
}

/** 
 * Get the value of the mechanical power parameter
 * @return value of the mechanical power 
 */
double Wsieg1Gov::getMechanicalPower()
{
  return Pmech1;
}

/** 
 * Get the value of the rotor speed deviation parameter
 * @return value of the rotor speed deviation 
 */
double Wsieg1Gov::getRotorSpeedDeviation()
{
  return w;
}

void Wsieg1Gov::setVcomp(double Vcomp)
{
}
   
/**
 * Set the value of the time step
 * @return value of the time step
 */
/*void Wsieg1Gov::setTimestep(double timestep)
{
}*/

/**
 * Set the value of the time increment 
 * @return value of the time increment
 */
/*void Wsieg1Gov::setTimeincrement(double timeincrement)
{
}*/

