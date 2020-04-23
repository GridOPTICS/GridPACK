/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wsieg1.cpp
 * @author Shuangshuang Jin
 * @author Shrirang Abhyankar
 * @Last modified:   01/02/20
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
  xLL = 0.0; 
  xGV = 0.0; 
  xT1 = 0.0; 
  xT2 = 0.0; 
  xT3 = 0.0; 
  xT4 = 0.0;
  dxLL = 0.0;
  dxGV = 0.0;
  dxT1 = 0.0;
  dxT2 = 0.0;
  dxT3 = 0.0;
  dxT4 = 0.0;
  xLLprev = 0.0;
  xGVprev = 0.0;
  xT1prev = 0.0;
  xT2prev = 0.0;
  xT3prev = 0.0;
  xT4prev = 0.0;
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

  nxgov = 6;
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
  data->getValue(BUS_NUMBER, &bid);
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

  //Set flags for differential or algebraic equations
  iseq_diff[0] = (T1 == 0 || T2 == 0)?0:1;
  iseq_diff[1] = 1;
  iseq_diff[2] = (T4 == 0)?0:1;
  iseq_diff[3] = (T5 == 0)?0:1;
  iseq_diff[4] = (T6 == 0)?0:1;
  iseq_diff[5] = (T7 == 0)?0:1;
}

/**
 * Initialize governor model before calculation
 * @param [output] values - array where initialized governor variables should be set
 */
void Wsieg1Gov::init(gridpack::ComplexType* values) 
{
  BaseGenModel* gen=getGenerator();
  double dw = gen->getRotorSpeedDeviation();
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
  xT4 = PGV;
  xT3 = PGV;
  xT2 = PGV;
  xT1 = PGV;
  xGV = PGV;

  Pref = PGV;

  // Note: (GV > Pmax) or (GV < Pmin) is an initial state violation
  if (Iblock == 1 && Pmin == 0) Pmin = PGV;
  if (Iblock == 2 && Pmax == 0) Pmax = PGV;
  if (Iblock == 3 && Pmin == 0) Pmin = PGV;
  if (Iblock == 3 && Pmax == 0) Pmax = PGV;
  if (iseq_diff[0]) xLL = K*dw * (1 - T2 / T1);
  else xLL = K*dw;

  values[0] = xLL;
  values[1] = xGV;
  values[2] = xT1;
  values[3] = xT2;
  values[4] = xT3;
  values[5] = xT4;
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
  return false;
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
  *nvar = nxgov;
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
    xLL = real(values[0]);
    xGV = real(values[1]);
    xT1 = real(values[2]);
    xT2 = real(values[3]);
    xT3 = real(values[4]);
    xT4 = real(values[5]);
  } else if(p_mode == XDOTVECTOBUS) {
    dxLL = real(values[0]);
    dxGV = real(values[1]);
    dxT1 = real(values[2]);
    dxT2 = real(values[3]);
    dxT3 = real(values[4]);
    dxT4 = real(values[5]);
  } else if (p_mode == XVECPRETOBUS) {
    xLLprev = real(values[0]);
    xGVprev = real(values[1]);
    xT1prev = real(values[2]);
    xT2prev = real(values[3]);
    xT3prev = real(values[4]);
    xT4prev = real(values[5]);
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
  double yLL;
  BaseGenModel* gen=getGenerator();
  double dw = gen->getRotorSpeedDeviation();

  // On fault (p_mode == FAULT_EVAL flag), the governor variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    if(iseq_diff[x1_idx]) values[x1_idx] = xLL - xLLprev; 
    else values[x1_idx] = -xLL + K*dw;
    yLL = xLL;

    values[x2_idx] = xGV - xGVprev;

    if(iseq_diff[x3_idx]) values[x3_idx] = xT1 - xT1prev;
    else values[x3_idx] = xGV - xT1;
      
    if(iseq_diff[x4_idx]) values[x4_idx] = xT2 - xT2prev;
    else values[x4_idx] = xT1 - xT2;

    if(iseq_diff[x5_idx]) values[x5_idx] = xT3 - xT3prev;
    else values[x5_idx] = xT2 - xT3;

    if(iseq_diff[x6_idx]) values[x6_idx] = xT4 - xT4prev;
    else values[x6_idx] = xT3 - xT4;

  } else if(p_mode == RESIDUAL_EVAL) {
    // State 1 xLL
    if(iseq_diff[x1_idx]) {
      values[x1_idx] = (-xLL + (1 - T2/T1)*K*dw)/T1 - dxLL;
      yLL = xLL + T2/T1*K*dw;
    } else {
      values[x1_idx] = -xLL + K*dw;
      yLL = xLL;
    }

    // State 2 xGV ignoring limits Uo, Uc
    //    values[x2_idx] = fmin(fmax(Uc,(Pref - xGV - yLL)/T3),Uo) -dxGV;
    values[x2_idx] = (Pref - xGV - yLL)/T3 - dxGV;
    
    // State 3 xT1
    if(iseq_diff[x3_idx]) values[x3_idx] = (xGV - xT1)/T4 - dxT1;
    else values[x3_idx] = xGV - xT1;

    // State 4 xT2
    if(iseq_diff[x4_idx]) values[x4_idx] = (xT1 - xT2)/T5 - dxT2;
    else values[x4_idx] = xT1 - xT2;

    // State 5 xT3
    if(iseq_diff[x5_idx]) values[x5_idx] = (xT2 - xT3)/T6 - dxT3;
    else values[x5_idx] = xT2 - xT3;

    // State 6 xT4
    if(iseq_diff[x6_idx]) values[x6_idx] = (xT3 - xT4)/T7 - dxT4;
    else values[x6_idx] = xT3 - xT4;
        
  }
  
  return true;
}

/**
 * Set Jacobian block
 * @param values a 2-d array of Jacobian block for the bus
 */
bool Wsieg1Gov::setJacobian(gridpack::ComplexType **values)
{
  int xLL_idx = offsetb;
  int xGV_idx = offsetb+1;
  int xT1_idx = offsetb+2;
  int xT2_idx = offsetb+3;
  int xT3_idx = offsetb+4;
  int xT4_idx = offsetb+5;
  double yLL;
  double dyLL_dxLL=0.0,dyLL_dxGV=0.0;
  double dyLL_dxT1=0.0,dyLL_dxT2=0.0;
  double dyLL_dxT3=0.0,dyLL_dxT4=0.0;
  double dyLL_ddw=0.0;
  int    dw_idx;

  dw_idx = getGenerator()->getRotorSpeedDeviationLocation();
  if(p_mode == FAULT_EVAL) {
    if(iseq_diff[0]) {
      values[xLL_idx][xLL_idx] = 1.0;
    } else {
      values[xLL_idx][xLL_idx] = -1.0;
      dw_idx = getGenerator()->getRotorSpeedDeviationLocation();
      values[dw_idx][xLL_idx] = K;
    }

    values[xGV_idx][xGV_idx] = 1.0;

    if(iseq_diff[2]) {
      values[xT1_idx][xT1_idx] = 1.0;
    } else {
      values[xGV_idx][xT1_idx] = 1.0;
      values[xT1_idx][xT1_idx] = -1.0;
    }

    if(iseq_diff[3]) {
      values[xT2_idx][xT2_idx] = 1.0;
    } else {
      values[xT1_idx][xT2_idx] = 1.0;
      values[xT2_idx][xT2_idx] = -1.0;
    }

    if(iseq_diff[4]) {
      values[xT3_idx][xT3_idx] = 1.0;
    } else {
      values[xT2_idx][xT3_idx] = 1.0;
      values[xT3_idx][xT3_idx] = -1.0;
    }

    if(iseq_diff[5]) {
      values[xT4_idx][xT4_idx] = 1.0;
    } else {
      values[xT3_idx][xT4_idx] = 1.0;
      values[xT4_idx][xT4_idx] = -1.0;
    }
  } else {
    if(iseq_diff[0]) {
      values[xLL_idx][xLL_idx] = -1.0/T1 - shift;

      values[dw_idx][xLL_idx] = (1 - T2/T1)*K/T1;
      dyLL_dxLL = 1.0;
      dyLL_ddw  = T2/T1*K;
    } else {
      values[xLL_idx][xLL_idx] = -1.0;

      values[dw_idx][xLL_idx] = K;
      dyLL_dxLL = 1.0;
    }

    values[xLL_idx][xGV_idx] = -dyLL_dxLL/T3;
    values[xGV_idx][xGV_idx] = -1.0/T3 - shift; 
    values[xT1_idx][xGV_idx] = -dyLL_dxT1/T3;
    values[xT2_idx][xGV_idx] = -dyLL_dxT2/T3;
    values[xT3_idx][xGV_idx] = -dyLL_dxT3/T3;
    values[xT4_idx][xGV_idx] = -dyLL_dxT4/T3;

    values[dw_idx][xGV_idx] = -dyLL_ddw/T3;
    
    if(iseq_diff[2]) {
      values[xGV_idx][xT1_idx] = 1/T4;
      values[xT1_idx][xT1_idx] = -1/T4 - shift;
    } else {
      values[xGV_idx][xT1_idx] = 1.0;
      values[xT1_idx][xT1_idx] = -1.0;
    }

    if(iseq_diff[3]) {
      values[xT1_idx][xT2_idx] = 1/T5;
      values[xT2_idx][xT2_idx] = -1/T5 - shift;
    } else {
      values[xT1_idx][xT2_idx] = 1.0;
      values[xT2_idx][xT2_idx] = -1.0;
    }

    if(iseq_diff[4]) {
      values[xT2_idx][xT3_idx] = 1/T6;
      values[xT3_idx][xT3_idx] = -1/T6 - shift;
    } else {
      values[xT2_idx][xT3_idx] = 1.0;
      values[xT3_idx][xT3_idx] = -1.0;
    }

    if(iseq_diff[5]) {
      values[xT3_idx][xT4_idx] = 1/T7;
      values[xT4_idx][xT4_idx] = -1/T7 - shift;
    } else {
      values[xT3_idx][xT4_idx] = 1.0;
      values[xT4_idx][xT4_idx] = -1.0;
    }
  }

  return true;
}


/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power 
 */
void Wsieg1Gov::setInitialMechanicalPower(double Pmech0)
{
  Pmech1 = Pmech0;
  Pmech2 = Pmech0;
}

/** 
 * Get the value of the mechanical power parameter
 * @return value of the mechanical power 
 */
double Wsieg1Gov::getMechanicalPower()
{
  Pmech1 = xT1 * K1 + xT2 * K3 + xT3 * K5 + xT4 * K7;
  Pmech2 = xT1 * K2 + xT2 * K4 + xT3 * K6 + xT4 * K8;

  return Pmech1;
}

/**
 * Partial derivatives of Mechanical Power Pmech w.r.t. governor variables
 * @param xgov_loc locations of governor variables
 * @param dPmech_dxgov partial derivatives of mechanical power Pmech w.r.t governor variables
*/
bool Wsieg1Gov::getMechanicalPowerPartialDerivatives(int *xgov_loc,double *dPmech_dxgov)
{
  int i;

  for(i=0; i < nxgov; i++) xgov_loc[i] = offsetb + i;

  dPmech_dxgov[0] = dPmech_dxgov[1] = 0.0;
  dPmech_dxgov[2] = K1; 
  dPmech_dxgov[3] = K3; 
  dPmech_dxgov[4] = K5; 
  dPmech_dxgov[5] = K7; 

  return true;
}

void Wsieg1Gov::setVcomp(double Vcomp)
{
}
