/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   exdc1.cpp
 * @author Shrirang Abhyankar 
 * @Last modified:   01/02/20
 *  
 * @brief  
 *
 *
 */

#include <exdc1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Exdc1Exc::Exdc1Exc(void)
{
  Vmeas = 0.0; 
  dVmeas = 0.0;
  xLL = 0.0; 
  dxLL = 0.0;
  VR  = 0.0;
  dVR = 0.0;
  Efd    = 0.0;
  dEfd   = 0.0;
  xf    = 0.0;
  dxf   = 0.0;
  TR = 0.0; 
  Vrmax = 0.0; 
  Vrmin = 0.0; 
  TC = 0.0; 
  TB = 0.0;
  KA = 0.0;
  TA = 0.0;
  KF = 0.0;
  TF = 0.0;
  SWITCH = 0;
}

Exdc1Exc::~Exdc1Exc(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to BaseExciterModel
 */
void Exdc1Exc::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  data->getValue(BUS_NUMBER, &bid);
  BaseExcModel::load(data,idx); // load parameters in base exciter model
  if (!data->getValue(EXCITER_TR, &TR, idx)) TR = 0.0; // TR
  if (!data->getValue(EXCITER_KA, &KA, idx)) KA = 0.0; // KA
  if (!data->getValue(EXCITER_TA, &TA, idx)) TA = 0.0; // TA
  if (!data->getValue(EXCITER_TB, &TB, idx)) TB = 0.0; // TB
  if (!data->getValue(EXCITER_TC, &TC, idx)) TC = 0.0; // TC
  if (!data->getValue(EXCITER_VRMAX, &Vrmax, idx)) Vrmax = 0.0; // Vrmax
  if (!data->getValue(EXCITER_VRMIN, &Vrmin, idx)) Vrmin = 0.0; // Vrmin
  if (!data->getValue(EXCITER_KE, &KE, idx)) KE = 0.0; // KE
  if (!data->getValue(EXCITER_TE, &TE, idx)) TE = 0.0; // TE
  if (!data->getValue(EXCITER_KF, &KF, idx)) KF = 0.0; // KF
  if (!data->getValue(EXCITER_TF1, &TF, idx)) TF = 0.0; // TF
  if (!data->getValue(EXCITER_SWITCH, &SWITCH, idx)) SWITCH = 0.0; // SWITCH
  if (!data->getValue(EXCITER_E1, &E1, idx)) E1 = 0.0; // E1
  if (!data->getValue(EXCITER_SE1, &SE1, idx)) SE1 = 0.0; // SE1
  if (!data->getValue(EXCITER_E2, &E2, idx)) E2 = 0.0; // E2
  if (!data->getValue(EXCITER_SE2, &SE2, idx)) SE2 = 0.0; // SE2

  // Set flags for differential or algebraic equations
  iseq_diff[0] = (TR == 0)?0:1;
  iseq_diff[1] = (TB == 0 || TC == 0)?0:1;
  iseq_diff[2] = (TA == 0)?0:1;
  iseq_diff[3] = 1; // TE is always > 0
  iseq_diff[4] = 1; //TF always > 0
}

/**
 * Initialize exciter model before calculation
 * @param [output] values - array where initialized exciter variables should be set
 */
void Exdc1Exc::init(gridpack::ComplexType* values) 
{
  double Ec = sqrt(VD*VD + VQ*VQ);
  double yLL;
  double Vf=0.0;

  // Efd is already set by the generator model 
  Vmeas    = Ec;
  xf       = -KF/TF*Efd;

  // Assuming saturation SE(Efd) = 0.0 during initialization
  if(KE != 0.0) VR  = Efd/KE;
  else VR = Efd;

  yLL     = VR/KA;

  Vf   = xf + KF/TF*Efd;
  Vref = yLL + Vmeas + Vf;

  if(iseq_diff[1]) xLL    = (1 - TC/TB)*(Vref - Vmeas - Vf);
  else xLL = Vref - Vmeas - Vf;

  values[0] = Vmeas;
  values[1] = xLL;
  values[2] = VR;
  values[3] = Efd;
  values[4] = xf;
}

/**
 * Write output from exciters to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Exdc1Exc::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out exciter state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Exdc1Exc::write(const char* signal, char* string)
{
}

/**
 *  Set the number of variables for this exciter model
 *  @param [output] number of variables for this model
 */
bool Exdc1Exc::vectorSize(int *nvar) const
{
  *nvar = 5;
  return true;
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto exciters
 * @param values array containing exciter state variables
*/
void Exdc1Exc::setValues(gridpack::ComplexType *values)
{  
  if(p_mode == XVECTOBUS) {
    Vmeas = real(values[0]);
    xLL = real(values[1]);
    VR = real(values[2]);
    Efd   = real(values[3]);
    xf   = real(values[4]);
  } else if(p_mode == XDOTVECTOBUS) {
    dVmeas = real(values[0]);
    dxLL = real(values[1]);
    dVR = real(values[2]);
    dEfd   = real(values[3]);
    dxf   = real(values[4]);
  } else if(p_mode == XVECPRETOBUS) {
    Vmeasprev = real(values[0]);
    xLLprev = real(values[1]);
    VRprev = real(values[2]);
    Efdprev   = real(values[3]);
    xfprev   = real(values[4]);
  }
}

/**
 * Return the values of the exciter vector block
 * @param values: pointer to vector values
 * @return: false if exciter does not contribute
 *        vector element
 */
bool Exdc1Exc::vectorValues(gridpack::ComplexType *values)
{
  int x1_idx = 0;
  int x2_idx = 1;
  int x3_idx = 2;
  int x4_idx = 3;
  int x5_idx = 4;
  double Ec,yLL,Vf;
  Ec = sqrt(VD*VD + VQ*VQ);
  
  Efd = Exdc1Exc::getFieldVoltage();
  // On fault (p_mode == FAULT_EVAL flag), the exciter variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    // State 1 Vmeas
    if(iseq_diff[0]) values[0] = Vmeas - Vmeasprev;
    else values[0] = -Vmeas + Ec;


    Vf = xf + KF/TF*Efd;

    // State 2 xLL
    if(iseq_diff[1]) {
      values[1] = xLL - xLLprev;
      yLL = xLL + TC/TB*(Vref - Vmeas - Vf);
    } else {
      values[1] = -xLL + Vref - Vmeas - Vf;
      yLL = xLL;
    }

    // State 3 VR
    if(iseq_diff[2]) values[2] = VR - VRprev;
    else values[2] = -VR + KA*yLL;

    // State 4 Efd
    values[3] = Efd - Efdprev;

    // State 5 xf
    values[4] = xf - xfprev;

  } else if(p_mode == RESIDUAL_EVAL) {
    // State 1 Vmeas
    if(iseq_diff[0]) values[0] = (-Vmeas + Ec)/TR - dVmeas;
    else values[0] = -Vmeas + Ec;

    // State 2 xLL
    Vf = xf + KF/TF*Efd;
    if(iseq_diff[1]) {
      values[1] = (-xLL + (1 - TC/TB)*(Vref - Vmeas - Vf))/TB - dxLL;
      yLL = xLL + TC/TB*(Vref - Vmeas - Vf);
    } else {
      values[1] = -xLL + Vref - Vmeas - Vf;
      yLL = xLL;
    }

    // State 3 VR
    if(iseq_diff[2]) values[2] = (-VR + KA*yLL)/TA - dVR;
    else values[2] = -VR + KA*yLL;

    // State 4 Efd // Need to add saturation here
    values[3] = (VR - KE*Efd)/TE - dEfd;

    // State 5 xf
    values[4] = (-xf - KF/TF*Efd)/TF - dxf;
    
  }
  
  return true;
}

/**
 * Return the matrix entries
 * @param [output] nval - number of values set
 * @param [output] row - row indics for matrix entries
 * @param [output] col - col indices for matrix entries
 * @param [output] values - matrix entries
 * return true when matrix entries set
 */
bool Exdc1Exc::matrixDiagEntries(int *nval,int *row, int *col, gridpack::ComplexType *values)
{
  int idx = 0;
  if(p_mode == FAULT_EVAL) { // SJin: put values 1 along diagonals, 0 along off diagonals
  // On fault (p_mode == FAULT_EVAL flag), the exciter variables are held constant. This is done by setting the diagonal matrix entries to 1.0 and all other entries to 0. The residual function values are already set to 0.0 in the vector values function. This results in the equation 1*dx = 0.0 such that dx = 0.0 and hence x does not get changed.
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
    *nval = idx;
  } /*else if(p_mode == DIG_DV) { // SJin: Jacobian matrix block Jgy
    // These are the partial derivatives of the exciter currents (see getCurrent function) w.r.t to the voltage variables VD and VQ
    
    *nval = idx;
  } else if(p_mode == DFG_DV) {  // SJin: Jacobian matrix block Jfyi
    // These are the partial derivatives of the exciter equations w.r.t variables VD and VQ  

    *nval = idx;
  } else if(p_mode == DIG_DX) { // SJin: Jacobian matrix block Jgx
    // These are the partial derivatives of the exciter currents (see getCurrent) w.r.t exciter variables

    *nval = idx;
  } else { // SJin: Jacobin matrix block Jfxi
    // Partials of exciter equations w.r.t exciter variables
 
    *nval = idx;
  }*/
  return true;
}

/**
 * Set the initial field voltage (at t = tstart) for the exciter
 * @param fldv value of the field voltage
 */
void Exdc1Exc::setInitialFieldVoltage(double fldv)
{
  Efd = fldv;
  //printf("Efd in Exdc1 = %f\n", Efd);
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double Exdc1Exc::getFieldVoltage()
{
  return Efd;
}

