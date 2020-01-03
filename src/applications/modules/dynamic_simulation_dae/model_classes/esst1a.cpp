/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   esst1a.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   10/18/19
 *  
 * @brief  
 *
 *
 */

#include <esst1a.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Esst1aExc::Esst1aExc(void)
{
  Vmeas = 0.0; 
  dVmeas = 0.0;
  xLL1 = 0.0; 
  dxLL1 = 0.0;
  xLL2  = 0.0;
  dxLL2 = 0.0;
  Va    = 0.0;
  dVa   = 0.0;
  xf    = 0.0;
  dxf   = 0.0;
  UEL = 0; 
  VOS = 0; 
  Tr = 0.0; 
  Vimax = 0.0; 
  Vimin = 0.0; 
  Tc = 0.0; 
  Tb = 0.0;
  Tc1 = 0.0;
  Tb1 = 0.0; 
  Ka = 0.0;
  Ta = 0.0;
  Vamax = 0.0;
  Vamin = 0.0;
  Vrmax = 0.0;
  Vrmin = 0.0;
  Kc = 0.0;
  Kf = 0.0;
  Tf = 0.0;
  Klr = 0.0;
  Ilr = 0.0;    
}

Esst1aExc::~Esst1aExc(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to BaseExciterModel
 */
void Esst1aExc::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  data->getValue(BUS_NUMBER, &bid);
  BaseExcModel::load(data,idx); // load parameters in base exciter model
  
  // load parameters for the model type
  data->getValue(EXCITER_UEL, &UEL, idx);
  data->getValue(EXCITER_VOS, &VOS, idx);
  data->getValue(EXCITER_TR, &Tr, idx);
  data->getValue(EXCITER_VIMAX, &Vimax, idx);
  data->getValue(EXCITER_VIMIN, &Vimin, idx);
  data->getValue(EXCITER_TC, &Tc, idx);
  data->getValue(EXCITER_TB, &Tb, idx);
  data->getValue(EXCITER_TC1, &Tc1, idx);
  data->getValue(EXCITER_TB1, &Tb1, idx);
  data->getValue(EXCITER_KA, &Ka, idx);
  data->getValue(EXCITER_TA, &Ta, idx);
  data->getValue(EXCITER_VAMAX, &Vamax, idx);
  data->getValue(EXCITER_VAMIN, &Vamin, idx);
  data->getValue(EXCITER_VRMAX, &Vrmax, idx);
  data->getValue(EXCITER_VRMIN, &Vrmin, idx);
  data->getValue(EXCITER_KC, &Kc, idx);
  data->getValue(EXCITER_KF, &Kf, idx);
  data->getValue(EXCITER_TF, &Tf, idx);
  data->getValue(EXCITER_KLR, &Klr, idx);
  data->getValue(EXCITER_ILR, &Ilr, idx);

  // Set flags for differential or algebraic equations
  iseq_diff[0] = (Tr == 0)?0:1;
  iseq_diff[1] = (Tb == 0 || Tc == 0)?0:1;
  iseq_diff[2] = (Tb1 == 0 || Tc1 == 0)?0:1;
  iseq_diff[3] = (Ta == 0)?0:1;
  iseq_diff[4] = 1; // Tf is always > 0

  //printf("UEL=%f,VOS=%f,Tr=%f,Vimax=%f,Vimin=%f,Tc=%f,Tb=%f,Tc1=%f,Tb1=%f,Ka=%f,Ta=%f,Vamax=%f,Vamin=%f,Vrmax=%f,Vrmin=%f,Kc=%f,Kf=%f,Tf=%f,Klr=%f,Ilr=%f\n",UEL,VOS,Tr,Vimax,Vimin,Tc,Tb,Tc1,Tb1,Ka,Ta,Vamax,Vamin,Vrmax,Vrmin,Kc,Kf,Tf,Klr,Ilr);

}

/**
 * Initialize exciter model before calculation
 * @param [output] values - array where initialized exciter variables should be set
 */
void Esst1aExc::init(gridpack::ComplexType* values) 
{
  double Ec = sqrt(VD*VD + VQ*VQ);
  double yLL2,yLL1;
  double Vf=0.0,Vfd;
  BaseGenModel *gen=getGenerator();
  double LadIfd = gen->getFieldCurrent();

  // Field voltage (Efd) and bus voltage (VD,VQ) are already set 
  // Need to set the initial values for all the state variables

  Ec = sqrt(VD*VD + VQ*VQ);
  Vfd = Klr*(LadIfd - Ilr); 
  Vmeas    = Ec;
  xf       = -Kf/Tf*Efd;
  Va       = Efd + Vfd;
  yLL2     = Va/Ka;
  yLL1     = yLL2;
  if(iseq_diff[2]) xLL2    = (1 - Tc1/Tb1)*yLL2;
  else xLL2 = yLL2;
  Vref     = yLL1 + Vmeas + Vf;
  if(iseq_diff[1]) xLL1    = (1 - Tc/Tb)*(Vref - Vmeas - Vf);
  else xLL1 = Vref - Vmeas - Vf;

  values[0] = Vmeas;
  values[1] = xLL1;
  values[2] = xLL2;
  values[3] = Va;
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
bool Esst1aExc::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out exciter state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Esst1aExc::write(const char* signal, char* string)
{
}

/**
 *  Set the number of variables for this exciter model
 *  @param [output] number of variables for this model
 */
bool Esst1aExc::vectorSize(int *nvar) const
{
  *nvar = 5;
  return true;
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto exciters
 * @param values array containing exciter state variables
*/
void Esst1aExc::setValues(gridpack::ComplexType *values)
{  
  if(p_mode == XVECTOBUS) {
    Vmeas = real(values[0]);
    xLL1 = real(values[1]);
    xLL2 = real(values[2]);
    Va   = real(values[3]);
    xf   = real(values[4]);
  } else if(p_mode == XDOTVECTOBUS) {
    dVmeas = real(values[0]);
    dxLL1 = real(values[1]);
    dxLL2 = real(values[2]);
    dVa   = real(values[3]);
    dxf   = real(values[4]);
  } else if(p_mode == XVECPRETOBUS) {
    Vmeasprev = real(values[0]);
    xLL1prev = real(values[1]);
    xLL2prev = real(values[2]);
    Vaprev   = real(values[3]);
    xfprev   = real(values[4]);
  }
}

/**
 * Return the values of the exciter vector block
 * @param values: pointer to vector values
 * @return: false if exciter does not contribute
 *        vector element
 */
bool Esst1aExc::vectorValues(gridpack::ComplexType *values)
{
  int x1_idx = 0;
  int x2_idx = 1;
  int x3_idx = 2;
  int x4_idx = 3;
  int x5_idx = 4;
  double Ec,yLL1,yLL2,Vf;
  Ec = sqrt(VD*VD + VQ*VQ);
  
  Efd = Esst1aExc::getFieldVoltage();
  // On fault (p_mode == FAULT_EVAL flag), the exciter variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    // State 1 Vmeas
    if(iseq_diff[0]) values[0] = Vmeas - Vmeasprev;
    else values[0] = -Vmeas + Ec;

    // State 2 xLL1
    Vf = xf + Kf/Tf*Efd;
    if(iseq_diff[1]) {
      values[1] = xLL1 - xLL1prev;
      yLL1 = xLL1 + Tc/Tb*(Vref - Vmeas - Vf);
    } else {
      values[1] = -xLL1 + Vref - Vmeas - Vf;
      yLL1 = xLL1;
    }

    // State 3 xLL2
    if(iseq_diff[2]) {
      values[2] = xLL2 - xLL2prev;
      yLL2 = xLL2 + Tc1/Tb1*yLL1;
    } else {
      values[2] = -xLL2 + yLL1;
      yLL2 = xLL2;
    }

    // State 4 Va
    if(iseq_diff[3]) values[3] = Va - Vaprev;
    else values[3] = -Va + Ka*yLL2;

    // State 5 xf
    values[4] = xf - xfprev;

  } else if(p_mode == RESIDUAL_EVAL) {
    //    if (bid == 1) printf("\t\t%f\t%f\t%f\t%f\t%f\n", Vmeas, xLL1, xLL2, Va, xf);    

    // State 1 Vmeas
    if(iseq_diff[0]) values[0] = (-Vmeas + Ec)/Tr - dVmeas;
    else values[0] = -Vmeas + Ec;

    // State 2 xLL1
    Vf = xf + Kf/Tf*Efd;
    if(iseq_diff[1]) {
      values[1] = (-xLL1 + (1 - Tc/Tb)*(Vref - Vmeas - Vf))/Tb - dxLL1;
      yLL1 = xLL1 + Tc/Tb*(Vref - Vmeas - Vf);
    } else {
      values[1] = -xLL1 + Vref - Vmeas - Vf;
      yLL1 = xLL1;
    }

    // State 3 xLL2
    if(iseq_diff[2]) {
      values[2] = (-xLL2 + (1 - Tc1/Tb1)*yLL1)/Tb1 - dxLL2;
      yLL2 = xLL2 + Tc1/Tb1*yLL1;
    } else {
      values[2] = -xLL2 + yLL1;
      yLL2 = xLL2;
    }

    // State 4 Va
    if(iseq_diff[3]) values[3] = (-Va + Ka*yLL2)/Ta - dVa;
    else values[3] = -Va + Ka*yLL2;

    // State 5 xf
    values[4] = (-xf - Kf/Tf*Efd)/Tf - dxf;
    
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
bool Esst1aExc::matrixDiagEntries(int *nval,int *row, int *col, gridpack::ComplexType *values)
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
void Esst1aExc::setInitialFieldVoltage(double fldv)
{
  Efd = fldv;
  //printf("Efd in Esst1a = %f\n", Efd);
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double Esst1aExc::getFieldVoltage()
{
  BaseGenModel* gen = getGenerator();
  LadIfd = gen->getFieldCurrent();
  return Va - Klr*(LadIfd - Ilr); // should be actually max(0,Klr*(LadIfd - Ilr));
}
