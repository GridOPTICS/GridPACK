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

#define TS_THRESHOLD 1

Esst1aExc::Esst1aExc(void)
{
  x1Va = 0.0; 
  x2Vcomp = 0.0; 
  x3LL1 = 0.0; 
  x4LL2 = 0.0; 
  x5Deriv = 0.0; 
  dx1Va = 0.0;
  dx2Vcomp = 0.0;
  dx3LL1 = 0.0;
  dx4LL2 = 0.0;
  dx5Deriv = 0.0;
  UEL = 0.0; 
  VOS = 0.0; 
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
  OptionToModifyLimitsForInitialStateLimitViolation = true;
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
  BaseExcModel::load(data,idx); // load parameters in base exciter model
  
  // load parameters for the model type
  if (!data->getValue(EXCITER_UEL, &UEL, idx)) UEL = 0.0; // TBD: UEL
  if (!data->getValue(EXCITER_VOS, &VOS, idx)) VOS = 0.0; // TBD: VOS
  if (!data->getValue(EXCITER_TR, &Tr, idx)) Tr = 0.0; // Tr
  if (!data->getValue(EXCITER_VIMAX, &Vimax, idx)) Vimax = 0.0; // TBD: Vimax
  if (!data->getValue(EXCITER_VIMIN, &Vimin, idx)) Vimin = 0.0; // TBD: Vimin
  if (!data->getValue(EXCITER_TC, &Tc, idx)) Tc = 0.0; // Tc
  if (!data->getValue(EXCITER_TB, &Tb, idx)) Tb = 0.0; // Tb
  if (!data->getValue(EXCITER_TC1, &Tc1, idx)) Tc1 = 0.0; // TBD: Tc1
  if (!data->getValue(EXCITER_TB1, &Tb1, idx)) Tb1 = 0.0; // TBD: Tb1
  if (!data->getValue(EXCITER_KA, &Ka, idx)) Ka = 0.0; // Ka
  if (!data->getValue(EXCITER_TA, &Ta, idx)) Ta = 0.0; // Ta
  if (!data->getValue(EXCITER_VAMAX, &Vamax, idx)) Vamax = 0.0; // TBD: Vamax
  if (!data->getValue(EXCITER_VAMIN, &Vamin, idx)) Vamin = 0.0; // TBD: Vamin
  if (!data->getValue(EXCITER_VRMAX, &Vrmax, idx)) Vrmax = 0.0; // Vrmax
  if (!data->getValue(EXCITER_VRMIN, &Vrmin, idx)) Vrmin = 0.0; // Vrmin
  if (!data->getValue(EXCITER_KC, &Kc, idx)) Kc = 0.0; // TBD: Kc
  if (!data->getValue(EXCITER_KF, &Kf, idx)) Kf = 0.0; // Kf
  if (!data->getValue(EXCITER_TF, &Tf, idx)) Tf = 0.0; // Tf
  if (!data->getValue(EXCITER_KLR, &Klr, idx)) Klr = 0.0; // TBD: Klr
  if (!data->getValue(EXCITER_ILR, &Ilr, idx)) Ilr = 0.0; // TBD: Ilr
  //if (!data->getValue(EXCITER_TA1, &Ta1, idx)) Ta1 = 0.0; // Ta1

  //printf("UEL=%f,VOS=%f,Tr=%f,Vimax=%f,Vimin=%f,Tc=%f,Tb=%f,Tc1=%f,Tb1=%f,Ka=%f,Ta=%f,Vamax=%f,Vamin=%f,Vrmax=%f,Vrmin=%f,Kc=%f,Kf=%f,Tf=%f,Klr=%f,Ilr=%f\n",UEL,VOS,Tr,Vimax,Vimin,Tc,Tb,Tc1,Tb1,Ka,Ta,Vamax,Vamin,Vrmax,Vrmin,Kc,Kf,Tf,Klr,Ilr);

  // Convert exciter parameters from machine base to MVA base
  /*p_H *= mbase/sbase;
  p_D *= mbase/sbase;
  p_Xdp /= mbase/sbase;*/

}

/**
 * Saturation function
 * @ param x
 */
double Esst1aExc::Sat(double x)
{

}

/**
 * Initialize exciter model before calculation
 * @param [output] values - array where initialized exciter variables should be set
 */
void Esst1aExc::init(gridpack::ComplexType* values) 
{
  //double Vm = sqrt(VD*VD + VQ*VQ); // SJin: voltage VD and VQ come from base_exc_model.hpp
  //double mag = Vm;
   //printf("esst1a: Efd = %f, Vm=%f\n", Efd,Vm);
  double mag = Vterm;
  //printf("esst1a: Efd = %f, mag=%f\n", Efd,mag);

  // Parameter cleanup
  // Just to make the code simpler below we will do the following cleanup 
  // on the input parameters and also store some variables.
  if (Tb == 0) Tc = 0;
  if (Tb1 == 0) Tc1 = 0;
  if (Tf == 0) Kf = 0;
  // Following is needed to avoid an algebraic loop in the model
  if (Kf != 0 && Ta < TS_THRESHOLD * ts) Ta = TS_THRESHOLD * ts;

  //Vterm = mag;
  //presentMag = mag;
  //Theta = ang;
  //presentAng = ang;
  //printf("esst1a init: Efd = %f\n", Efd);
  //State 1
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (Efd > (Vterm * Vrmax - Kc * LadIfd)) Vrmax = (Efd + Kc * LadIfd) / Vterm;
    if (Efd < (Vterm * Vrmin)) Vrmax = Efd / Vterm;
  }
  //printf("/////Efd=%f,Kc=%f,LadIfd=%f,Vterm=%f, Vrmax=%f\n",Efd,Kc,LadIfd,Vterm,Vrmax);
  double Va;
  if (LadIfd > Ilr) Va = Efd + Klr * (LadIfd - Ilr);
  else Va = Efd;
  // Check for initial limit violations - should not happen!
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (Va > Vamax) Vamax = Va;
    if (Va < Vamin) Vamin = Va;
  }
  x1Va = Va;
  //State 4
  double TempLL = Va / Ka;
  if (Tb1 == 0) x4LL2 = TempLL;
  else if (Tb1 < TS_THRESHOLD * ts) x4LL2 = 0;
  else x4LL2 = TempLL * (1 - Tc1 / Tb1);
  //State 3
  if (Tb == 0) x3LL1 = TempLL;
  else if (Tb < TS_THRESHOLD * ts) x3LL1 = 0;
  else x3LL1 = TempLL * (1 - Tc / Tb);
  //State 2
  if (OptionToModifyLimitsForInitialStateLimitViolation) {
    if (TempLL > Vimax) Vimax = TempLL;
    if (TempLL < Vimin) Vimin = TempLL;
  }
  //printf ("Esst1a ini() Vcomp = %f \n", Vcomp);
  x2Vcomp = Vcomp;
  //State 5
  x5Deriv = 0;
  Vref = Vcomp + TempLL;

  //printf("esst1a init:  %f\t%f\t%f\t%f\t%f\t%f\n", x1Va, x2Vcomp, x3LL1, x4LL2, x5Deriv, Vref);
  values[0] = x1Va;
  values[1] = x2Vcomp;
  values[2] = x3LL1;
  values[3] = x4LL2;
  values[4] = x5Deriv;
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
    x1Va = real(values[0]);
    x2Vcomp = real(values[1]);
    x3LL1 = real(values[2]);
    x4LL2 = real(values[3]);
    x5Deriv = real(values[4]);
  } else if(p_mode == XDOTVECTOBUS) {
    dx1Va = real(values[0]);
    dx2Vcomp = real(values[1]);
    dx3LL1 = real(values[2]);
    dx4LL2 = real(values[3]);
    dx5Deriv = real(values[4]);
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
  // On fault (p_mode == FAULT_EVAL flag), the exciter variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    //values[delta_idx] = values[dw_idx] = 0.0;
    values[x1_idx] = 0.0;
    values[x2_idx] = 0.0;
    values[x3_idx] = 0.0;
    values[x4_idx] = 0.0;
    values[x5_idx] = 0.0;
  } else if(p_mode == RESIDUAL_EVAL) {
    //printf("\n esst1a: what's the initial values for the first iteration?\n");
    //printf("essti1 in values: %f\t%f\t%f\t%f\t%f\t%f\n", x1Va, x2Vcomp, x3LL1, x4LL2, x5Deriv, Vref);
    // Exciter equations
    //printf("...........%f\t%f\t%f\t%f\t%f\n", x1Va, x2Vcomp, x3LL1, x4LL2, x5Deriv);
    //printf(".........Vterminal = %f\n", Vterminal);
    /*bool flag2 = false; // set flags for residual function condisition to be used for Jacobian matrix calculation 
    bool flag3 = false;
    bool flag4 = false;
    bool flag5 = false;*/
    //printf("HIHIHI from esst1a\n"); 
    // State 2
    if (Tr < TS_THRESHOLD * t_inc) {
      x2Vcomp = Vcomp; // Must propogate input value instantaneously
      values[x2_idx] = - dx2Vcomp;
    } else {
      values[x2_idx] = 1 / Tr * (Vcomp - x2Vcomp) - dx2Vcomp;
    }
    // State 5
    double TempIn;
    if (Kf > 0) {
      TempIn = Klr * (LadIfd - Ilr);
      if (TempIn < 0) TempIn = 0;
      // Note: if Ta = 0, then we would have an algebriac loop here which
      // could cause numerical troubles. Thus we enforce Ta > 0 if (Kf > 0 and Tf > 0)
      TempIn = x1Va - TempIn;
      //if ("VOS at Output") TempIn = TempIn + Vstab; // TBD: "VOS at Output"???
      // Ignore Over and Under Excitation Limit for now
      double UseTf;
      if (Tf > TS_THRESHOLD * t_inc) UseTf = Tf;
      else UseTf = TS_THRESHOLD * t_inc;
      values[x5_idx] = (TempIn * (-Kf / UseTf) - x5Deriv) / UseTf - dx5Deriv;
      TempIn = TempIn * Kf / UseTf + x5Deriv; // Input into summation block
    } else {
      values[x5_idx] = - dx5Deriv;
      TempIn = 0; // Input into summation block;
    }
    // State 3
    //printf("x2Vcomp=%f,TempIn=%f,Vref=%f\n",x2Vcomp,TempIn,Vref);
    TempIn = - x2Vcomp - TempIn + Vref;
    //if ("VOS at Input") TempIn = TempIn + Vstab; // TBD: "VOS at Input"???
    if (TempIn > Vimax) TempIn = Vimax;
    if (TempIn < Vimin) TempIn = Vimin;
    // Ignore Under Excitation Limit
    if (Tb < TS_THRESHOLD * t_inc) {
      values[x3_idx] = - dx3LL1;
      TempIn = TempIn; // output of first lead lag - just pass input
      //printf("TempIn 3333.1\n");
    } else {
      values[x3_idx] = (TempIn * (1 - Tc / Tb) - x3LL1)/Tb - dx3LL1;
      TempIn = TempIn * Tc / Tb + x3LL1; // output of first lead lag
      //printf("TempIn 3333.2\n");
    }
    // State 4
    if (Tb1 < TS_THRESHOLD * t_inc) {
    //    printf ("entering Tb < 4h\n");
      values[x4_idx] = - dx4LL2;
      TempIn = TempIn; // output of second lead lag - just pass input
    } else {
      values[x4_idx] = (TempIn * (1 - Tc1 / Tb1) - x4LL2)/Tb1 - dx4LL2;
      TempIn = TempIn * Tc1 / Tb1 + x4LL2; // output of second lead lag
    }
    // State 1
    if (Ta < TS_THRESHOLD * t_inc) x1Va = Ka * TempIn;
    if (x1Va > Vamax) x1Va = Vamax; 
    if (x1Va < Vamin) x1Va = Vamin; 
    if (Ta < TS_THRESHOLD * t_inc) values[x1_idx] = - dx1Va;
    else {
       values[x1_idx] = (Ka * TempIn - x1Va) / Ta - dx1Va;
    }
    if (dx1Va > 0 && x1Va >= Vamax) values[x1_idx] = - dx1Va;
    if (dx1Va < 0 && x1Va <= Vamin) values[x1_idx] = - dx1Va; 
     
    // Update Efd
    double Temp = Klr * (LadIfd - Ilr);
    if (Temp < 0) Temp = 0;
    Temp = x1Va - Temp;
    // Ignore Over and Under Excitation Limit for now
    if (Temp > (Vterm * Vrmax - Kc * LadIfd)) {
       Temp = Vterm * Vrmax - Kc * LadIfd;
    }
    if (Temp < (Vterm * Vrmin)) {
       Temp = Vterm * Vrmin;
    }
    Efd = Temp;
    //printf("x1Va = %f, Temp=%f\n", x1Va, Temp);

    /*values[x1_idx] = 21;
    values[x2_idx] = 22;
    values[x3_idx] = 23;
    values[x4_idx] = 24;
    values[x5_idx] = 25;*/

    //printf("esst1a: %f\t%f\t%f\t%f\t%f\n", real(values[x1_idx]),real(values[x2_idx]),real(values[x3_idx]),real(values[x4_idx]),real(values[x5_idx]));
    //printf("esst1aidx: %d %d %d %d %d\n",x1_idx,x2_idx,x3_idx,x4_idx,x5_idx);
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
 * Set the field voltage parameter inside the exciter
 * @param fldv value of the field voltage
 */
void Esst1aExc::setFieldVoltage(double fldv)
{
  Efd = fldv;
  //printf("Efd in Esst1a = %f\n", Efd);
}

/**
 * Set the field current parameter inside the exciter
 * @param fldc value of the field current
 */
void Esst1aExc::setFieldCurrent(double fldc)
{
  LadIfd = fldc;
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double Esst1aExc::getFieldVoltage()
{
  return Efd;
}

/** 
 * Get the value of the field current parameter
 * @return value of field current
 */
double Esst1aExc::getFieldCurrent()
{
  return 0.0;
}

/** 
 * Set the value of the Vterminal
 * @return value of field current
 */
void Esst1aExc::setVterminal(double mag)
{
  Vterm = mag;
}

void Esst1aExc::setVcomp(double vtmp)
{
  Vcomp = vtmp;
  //printf("Vcomp has been set to %f\n", Vcomp);
}
   
/**
 * Set the value of the time step
 * @return value of the time step
 */
/*void Esst1aExc::setTimestep(double timestep)
{
  ts = timestep;
  t_inc = timestep;
}*/

/**
 * Set the value of the time increment 
 * @return value of the time increment
 */
//void Esst1aExc::setTimeincrement(double timeincrement)
//{
//}

