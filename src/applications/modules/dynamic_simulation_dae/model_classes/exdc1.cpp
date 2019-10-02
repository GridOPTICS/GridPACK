/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   exdc1.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   08/23/19
 *  
 * @brief  
 *
 *
 */

#include <exdc1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

#define TS_THRESHOLD 1

Exdc1Exc::Exdc1Exc(void)
{
  x1 = 0.0; 
  x2 = 0.0; 
  x3 = 0.0; 
  x4 = 0.0; 
  x5 = 0.0; 
  dx1 = 0.0;
  dx2 = 0.0;
  dx3 = 0.0;
  dx4 = 0.0;
  dx5 = 0.0;
  TR = 0.0; 
  KA = 0.0; 
  TA = 0.0; 
  TB = 0.0; 
  TC = 0.0; 
  Vrmax = 0.0; 
  Vrmin = 0.0;
  KE = 0.0;
  TE = 0.0; 
  KF = 0.0;
  TF = 0.0;
  SWITCH = 0.0;
  E1 = 0.0;
  SE1 = 0.0;
  E2 = 0.0;
  SE2 = 0.0;
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
  BaseExcModel::load(data,idx); // load parameters in base exciter model
  
  // load parameters for the model type
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
  //if (!data->getValue(EXCITER_TF1, &TF1, idx)) TF1 = 0.0; // TF1
  if (!data->getValue(EXCITER_TF1, &TF, idx)) TF = 0.0; // TF
  //printf("load TF = %f\n", TF);
  if (!data->getValue(EXCITER_SWITCH, &SWITCH, idx)) SWITCH = 0.0; // SWITCH
  if (!data->getValue(EXCITER_E1, &E1, idx)) E1 = 0.0; // E1
  if (!data->getValue(EXCITER_SE1, &SE1, idx)) SE1 = 0.0; // SE1
  if (!data->getValue(EXCITER_E2, &E2, idx)) E2 = 0.0; // E2
  if (!data->getValue(EXCITER_SE2, &SE2, idx)) SE2 = 0.0; // SE2

  //printf("TR=%f,KA=%f,TA=%f,TB=%f,TC=%f,Vrmax=%f,Vrmin=%f,KE=%f,TE=%f,KF=%f,TF=%f,SWITCH=%f,E1=%f,SE1=%f,E2=%f,SE2=%f\n",TR,KA,TA,TB,TC,Vrmax,Vrmin,KE,TE,KF,TF,SWITCH,E1,SE1,E2,SE2);

  // Convert exciter parameters from machine base to MVA base
  /*p_H *= mbase/sbase;
  p_D *= mbase/sbase;
  p_Xdp /= mbase/sbase;*/

}

/**
 * Saturation function
 * @ param x
 */
double Exdc1Exc::Sat(double x)
{
    B = log(SE2 / SE1)/(E2 - E1); //SJin: nan value when SE1 = 0 (read from EXDC1 input data)
    A = SE1 / exp(B * E1);
    return A * exp(B * x);
}

/**
 * Initialize exciter model before calculation
 * @param [output] values - array where initialized exciter variables should be set
 */
void Exdc1Exc::init(gridpack::ComplexType* values) 
{
  double Vm = sqrt(VD*VD + VQ*VQ); // SJin: voltage VD and VQ come from base_exc_model.hpp
  double mag = Vm;

  ///printf("exdc1: Efd = %f\n", Efd);
  x1 = Efd;
  //x4 = Efd * (KE + Sat(Efd));
  x4 = Efd*KE; // SJin: remove Sat function temporarily
  if (TB > (TS_THRESHOLD * ts)) // get ts from setTimestep() call 
    x3 = (x4 / KA) * (1 - TC / TB); // SJin: x4 is Vr 
  else
    x3 = x4 / KA;
  x2 = mag; // SJin: mag is Vterminal 
  //printf("KF = %f, TF = %f, x1 = %f, ts = %f\n", KF, TF, x1, ts);
  if (TF > (TS_THRESHOLD * ts)) 
    x5 = x1 * (KF / TF); // SJin: x1 is Ve
  else
    x5 = 0.0; // ignore this state
  Vref = mag + x4 / KA;
  //printf("Vref = %f\n", Vref);
  printf("exdc1 init:  %f\t%f\t%f\t%f\t%f\n", x1, x2, x3, x4, x5); 
  values[0] = x1;
  values[1] = x2;
  values[2] = x3;
  values[3] = x4;
  values[4] = x5;
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
}

double Exdc1Exc::getAngle(void)
{
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
    x1 = real(values[0]);
    x2 = real(values[1]);
    x3 = real(values[2]);
    x4 = real(values[3]);
    x5 = real(values[4]);
  } else if(p_mode == XDOTVECTOBUS) {
    dx1 = real(values[0]);
    dx2 = real(values[1]);
    dx3 = real(values[2]);
    dx4 = real(values[3]);
    dx5 = real(values[4]);
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
  // On fault (p_mode == FAULT_EVAL flag), the exciter variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    //values[delta_idx] = values[dw_idx] = 0.0;
    values[x1_idx] = 0.0;
    values[x2_idx] = 0.0;
    values[x3_idx] = 0.0;
    values[x4_idx] = 0.0;
    values[x5_idx] = 0.0;
  } else if(p_mode == RESIDUAL_EVAL) {
    // Exciter equations
    //printf("...........%f\t%f\t%f\t%f\t%f\n", x1, x2, x3, x4, x5);
    //printf(".........Vterminal = %f\n", Vterminal);
    bool flag2 = false; // set flags for residual function condisition to be used for Jacobian matrix calculation 
    bool flag3 = false;
    bool flag4 = false;
    bool flag5 = false;
    double Feedback;
    // State 2
    //printf("TR = %f\n", TR);
    if (TR > TS_THRESHOLD * t_inc) { 
      flag2 = true;
      values[x2_idx] = (Vterminal - x2) / TR - dx2; // SJin: x2 is Vsens 
    } else x2 = Vterminal; // propogate state immediately
    // State 5
    if (TF > TS_THRESHOLD * t_inc) {
      flag5 = true;
      values[x5_idx] = (x1 * KF / TF - x5) / TF - dx5;
      Feedback = x1 * KF / TF - x5;
    } else {
      x5 = 0;
      Feedback = 0;
    }
    //printf("TF = %f, t_inc = %f, x5 = %f, dx5 = %f\n", TF, t_inc, x5, dx5);
    // State 3
    double Vstab = 0.0; // SJin: Output from PSS, set to 0.0 for now.
    double LeadLagIN = Vref - x2 + Vstab - Feedback;
    //printf("Vref = %f, TB = %f, TC = %f\n", Vref, TB, TC); 
    double LeadLagOUT;
    //printf("LeadLagIN = %f, TC = %f, TB = %f, x3 = %f\n", LeadLagIN, TC, TB, x3);
    if (TB > (TS_THRESHOLD * t_inc)) {
      flag3 = true;
      values[x3_idx] = (LeadLagIN * (1 - TC / TB) - x3) / TB - dx3;
      LeadLagOUT = LeadLagIN * TC / TB + x3;
    } else 
      LeadLagOUT = LeadLagIN;
    // State 4
    //printf("x4 = %f, Vrmax = %f, Vrmin = %f, TA = %f, LeadLagOUT = %f, KA = %f\n", x4, Vrmax, Vrmin, TA, LeadLagOUT, KA);
    if (x4 > Vrmax) x4 = Vrmax;
    if (x4 < Vrmin) x4 = Vrmin;
    if (TA > (TS_THRESHOLD * t_inc)) { 
      flag4 = true;
      values[x4_idx] = (LeadLagOUT * KA - x4) / TA - dx4;
    } else {
      values[x4_idx] = -dx4;
      x4 = LeadLagOUT * KA;
    }
    if (dx4 > 0 && x4 >= Vrmax) { 
      values[x4_idx] = -dx4;
    } 
    if (dx4 < 0 && x4 <= Vrmin) { 
      values[x4_idx] = -dx4;
    }
    // State 1
    values[x1_idx] = x4 - x1 * (KE + Sat(x1)) - dx1;
      
    // Update Efd
    Efd = x1 * (1 + w);

  }
  
  return true;
}

/**
 * Return the exciter current injection (in rectangular form) 
 * @param [output] IGD - real part of the exciter current // SJin: match to Ir
 * @param [output] IGQ - imaginary part of the exciter current // SJin: match to Ii 
*/
void Exdc1Exc::getCurrent(double *IGD, double *IGQ)
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
    if (flag2) {
      row[idx] = 1; col[idx] = 0;
      values[idx] = 1/Vterminal*VD/TR; 
      idx++;
      row[idx] = 1; col[idx] = 1;
      values[idx] = 1/Vterminal*VQ/TR; 
      idx++;
    } else {
      if (flag3) {
        row[idx] = 2; col[idx] = 0;
        values[idx] = 1/Vterminal*VD/TR*(-(1-TC/TB)/TB); 
        idx++;
        row[idx] = 2; col[idx] = 1;
        values[idx] = 1/Vterminal*VQ/TR*(-(1-TC/TB)/TB); 
        idx++;
      } 
      if (flag4) {
        if (flag3) {
          row[idx] = 3; col[idx] = 0;
          values[idx] = 1/Vterminal*VD/TR*(-TC/TB*KA/TA); 
          idx++;
          row[idx] = 3; col[idx] = 1;
          values[idx] = 1/Vterminal*VQ/TR*(-TC/TB*KA/TA); 
          idx++;
        } else {
          row[idx] = 3; col[idx] = 0;
          values[idx] = 1/Vterminal*VD/TR*(-KA/TA); 
          idx++;
          row[idx] = 3; col[idx] = 1;
          values[idx] = 1/Vterminal*VQ/TR*(-KA/TA); 
          idx++;
        }
      }
    }   

    *nval = idx;
  } else if(p_mode == DIG_DX) { // SJin: Jacobian matrix block Jgx
    // These are the partial derivatives of the exciter currents (see getCurrent) w.r.t exciter variables

    *nval = idx;
  } else { // SJin: Jacobin matrix block Jfxi
    // Partials of exciter equations w.r.t exciter variables
    // Calculate directives of f1/xi
    row[idx] = 0; col[idx] = 0;
    values[idx] = -shift + (-KE-A*exp(B*x1)-A*B*x1*exp(x1));
    idx++;
    row[idx] = 0; col[idx] = 3;
    values[idx] = 1;
    idx++;

    // Calculate directives of f2/xi
    if (flag2) {
      row[idx] = 1; col[idx] = 1;
      values[idx] = -shift -1/TR;
      idx++;
    } else {
      row[idx] = 1; col[idx] = 1;
      values[idx] = -shift;
      idx++;
    }

    // Calculate directives of f3/xi
    if (flag3) {
      if (flag5) {
        row[idx] = 2; col[idx] = 0;
        values[idx] = -KF/TF*(1-TC/TB)/TB;
        idx++;
        if (flag2) {
          row[idx] = 2; col[idx] = 1;
          values[idx] = -(1-TC/TB)/TB;
          idx++;
        }
        row[idx] = 2; col[idx] = 2;
        values[idx] = -shift -1/TB;
        idx++;
        row[idx] = 2; col[idx] = 4;
        values[idx] = (1-TC/TB)/TB;
        idx++;
      } else {
        if (flag2) {
          row[idx] = 2; col[idx] = 1;
          values[idx] = -(1-TC/TB)/TB;
          idx++;
        }
        row[idx] = 2; col[idx] = 2;
        values[idx] = -shift -1/TB;
        idx++;
      } 
    } else {
      row[idx] = 2; col[idx] = 2;
      values[idx] = -shift;
      idx++;
    }

    // Calculate directives of f4/xi
    if (flag4) {
      if (flag3) {
        if (flag5) {
          row[idx] = 3; col[idx] = 0;
          values[idx] = -KF/TF*TC/TB*KA/TA;
          idx++;
          if (flag2) {
            row[idx] = 3; col[idx] = 1;
            values[idx] = -TC/TB*KA/TA;
            idx++;
          }
          row[idx] = 3; col[idx] = 2;
          values[idx] = KA/TA;
          idx++;
          row[idx] = 3; col[idx] = 3;
          values[idx] = -shift -1/TA;
          idx++;
          row[idx] = 3; col[idx] = 4;
          values[idx] = TC/TB*KA/TA;
          idx++;
        } else { // flag5 = false
          if (flag2) {
            row[idx] = 3; col[idx] = 1;
            values[idx] = -TC/TB*KA/TA;
            idx++;
          }
          row[idx] = 3; col[idx] = 2;
          values[idx] = KA/TA;
          idx++;
          row[idx] = 3; col[idx] = 3;
          values[idx] = -shift -1/TA;
          idx++;
        } // end of flag5
      } else { // flag3 = false
        if (flag5) {
          row[idx] = 3; col[idx] = 0;
          values[idx] = -KF/TF*KA/TA;
          idx++;
          if (flag2) {
            row[idx] = 3; col[idx] = 1;
            values[idx] = -KA/TA;
            idx++;
          }
          row[idx] = 3; col[idx] = 3;
          values[idx] = -shift -1/TA;
          idx++;
          row[idx] = 3; col[idx] = 4;
          values[idx] = KA/TA;
          idx++;
        } else { // flag5 = false
          if (flag2) {
            row[idx] = 3; col[idx] = 1;
            values[idx] = -KA/TA;
            idx++;
          }
          row[idx] = 3; col[idx] = 3;
          values[idx] = -1/TA;
          idx++;
        } // end of flag5 
      } // end of flag3      
    } else { // flag4 = false;
      row[idx] = 3; col[idx] = 3;
      values[idx] = -shift;
      idx++;
    } // end of flag4

    // Calculate directives of f5/xi
    if (flag5) {
      row[idx] = 4; col[idx] = 0;
      values[idx] = KF/TF/TF;
      idx++;
      row[idx] = 4; col[idx] = 4;
      values[idx] = -shift -1/TF;
      idx++;
    } else {
      row[idx] = 4; col[idx] = 4;
      values[idx] = -shift;
      idx++;
    }
 
    *nval = idx;
  }*/
  return true;
}

/**
 * Set the field voltage parameter inside the exciter
 * @param fldv value of the field voltage
 */
void Exdc1Exc::setFieldVoltage(double fldv)
{
  Efd = fldv;
  printf("Efd in EXDC1 = %f\n", Efd);
}

/**
 * Set the field current parameter inside the exciter
 * @param fldc value of the field current
 */
void Exdc1Exc::setFieldCurrent(double fldc)
{
  LadIfd = fldc;
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double Exdc1Exc::getFieldVoltage()
{
  return Efd;
}

/** 
 * Get the value of the field current parameter
 * @return value of field current
 */
double Exdc1Exc::getFieldCurrent()
{
  return 0.0;
}

/** 
 * Set the value of the Vterminal
 * @return value of field current
 */
void Exdc1Exc::setVterminal(double mag)
{
  Vterminal = mag;
}

/** 
 * Set the value of the omega
 * @return value of field current
 */
void Exdc1Exc::setOmega(double omega)
{
  w = omega;
}

void Exdc1Exc::setVcomp(double Vcomp)
{
}
   
/**
 * Set the value of the time step
 * @return value of the time step
 */
/*void Exdc1Exc::setTimestep(double timestep)
{
}*/

/**
 * Set the value of the time increment 
 * @return value of the time increment
 */
/*void Exdc1Exc::setTimeincrement(double timeincrement)
{
}*/

