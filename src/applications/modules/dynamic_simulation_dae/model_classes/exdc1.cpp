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

  nxexc = 5;
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
  *nvar = nxexc;
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
 * Set Jacobian block
 * @param values a 2-d array of Jacobian block for the bus
 */
bool Exdc1Exc::setJacobian(gridpack::ComplexType **values)
{
  int Vmeas_idx = offsetb;
  int xLL_idx   = offsetb+1;
  int VR_idx    = offsetb+2;
  int Efd_idx   = offsetb+3;
  int xf_idx    = offsetb+4;
  int VD_idx    = 0;
  int VQ_idx    = 1;
  double dVf_dxf  = 1.0;
  double dVf_dEfd = KF/TF;
  double dyLL_dVmeas=0.0,dyLL_dxLL=0.0;
  double dyLL_dVR=0.0,dyLL_dEfd=0.0;
  double dyLL_dxf=0.0;
  double Ec,yLL,Vf;

  Ec = sqrt(VD*VD + VQ*VQ);

  double dEc_dVD = VD/Ec;
  double dEc_dVQ = VQ/Ec;

  Vf = xf + KF/TF*Efd;

  if(p_mode == FAULT_EVAL) {
    // dEq.1_dX
    if(iseq_diff[0]) {
      values[Vmeas_idx][Vmeas_idx] = 1.0;
    } else {
      values[Vmeas_idx][Vmeas_idx] = -1.0;
      values[VD_idx][Vmeas_idx] = dEc_dVD;
      values[VQ_idx][Vmeas_idx]  = dEc_dVQ;
    }

    // dEq.2_dX
    if(iseq_diff[1]) {
      values[xLL_idx][xLL_idx] = 1.0;
      yLL = xLL + TC/TB*(Vref - Vmeas - Vf);
      dyLL_dxLL  = 1.0;
      dyLL_dVmeas = -TC/TB;
      dyLL_dxf    = -TC/TB*dVf_dxf;
    } else {
      values[Vmeas_idx][xLL_idx] = -1.0;
      values[xLL_idx][xLL_idx]  = -1.0;
      values[xf_idx][xLL_idx]    = -dVf_dxf;
      yLL = xLL;
      dyLL_dxLL = 1.0;
    }

    if(iseq_diff[2]) {
      values[VR_idx][VR_idx] = 1.0;
    } else {
      values[Vmeas_idx][VR_idx] = KA*dyLL_dVmeas;
      values[xLL_idx][VR_idx]   = KA*dyLL_dxLL;
      values[VR_idx][VR_idx]    = -1.0 + KA*dyLL_dVR;
      values[Efd_idx][VR_idx]   = KA*dyLL_dEfd;
      values[xf_idx][VR_idx]    = KA*dyLL_dxf;
    }
    
    values[Efd_idx][Efd_idx] = 1.0;
    values[xf_idx][xf_idx]   = 1.0;
  } else {
    // dEq.1_dX
    if(iseq_diff[0]) {
      values[Vmeas_idx][Vmeas_idx] = -1.0/TR - shift;
      values[VD_idx][Vmeas_idx]    = dEc_dVD/TR;
      values[VQ_idx][Vmeas_idx]    = dEc_dVQ/TR;
    } else {
      values[Vmeas_idx][Vmeas_idx] = -1.0;
      values[VD_idx][Vmeas_idx] = dEc_dVD;
      values[VQ_idx][Vmeas_idx]  = dEc_dVQ;
    }

    // dEq.2_dX
    if(iseq_diff[1]) {
      values[Vmeas_idx][xLL_idx] = ((1 - TC/TB)*-1.0)/TB;
      values[xLL_idx][xLL_idx] = -1.0/TB - shift;
      values[Efd_idx][xLL_idx] = ((1 - TC/TB)*-KF/TF)/TB;
      values[xf_idx][xLL_idx]  = (1 - TC/TB)/TB;

      yLL = xLL + TC/TB*(Vref - Vmeas - Vf);
      dyLL_dxLL  = 1.0;
      dyLL_dVmeas = -TC/TB;
      dyLL_dxf    = -TC/TB*dVf_dxf;
    } else {
      values[Vmeas_idx][xLL_idx] = -1.0;
      values[xLL_idx][xLL_idx]  = -1.0;
      values[xf_idx][xLL_idx]    = -dVf_dxf;
      yLL = xLL;
      dyLL_dxLL = 1.0;
    }

    if(iseq_diff[2]) {
      values[Vmeas_idx][VR_idx] = KA*dyLL_dVmeas/TA;
      values[xLL_idx][VR_idx]   = KA*dyLL_dxLL/TA;
      values[VR_idx][VR_idx]    = (-1.0 + KA*dyLL_dVR)/TA - shift;
      values[Efd_idx][VR_idx]   = KA*dyLL_dEfd/TA;
      values[xf_idx][VR_idx]    = KA*dyLL_dxf/TA;
    } else {
      values[Vmeas_idx][VR_idx] = KA*dyLL_dVmeas;
      values[xLL_idx][VR_idx]   = KA*dyLL_dxLL;
      values[VR_idx][VR_idx]    = -1.0 + KA*dyLL_dVR;
      values[Efd_idx][VR_idx]   = KA*dyLL_dEfd;
      values[xf_idx][VR_idx]    = KA*dyLL_dxf;
    }
    
    values[VR_idx][Efd_idx]  = 1/TE;
    values[Efd_idx][Efd_idx] = -KE/TE - shift;

    values[Efd_idx][xf_idx]  = -KF/(TF*TF);
    values[xf_idx][xf_idx]   = -1/TF - shift;
  }

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

bool Exdc1Exc::getFieldVoltagePartialDerivatives(int *xexc_loc,double *dEfd_dxexc,double *dEfd_dxgen)
{
  int i,nxgen;

  for(i=0; i < nxexc; i++) xexc_loc[i] = offsetb+i;

  for(i=0; i < nxexc; i++) dEfd_dxexc[i] = 0.0;
  dEfd_dxexc[3] = 1.0;

  getGenerator()->vectorSize(&nxgen);
  for(i=0; i < nxgen; i++) dEfd_dxgen[i] = 0.0;

  return true;
}

