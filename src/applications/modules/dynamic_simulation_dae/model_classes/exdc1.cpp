/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   exdc1.cpp
 * @author Shrirang Abhyankar 
 * @Last modified:   04/24/20
 *  
 * @brief EXDC1 exciter mdoel implementation 
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
  VRmax = 0.0; 
  VRmin = 0.0; 
  TC = 0.0; 
  TB = 0.0;
  KA = 0.0;
  TA = 0.0;
  KF = 0.0;
  TF = 0.0;
  SWITCH = 0;
  Efdthresh = 999;
  satA = satB = 0.0;
  
  VR_at_min = VR_at_max = false;

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
  BaseExcModel::load(data,idx); // load parameters in base exciter model
  if (!data->getValue(EXCITER_TR, &TR, idx)) TR = 0.0; // TR
  if (!data->getValue(EXCITER_KA, &KA, idx)) KA = 0.0; // KA
  if (!data->getValue(EXCITER_TA, &TA, idx)) TA = 0.0; // TA
  if (!data->getValue(EXCITER_TB, &TB, idx)) TB = 0.0; // TB
  if (!data->getValue(EXCITER_TC, &TC, idx)) TC = 0.0; // TC
  if (!data->getValue(EXCITER_VRMAX, &VRmax, idx)) VRmax = 0.0; // VRmax
  if (!data->getValue(EXCITER_VRMIN, &VRmin, idx)) VRmin = 0.0; // VRmin
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

  // Calculate saturation curve parameters
  if(SE2 != 0 && E2 != 0.0) {
    double alpha = std::sqrt(SE1*E1/SE2*E2);
    satA = (alpha*E2 - E1)/(alpha - 1.0);
    satB = SE1*E1/((E1 - satA)*(E1 - satA));
    Efdthresh = satA;
  }
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
  double SE=0.0;

  // Efd is already set by the generator model 
  Vmeas    = Ec;
  xf       = -KF/TF*Efd;

  if(Efd > Efdthresh) SE = satB*(Efd - satA)*(Efd - satA)/Efd;
  VR = (SE + KE)*Efd;

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
  double Ec,yLL,Vf,SE=0.0;
  Ec = sqrt(VD*VD + VQ*VQ);
  
  // On fault (p_mode == FAULT_EVAL flag), the exciter variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    // Vmeas equation
    if(iseq_diff[0]) values[0] = Vmeas - Vmeasprev;
    else values[0] = -Vmeas + Ec;


    Vf = xf + KF/TF*Efd;

    // xLL equation
    if(iseq_diff[1]) {
      values[1] = xLL - xLLprev;
      yLL = xLL + TC/TB*(Vref - Vmeas - Vf);
    } else {
      values[1] = -xLL + Vref - Vmeas - Vf;
      yLL = xLL;
    }

    // VR equation
    if(VR_at_min) {
      values[2] = VR - VRmin;
    } else if(VR_at_max) {
      values[2] = VR - VRmax;
    } else {
      if(iseq_diff[2]) values[2] = VR - VRprev;
      else values[2] = -VR + KA*yLL;
    }
    // Efd equation
    values[3] = Efd - Efdprev;

    // xf equation
    values[4] = xf - xfprev;

  } else if(p_mode == RESIDUAL_EVAL) {
    // Vmeas equation
    if(iseq_diff[0]) values[0] = (-Vmeas + Ec)/TR - dVmeas;
    else values[0] = -Vmeas + Ec;

    // xLL equation
    Vf = xf + KF/TF*Efd;
    if(iseq_diff[1]) {
      values[1] = (-xLL + (1 - TC/TB)*(Vref - Vmeas - Vf))/TB - dxLL;
      yLL = xLL + TC/TB*(Vref - Vmeas - Vf);
    } else {
      values[1] = -xLL + Vref - Vmeas - Vf;
      yLL = xLL;
    }

    // VR equation
    if(VR_at_min) {
      values[2] = VR - VRmin;
    } else if(VR_at_max) {
      values[2] = VR - VRmax;
    } else {
      if(iseq_diff[2]) values[2] = (-VR + KA*yLL)/TA - dVR;
      else values[2] = -VR + KA*yLL;
    }

    if(Efd > Efdthresh) SE = satB*(Efd - satA)*(Efd - satA)/Efd;
    values[3] = (VR - (SE + KE)*Efd)/TE - dEfd;

    // xf equation
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
  double Ec,yLL,Vf,SE=0.0;

  Ec = sqrt(VD*VD + VQ*VQ);

  double dEc_dVD = VD/Ec;
  double dEc_dVQ = VQ/Ec;

  Vf = xf + KF/TF*Efd;

  if(p_mode == FAULT_EVAL) {
    // Partial derivatives of Vmeas equation
    if(iseq_diff[0]) {
      values[Vmeas_idx][Vmeas_idx] = 1.0;
    } else {
      values[Vmeas_idx][Vmeas_idx] = -1.0;
      values[VD_idx][Vmeas_idx] = dEc_dVD;
      values[VQ_idx][Vmeas_idx]  = dEc_dVQ;
    }

    // Partial derivatives of xLL equation
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

    // Partial derivatives of VR equation
    if(VR_at_min || VR_at_max) {
      values[VR_idx][VR_idx] = 1.0;
    } else {
      if(iseq_diff[2]) {
	values[VR_idx][VR_idx] = 1.0;
      } else {
	values[Vmeas_idx][VR_idx] = KA*dyLL_dVmeas;
	values[xLL_idx][VR_idx]   = KA*dyLL_dxLL;
	values[VR_idx][VR_idx]    = -1.0 + KA*dyLL_dVR;
	values[Efd_idx][VR_idx]   = KA*dyLL_dEfd;
	values[xf_idx][VR_idx]    = KA*dyLL_dxf;
      }
    }
    
    // Partial derivatives of Efd equation
    values[Efd_idx][Efd_idx] = 1.0;

    // Partial derivatives of xf equation
    values[xf_idx][xf_idx]   = 1.0;
  } else {
    // Partial derivatives of Vmeas equation
    if(iseq_diff[0]) {
      values[Vmeas_idx][Vmeas_idx] = -1.0/TR - shift;
      values[VD_idx][Vmeas_idx]    = dEc_dVD/TR;
      values[VQ_idx][Vmeas_idx]    = dEc_dVQ/TR;
    } else {
      values[Vmeas_idx][Vmeas_idx] = -1.0;
      values[VD_idx][Vmeas_idx] = dEc_dVD;
      values[VQ_idx][Vmeas_idx]  = dEc_dVQ;
    }

    // Partial derivatives of xLL equation
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

    // Partial derivatives of VR equation
    if(VR_at_min || VR_at_max) {
      values[VR_idx][VR_idx] = 1.0;
    } else {
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
    }
    
    // Partial derivatives of Efd equation
    double dSE_dEfd = 0.0;
    if(Efd > Efdthresh) dSE_dEfd = 2*satB*(Efd - satA);
    values[VR_idx][Efd_idx]  = 1/TE;
    values[Efd_idx][Efd_idx] = -(dSE_dEfd + KE)/TE - shift;

    // Partial derivatives of xf equation
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

/**
 * Update the event function values
 */
void Exdc1Exc::eventFunction(const double&t,gridpack::ComplexType *state,std::vector<std::complex<double> >& evalues)
{
  int offset    = getLocalOffset();
  int Vmeas_idx = offset;
  int xLL_idx   = offset+1;
  int VR_idx    = offset+2;
  int Efd_idx   = offset+3;
  int xf_idx    = offset+4;

  Vmeas = real(state[Vmeas_idx]);
  xLL   = real(state[xLL_idx]);
  VR    = real(state[VR_idx]);
  Efd   = real(state[Efd_idx]);
  xf    = real(state[xf_idx]);

  /* Only considering limits on VR */
  double Vf,yLL,dVR_dt;
  Vf = xf + KF/TF*Efd;
  if(iseq_diff[1]) {
    yLL = xLL + TC/TB*(Vref - Vmeas - Vf);
  } else {
    yLL = xLL;
  }
  dVR_dt = (-VR + KA*yLL)/TA;

  /* Limits on VR */
  if(!VR_at_min) {
    evalues[0] = VR - VRmin;
  } else {
    evalues[0] = -dVR_dt; /* Release when derivative reaches 0 */
  }

  if(!VR_at_max) {
    evalues[1] = VRmax - VR;
    //    printf("VR = %f\n",VR);
  } else {
    evalues[1] = dVR_dt; /* Release when derivative reaches 0 */
    //    printf("VR = %f, dVR_dt = %f\n",VR,dVR_dt);
  }
} 

/**
 * Reset limiter flags after a network resolve
 */
void Exdc1Exc::resetEventFlags()
{
  /* Note that the states are already pushed onto the network, so we can access these
     directly
  */
  double Vf,yLL,dVR_dt;
  Vf = xf + KF/TF*Efd;
  if(iseq_diff[1]) {
    yLL = xLL + TC/TB*(Vref - Vmeas - Vf);
  } else {
    yLL = xLL;
  }
  dVR_dt = (-VR + KA*yLL)/TA;

  if(!VR_at_min) {
    if(VR - VRmin < 0) VR_at_min = true;
  } else {
    if(dVR_dt > 0) VR_at_min = false; /* Release */
  }

  if(!VR_at_max) {
    if(VRmax - VR < 0) VR_at_max = true;
  } else {
    if(dVR_dt < 0) VR_at_max = false; /* Release */
  }
}

/**
 * Event handler
 */
void Exdc1Exc::eventHandlerFunction(const bool *triggered, const double& t, gridpack::ComplexType *state)
{
  int offset    = getLocalOffset();
  int Vmeas_idx = offset;
  int xLL_idx   = offset+1;
  int VR_idx    = offset+2;
  int Efd_idx   = offset+3;
  int xf_idx    = offset+4;

  Vmeas = real(state[Vmeas_idx]);
  xLL   = real(state[xLL_idx]);
  VR    = real(state[VR_idx]);
  Efd   = real(state[Efd_idx]);
  xf    = real(state[xf_idx]);

  /* Only considering limits on VR */
  double Vf,yLL,dVR_dt;
  Vf = xf + KF/TF*Efd;
  if(iseq_diff[1]) {
    yLL = xLL + TC/TB*(Vref - Vmeas - Vf);
  } else {
    yLL = xLL;
  }
  dVR_dt = (-VR + KA*yLL)/TA;

  if(triggered[0]) {
    if(!VR_at_min && dVR_dt < 0) {
      /* Hold VR at VRmin */
      VR_at_min = true;
    } else {
      /* Release */
      VR_at_min = false;
    }
  }

  if(triggered[1]) {
    if(!VR_at_max && dVR_dt > 0) {
      /* Hold VR at Vamax */
      VR_at_max = true;
    } else {
      /* Release */
      VR_at_max = false;
    }
  }
}

/**
 * Set event
 */
void Exdc1Exc::setEvent(gridpack::math::DAESolver::EventManagerPtr eman)
{
  gridpack::math::DAESolver::EventPtr e(new Exdc1ExcEvent(this));

  eman->add(e);
}

void Exdc1ExcEvent::p_update(const double& t,gridpack::ComplexType *state)
{
  p_exc->eventFunction(t,state,p_current);
}

void Exdc1ExcEvent::p_handle(const bool *triggered, const double& t, gridpack::ComplexType *state)
{
  p_exc->eventHandlerFunction(triggered,t,state);
}
