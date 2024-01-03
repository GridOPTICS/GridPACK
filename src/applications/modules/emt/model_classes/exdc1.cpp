/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   exdc1.cpp
 *  
 * @brief EXDC1 exciter mdoel implementation 
 *
 *
 */

#include <exdc1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Exdc1::Exdc1(void)
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

Exdc1::~Exdc1(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to BaseExciterModel
 */
void Exdc1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTExcModel::load(data,idx); // load parameters in base exciter model
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
  //  iseq_diff[0] = (TR == 0)?0:1;
  //  iseq_diff[1] = (TB == 0 || TC == 0)?0:1;
  //  iseq_diff[2] = (TA == 0)?0:1;
  //  iseq_diff[3] = 1; // TE is always > 0
  //  iseq_diff[4] = 1; //TF always > 0

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
void Exdc1::init(gridpack::ComplexType* xin) 
{
  gridpack::ComplexType *x = xin+offsetb; // exciter array starts from this location
  double Ec = sqrt(VD*VD + VQ*VQ);
  double yLL;
  double Vf=0.0;
  double SE=0.0;

  // Efd is already set by the generator model
  Efd = Efd0;
  
  Vmeas    = Ec;
  xf       = -KF/TF*Efd;

  if(Efd > Efdthresh) SE = satB*(Efd - satA)*(Efd - satA)/Efd;
  VR = (SE + KE)*Efd;

  yLL     = VR/KA;

  Vf   = xf + KF/TF*Efd;
  Vref = yLL + Vmeas + Vf;

  if(TB != 0 && TC != 0) xLL = (1 - TC/TB)*(Vref - Vmeas - Vf);
  else xLL = Vref - Vmeas - Vf;

  x[0] = Vmeas;
  x[1] = xLL;
  x[2] = VR;
  x[3] = Efd;
  x[4] = xf;
}

/**
 * Write output from exciters to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Exdc1::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out exciter state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Exdc1::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto exciters
 * @param values array containing exciter state variables
*/
void Exdc1::setValues(gridpack::ComplexType *val)
{
  gridpack::ComplexType *values = val+offsetb; // exciter array starts from this location

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
  }
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Exdc1::vectorGetValues(gridpack::ComplexType *values)
{
  gridpack::ComplexType *f = values+offsetb; // exciter array starts from this location

  double Ec,yLL,Vf,SE=0.0;

  double vabc[3],vdq0[3];

  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;

  double delta = getGenerator()->getAngle();
  
  abc2dq0(vabc,p_time,delta,vdq0);
  double Vd, Vq;
  Vd = vdq0[0]; Vq = vdq0[1];
  
  Ec = sqrt(Vd*Vd + Vq*Vq);
  
  if(p_mode == RESIDUAL_EVAL) {
    // Vmeas equation
    if(TR != 0) f[0] = (-Vmeas + Ec)/TR - dVmeas;
    else f[0] = -Vmeas + Ec;

    // xLL equation
    Vf = xf + KF/TF*Efd;
    if(TB != 0 && TC != 0) {
      f[1] = (-xLL + (1 - TC/TB)*(Vref - Vmeas - Vf))/TB - dxLL;
      yLL = xLL + TC/TB*(Vref - Vmeas - Vf);
    } else {
      f[1] = -xLL + Vref - Vmeas - Vf;
      yLL = xLL;
    }

    // VR equation
    if(VR_at_min) {
      f[2] = VR - VRmin;
    } else if(VR_at_max) {
      f[2] = VR - VRmax;
    } else {
      if(TA != 0) f[2] = (-VR + KA*yLL)/TA - dVR;
      else f[2] = -VR + KA*yLL;
    }

    if(Efd > Efdthresh) SE = satB*(Efd - satA)*(Efd - satA)/Efd;
    f[3] = (VR - (SE + KE)*Efd)/TE - dEfd;

    // xf equation
    f[4] = (-xf - KF/TF*Efd)/TF - dxf;    
  }
}

/**
 * Set Jacobian block
 * @param values a 2-d array of Jacobian block for the bus
 */
bool Exdc1::setJacobian(gridpack::ComplexType **values)
{
  return false;
}

/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values
 */
int Exdc1::matrixNumValues()
{
  return 0;
}

  /**
   * Return values from Jacobian matrix
   * @param nvals: number of values to be inserted
   * @param values: pointer to matrix block values
   * @param rows: pointer to matrix block rows
   * @param cols: pointer to matrix block cols
   */
void Exdc1::matrixGetValues(int *nvals, gridpack::ComplexType *values, int *rows, int *cols)
{

}


/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double Exdc1::getFieldVoltage()
{
  return Efd;
}

bool Exdc1::getFieldVoltagePartialDerivatives(int *xexc_loc,double *dEfd_dxexc,double *dEfd_dxgen)
{
  return false;
}

/**
 * Update the event function values
 */
void Exdc1::eventFunction(const double&t,gridpack::ComplexType *state,std::vector<std::complex<double> >& evalues)
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
  if(TB != 0 && TC != 0) {
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
void Exdc1::resetEventFlags()
{
  /* Note that the states are already pushed onto the network, so we can access these
     directly
  */
  double Vf,yLL,dVR_dt;
  Vf = xf + KF/TF*Efd;
  if(TB != 0 && TC != 0) {
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
void Exdc1::eventHandlerFunction(const bool *triggered, const double& t, gridpack::ComplexType *state)
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
  if(TB != 0 && TC != 0) {
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
void Exdc1::setEvent(gridpack::math::DAESolver::EventManagerPtr eman)
{
  gridpack::math::DAESolver::EventPtr e(new Exdc1Event(this));

  eman->add(e);
}

void Exdc1Event::p_update(const double& t,gridpack::ComplexType *state)
{
  p_exc->eventFunction(t,state,p_current);
}

void Exdc1Event::p_handle(const bool *triggered, const double& t, gridpack::ComplexType *state)
{
  p_exc->eventHandlerFunction(triggered,t,state);
}
