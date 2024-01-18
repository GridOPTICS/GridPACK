/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   sexs.cpp
 * @author Shrirang Abhyankar
 * @Added:   Nov 6, 2022
 * @Modified by Shuangshuang Jin: Jan 18, 2024
 * 
 * @brief  
 * 
 * 
 */

#include <sexs.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Sexs::Sexs(void)
{
  x1 = 0.0;
  dx1 = 0.0;
  x2 = 0.0;
  dx2 = 0.0;
  TA_OVER_TB = 0.0;
  TB = 0.0;
  K = 0.0;
  TE = 0.0;
  EMIN = 0.0;
  EMAX = 0.0;
  nxexc = 2;
}

Sexs::~Sexs(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to BaseExciterModel
 */
void Sexs::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTExcModel::load(data,idx); // load parameters in base sex model
  if (!data->getValue(EXCITER_TA_OVER_TB, &TA_OVER_TB, idx)) TA_OVER_TB = 0.0; // TA_OVER_TB
  if (!data->getValue(EXCITER_TB, &TB, idx)) TB = 0.0; // TB
  if (!data->getValue(EXCITER_K, &K, idx))   K  = 0.0; // K
  if (!data->getValue(EXCITER_TE, &TE, idx)) TE = 0.0; // TE
  if (!data->getValue(EXCITER_EMIN, &EMIN, idx)) EMIN = 0.0; // EMIN
  if (!data->getValue(EXCITER_EMAX, &EMAX, idx)) EMAX = 0.0; // EMAX
  
  TA = TA_OVER_TB*TB;

  // Set parameters for the first block
  /*leadlagblock.setparams(TA,TB);

  zero_TE = false;
  if(fabs(TE) < 1e-6) {
    zero_TE = true;
  }
  
  // Set parameters for the second block
  if(!zero_TE) {
    filterblock.setparams(K,TE,EMIN,EMAX,-1000.0,1000);
  } else {
    gainblock.setparams(K,EMIN,EMAX);
  }*/

  Vs = 0.0; 
}

/**
 * Initialize exciter model before calculation
 * @param [output] values - array where initialized exciter variables should be set
 */
void Sexs::init(gridpack::ComplexType *values)
{
  //gridpack::ComplexType *x = xin+offsetb; // exciter array starts from this location
  
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
/*void gridpack::dynamic_simulation::SexsModel::init(double Vm, double ang, double ts)
{
  double y1,u1;

  // For initialization, we are given the output for the model Efd.
  // To initialize the model blocks, we need to go backwards starting
  // from initializing second block and then the first one and then
  // calculating the model input Vref
  
  // Initialize second block
  if(!zero_TE) {
    y1 = filterblock.init_given_y(Efd);
  } else {
    y1 = std::min(EMAX,std::max(EMIN,Efd/K));
  }

  // Initialize first block
  u1 = leadlagblock.init_given_y(y1); 

  // Note: Vm is same as Ec
  Vref = Vm + u1 - Vs;   // Voltage reference initial value

}*/

/**
 * Write output from exciters to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Sexs::serialWrite(char *string, const int bufsize, const char *signal)
{
    return false;
}

/**
 * Write out exciter state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Sexs::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto exciters
 * @param values array containing exciter state variables
*/
void Sexs::setValues(gridpack::ComplexType *val)
{
  gridpack::ComplexType *values = val+offsetb; // exciter array starts from this location

  if(p_mode == XVECTOBUS) {
    x1 = real(values[0]);
    x2 = real(values[1]);
  } else if(p_mode == XDOTVECTOBUS) {
    dx1 = real(values[0]);
    dx2 = real(values[1]);
  }
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Sexs::vectorGetValues(gridpack::ComplexType *values)
{
  gridpack::ComplexType *f = values+offsetb; // exciter array starts from this location

  /*double Ec,yLL,Vf,SE=0.0;

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
  }*/
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
/*void gridpack::dynamic_simulation::SexsModel::predictor(double t_inc, bool flag)
{
  double u1,y1;

  u1 = Vref + Vs - Ec;

  // Calculate first block output, last input flag = true
  // tells the block to do the state update (predictor update)
  y1 = leadlagblock.getoutput(u1,t_inc,PREDICTOR,true);

  // Calculate second block output, last input flag = true
  // tells the block to do the state update (predictor update)
  if(!zero_TE) {
    Efd = filterblock.getoutput(y1,t_inc,PREDICTOR,true);
  } else {
    Efd = gainblock.getoutput(y1);
  }
}*/

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
/*oid gridpack::dynamic_simulation::SexsModel::corrector(double t_inc, bool flag)
{
  double u1,y1;

  u1 = Vref + Vs - Ec;

  // Calculate first block output, last input flag = true
  // tells the block to do the state update (corrector update)
  y1 = leadlagblock.getoutput(u1,t_inc,CORRECTOR,true);

  // Calculate second block output, last input flag = true
  // tells the block to do the state update (corrector update)
  if(!zero_TE) {
    Efd = filterblock.getoutput(y1,t_inc,CORRECTOR,true);
  } else {
    Efd = gainblock.getoutput(y1);
  }
}*/

/**
 * Set Jacobian block
 * @param values a 2-d array of Jacobian block for the bus
 */
bool Sexs::setJacobian(gridpack::ComplexType **values)
{
  return false;
}

/**
   Non-zero pattern of the Jacobian (x denotes non-zero entry)
         Vmeas    xLL    VR    Efd    xf    delta    va    vb    vc
 eq.0 |    x                                  x      x     x     x
 eq.1 |    x      x             x     x
 eq.2 |    x      x      x      x     x
 eq.3 |                  x      x
 eq.4 |                         x     x

 Number of non-zeros = 5 + 4 + 5 + 2 + 2 = 18 
 * Get number of matrix values contributed by exciter
 * @return number of matrix values
 */
int Sexs::matrixNumValues()
{
  return 18;
}

/**
   * Return values from Jacobian matrix
   * @param nvals: number of values to be inserted
   * @param values: pointer to matrix block values
   * @param rows: pointer to matrix block rows
   * @param cols: pointer to matrix block cols
   */
void Sexs::matrixGetValues(int *nvals, gridpack::ComplexType *values, int *rows, int *cols)
{
  /*int ctr = 0;

  int Vmeas_idx = p_gloc;
  int xLL_idx   = p_gloc + 1;
  int VR_idx    = p_gloc + 2;
  int Efd_idx   = p_gloc + 3;
  int xf_idx    = p_gloc + 4;
  int delta_idx;
  int va_idx    = p_glocvoltage;
  int vb_idx    = p_glocvoltage+1;
  int vc_idx    = p_glocvoltage+2;

  double delta = getGenerator()->getAngle(&delta_idx);
  double Tdq0[3][3];
  double dTdq0ddelta[3][3];

  getTdq0(p_time,delta,Tdq0);
  getdTdq0dtheta(p_time,delta,dTdq0ddelta);

  double Ec, Vd, Vq, V0;
  double vabc[3],vdq0[3];
  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;

  abc2dq0(vabc,p_time,delta,vdq0);
  Vd = vdq0[0]; Vq = vdq0[1]; V0 = vdq0[2];

  Ec = sqrt(Vd*Vd + Vq*Vq);

  double dEc_dVd,dEc_dVq;
  dEc_dVd = Vd/Ec; dEc_dVq = Vq/Ec;

  double dVd_dvabc[3], dVq_dvabc[3];
  dVd_dvabc[0] = Tdq0[0][0]; dVd_dvabc[1] = Tdq0[0][1]; dVd_dvabc[2] = Tdq0[0][2];
  dVq_dvabc[0] = Tdq0[1][0]; dVq_dvabc[1] = Tdq0[1][1]; dVq_dvabc[2] = Tdq0[1][2];

  double dVd_ddelta = dTdq0ddelta[0][0]*vabc[0] + dTdq0ddelta[0][1]*vabc[1] + dTdq0ddelta[0][2]*vabc[2];
  double dVq_ddelta = dTdq0ddelta[1][0]*vabc[0] + dTdq0ddelta[1][1]*vabc[1] + dTdq0ddelta[1][2]*vabc[2];
  

  rows[ctr] = Vmeas_idx;  cols[ctr] = Vmeas_idx;
  rows[ctr+1] = Vmeas_idx; cols[ctr+1] = delta_idx;
  rows[ctr+2] = Vmeas_idx; cols[ctr+2] = va_idx;
  rows[ctr+3] = Vmeas_idx; cols[ctr+3] = vb_idx;
  rows[ctr+4] = Vmeas_idx; cols[ctr+4] = vc_idx;
  
  if(TR != 0) {
    values[ctr]   = -1.0/TR - shift;
    values[ctr+1] = (dEc_dVd*dVd_ddelta + dEc_dVq*dVq_ddelta)/TR;
    values[ctr+2] = (dEc_dVd*dVd_dvabc[0] + dEc_dVq*dVq_dvabc[0])/TR;
    values[ctr+3] = (dEc_dVd*dVd_dvabc[1] + dEc_dVq*dVq_dvabc[1])/TR;
    values[ctr+4] = (dEc_dVd*dVd_dvabc[2] + dEc_dVq*dVq_dvabc[2])/TR;
  } else {
    values[ctr]   = -1.0;
    values[ctr+1] = (dEc_dVd*dVd_ddelta + dEc_dVq*dVq_ddelta);
    values[ctr+2] = (dEc_dVd*dVd_dvabc[0] + dEc_dVq*dVq_dvabc[0]);
    values[ctr+3] = (dEc_dVd*dVd_dvabc[1] + dEc_dVq*dVq_dvabc[1]);
    values[ctr+4] = (dEc_dVd*dVd_dvabc[2] + dEc_dVq*dVq_dvabc[2]);
  }

  ctr += 5;

  rows[ctr] = xLL_idx; cols[ctr] = Vmeas_idx;
  rows[ctr+1] = xLL_idx; cols[ctr+1] = xLL_idx;
  rows[ctr+2] = xLL_idx; cols[ctr+2] = Efd_idx;
  rows[ctr+3] = xLL_idx; cols[ctr+3] = xf_idx;

  double dVf_dxf = 1.0;
  double dVf_dEfd = KF/TF;
  double param = (1 - TC/TB);
  double dyLL_dVmeas = 0.0, dyLL_dxLL = 0.0;
  double dyLL_dEfd = 0.0, dyLL_dxf = 0.0;
  if(TB != 0 && TC != 0) {
    values[ctr] = (param*-1)/TB;
    values[ctr+1] = -1.0/TB - shift;
    values[ctr+2] = (param*(-dVf_dEfd))/TB;
    values[ctr+3] = (param*(-dVf_dxf))/TB;

    dyLL_dVmeas = TC/TB*-1.0;
    dyLL_dxLL = 1.0;
    dyLL_dEfd = TC/TB*-KF/TF;
    dyLL_dxf  = TC/TB*-1.0;
  } else {
    values[ctr] = -1.0;
    values[ctr+1] = -1.0;
    values[ctr+2] = -dVf_dEfd;
    values[ctr+3] = -dVf_dxf;

    dyLL_dxLL = 1.0;
  }

  ctr += 4;

  rows[ctr] = VR_idx;  cols[ctr] = Vmeas_idx;
  rows[ctr+1] = VR_idx; cols[ctr+1] = xLL_idx;
  rows[ctr+2] = VR_idx; cols[ctr+2] = VR_idx;
  rows[ctr+3] = VR_idx; cols[ctr+3] = Efd_idx;
  rows[ctr+4] = VR_idx; cols[ctr+4] = xf_idx;

  values[ctr] = values[ctr+1] = values[ctr+2] = values[ctr+3] = values[ctr+4] = 0.0;

  if(VR_at_min || VR_at_max) {
    values[ctr+2] = 1.0;
  } else {
    if(TA != 0) {
      values[ctr] =   (KA*dyLL_dVmeas)/TA;
      values[ctr+1] = (KA*dyLL_dxLL)/TA;
      values[ctr+2] = -1.0/TA - shift;
      values[ctr+3] = (KA*dyLL_dEfd)/TA;
      values[ctr+4] = (KA*dyLL_dxf)/TA;
    } else {
      values[ctr] =   (KA*dyLL_dVmeas);
      values[ctr+1] = (KA*dyLL_dxLL);
      values[ctr+2] = -1.0/TA;
      values[ctr+3] = (KA*dyLL_dEfd)/TA;
      values[ctr+4] = (KA*dyLL_dxf)/TA;
    }
  }
      
  ctr += 5;

  double SE = 0.0;
  double dSE_dEfd = 0.0;
  if(Efd > Efdthresh) {
    SE = satB*(Efd - satA)*(Efd - satA)/Efd;
    dSE_dEfd = 2*satB*(Efd - satA)/Efd - SE/Efd;
  }
  rows[ctr] = Efd_idx; cols[ctr] = VR_idx;
  rows[ctr+1] = Efd_idx; cols[ctr+1] = Efd_idx;

  values[ctr] = (1.0)/TE;

  values[ctr+1] = (-SE - dSE_dEfd*Efd - KE)/TE - shift;
  
  ctr += 2;

  rows[ctr] = xf_idx; cols[ctr] = Efd_idx;
  rows[ctr+1] = xf_idx; cols[ctr+1] = xf_idx;

  values[ctr]   = (-KF/TF)/TF;
  values[ctr+1] = (-1.0)/TF - shift;

  ctr += 2;
  
  *nvals = ctr;*/
		     
}


/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double Sexs::getFieldVoltage()
{
  return Efd;
}

/** 
 * Get the value of the field voltage parameter
 * and its global location
 * @return value of field voltage
 */
double Sexs::getFieldVoltage(int *Efd_gloc)
{
  *Efd_gloc = p_gloc + 3;
  return Efd;
}

bool Sexs::getFieldVoltagePartialDerivatives(int *xexc_loc,double *dEfd_dxexc,double *dEfd_dxgen)
{
  return false;
}

/**
 * Update the event function values
 */
void Sexs::eventFunction(const double&t,gridpack::ComplexType *state,std::vector<std::complex<double> >& evalues)
{
  /*int offset    = getLocalOffset();
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

  // Only considering limits on VR 
  double Vf,yLL,dVR_dt;
  Vf = xf + KF/TF*Efd;
  if(TB != 0 && TC != 0) {
    yLL = xLL + TC/TB*(Vref - Vmeas - Vf);
  } else {
    yLL = xLL;
  }
  dVR_dt = (-VR + KA*yLL)/TA;

  // Limits on VR 
  if(!VR_at_min) {
    evalues[0] = VR - VRmin;
  } else {
    evalues[0] = -dVR_dt; // Release when derivative reaches 0 
  }

  if(!VR_at_max) {
    evalues[1] = VRmax - VR;
    //    printf("VR = %f\n",VR);
  } else {
    evalues[1] = dVR_dt; // Release when derivative reaches 0 
    //    printf("VR = %f, dVR_dt = %f\n",VR,dVR_dt);
  }*/
} 

/**
 * Reset limiter flags after a network resolve
 */
void Sexs::resetEventFlags()
{
  /* Note that the states are already pushed onto the network, so we can access these
     directly
  */
  /*double Vf,yLL,dVR_dt;
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
    if(dVR_dt > 0) VR_at_min = false; // Release 
  }

  if(!VR_at_max) {
    if(VRmax - VR < 0) VR_at_max = true;
  } else {
    if(dVR_dt < 0) VR_at_max = false; // Release 
  }*/
}

/**
 * Event handler
 */
void Sexs::eventHandlerFunction(const bool *triggered, const double& t, gridpack::ComplexType *state)
{
  /*int offset    = getLocalOffset();
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

  // Only considering limits on VR 
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
      // Hold VR at VRmin 
      VR_at_min = true;
    } else {
      // Release 
      VR_at_min = false;
    }
  }

  if(triggered[1]) {
    if(!VR_at_max && dVR_dt > 0) {
      // Hold VR at Vamax 
      VR_at_max = true;
    } else {
      // Release 
      VR_at_max = false;
    }
  }*/
}

/**
 * Set event
 */
void Sexs::setEvent(gridpack::math::DAESolver::EventManagerPtr eman)
{
  gridpack::math::DAESolver::EventPtr e(new SexsEvent(this));

  eman->add(e);
}

void SexsEvent::p_update(const double& t,gridpack::ComplexType *state)
{
  p_exc->eventFunction(t,state,p_current);
}

void SexsEvent::p_handle(const bool *triggered, const double& t, gridpack::ComplexType *state)
{
  p_exc->eventHandlerFunction(triggered,t,state);
}
