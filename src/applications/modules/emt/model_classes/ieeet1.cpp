/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ieeet1.cpp
 *  
 * @brief IEEET1 exciter mdoel implementation 
 *
 *
 */

#include <ieeet1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Ieeet1::Ieeet1(void)
{
  Vmeas = 0.0; 
  dVmeas = 0.0;
  VR  = 0.0;
  dVR = 0.0;
  Efd    = 0.0;
  dEfd   = 0.0;
  xf    = 0.0;
  dxf   = 0.0;
  TR = 0.0; 
  VRmax = 0.0; 
  VRmin = 0.0; 
  KA = 0.0;
  TA = 0.0;
  KF = 0.0;
  TF = 0.0;
  SWITCH = 0;
  Efdthresh = 999;
  satA = satB = 0.0;
  Vs = 0.0;
  
  VR_at_min = VR_at_max = false;

  has_Sat = true;
  zero_TA = false;
  zero_TR = false;

  nxexc = 4;
}

Ieeet1::~Ieeet1(void)
{
}

void Ieeet1::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxexc = 0;
  *nvar = nxexc;
}

void Ieeet1::preStep(double time, double timestep)
{
  if(integrationtype != EXPLICIT) return;

  double vabc[3],vdq0[3];

  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;

  double delta = getGenerator()->getAngle();
  
  abc2dq0(vabc,p_time,delta,vdq0);
  double Vd, Vq;
  Vd = vdq0[0]; Vq = vdq0[1];
  
  Ec = sqrt(Vd*Vd + Vq*Vq);
  
  if(!zero_TR) {
    Vmeas = Vmeas_blk.getoutput(Ec,timestep,true);
  } else {
    Vmeas = Ec;
  }

  double Verr = Vref - Vmeas + Vs;

  VF = Feedback_blk.getoutput(Efd,timestep,true);

  double Regulator_blk_in = Verr - VF;

  if(zero_TA) {
    VR = Regulator_gain_blk.getoutput(Regulator_blk_in);
  } else {
    VR = Regulator_blk.getoutput(Regulator_blk_in,timestep,true);
  }

  double SE = 0.0;
  if(has_Sat) {
    if(Efd > Efdthresh) SE = satB*(Efd - satA)*(Efd - satA)/Efd;
  }

  double output_blk_in = VR - (KE + SE)*Efd;

  Efd = Output_blk.getoutput(output_blk_in,timestep,true);
}

void Ieeet1::postStep(double time)
{

}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to BaseExciterModel
 */
void Ieeet1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTExcModel::load(data,idx); // load parameters in base exciter model
  if (!data->getValue(EXCITER_TR, &TR, idx)) TR = 0.0; // TR
  if (!data->getValue(EXCITER_KA, &KA, idx)) KA = 0.0; // KA
  if (!data->getValue(EXCITER_TA, &TA, idx)) TA = 0.0; // TA
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

  if(fabs(SE1*SE2) < 1e-6) has_Sat = false;
  if(fabs(TA) < 1e-6) zero_TA = true;
  if(fabs(TR) < 1e-6) zero_TR = true;

  if(has_Sat) {
    double alpha = std::sqrt(SE1*E1/(SE2*E2));
    satA = (alpha*E2 - E1)/(alpha - 1.0);
    satB = SE1*E1/((E1 - satA)*(E1 - satA));
    Efdthresh = satA;
  }

  if(integrationtype != IMPLICIT) {
      // Set up blocks
    if(!zero_TR) {
      Vmeas_blk.setparams(1.0,TR);
    }

    if(!zero_TA) {
      Regulator_blk.setparams(KA,TA,VRmin,VRmax,-1000.0,1000.0);
    } else {
      Regulator_gain_blk.setparams(KA,VRmin,VRmax);
    }

    double a[2],b[2];
    a[0] = TF; a[1] = 1.0;
    b[0] = KF;  b[1] = 0.0;
    Feedback_blk.setcoeffs(a,b);

    Output_blk.setparams(TE);
  }
}

/**
 * Initialize exciter model before calculation
 * @param [output] values - array where initialized exciter variables should be set
 */
void Ieeet1::init(gridpack::RealType* xin) 
{
  gridpack::RealType *x = xin+offsetb; // exciter array starts from this location
  double Ec = sqrt(VD*VD + VQ*VQ);
  double VRin;
  double Vf=0.0;
  double SE=0.0;

  // Get initial field voltage
  Efd = getInitialFieldVoltage();

  if(integrationtype != IMPLICIT) {
    // Initialization for explicit integration
    VF = Feedback_blk.init_given_u(Efd);

    double output_blk_in;
    output_blk_in = Output_blk.init_given_y(Efd);
    
    double sat_signal = KE*Efd;

    if(has_Sat) {
      if(Efd > Efdthresh) SE = satB*(Efd - satA)*(Efd - satA)/Efd;
    }

    VR = output_blk_in + (KE + SE)*Efd;

    double Regulator_blk_in;
    if(zero_TA) {
      Regulator_blk_in = VR/KA;
    } else {
      Regulator_blk_in = Regulator_blk.init_given_y(VR);
    }
    
    double Verr;
    Verr = Regulator_blk_in + VF;
    
    if(!zero_TR) {
      Vmeas = Vmeas_blk.init_given_u(Ec);
    } else Vmeas = Ec;
    
    Vref = Verr + Vmeas - Vs;

  } else {
    Vmeas    = Ec;
    xf       = -KF/TF*Efd;
    
    if(Efd > Efdthresh) SE = satB*(Efd - satA)*(Efd - satA)/Efd;
    VR = (SE + KE)*Efd;
    
    VRin     = VR/KA;

    Vref = VRin + Vmeas + Vf;

    x[0] = Vmeas;
    x[1] = VR;
    x[2] = Efd;
    x[3] = xf;
  }
}

/**
 * Write output from exciters to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Ieeet1::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out exciter state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Ieeet1::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto exciters
 * @param values array containing exciter state variables
*/
void Ieeet1::setValues(gridpack::RealType *val)
{
  gridpack::RealType *values = val+offsetb; // exciter array starts from this location

  if(integrationtype == EXPLICIT) return;

  if(p_mode == XVECTOBUS) {
    Vmeas = values[0];
    VR    = values[1];
    Efd   = values[2];
    xf    = values[3];
  } else if(p_mode == XDOTVECTOBUS) {
    dVmeas = values[0];
    dVR    = values[1];
    dEfd   = values[2];
    dxf    = values[3];
  }
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Ieeet1::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // exciter array starts from this location

  if(integrationtype == EXPLICIT) return;

  double Ec,yLL,Vf,SE=0.0, VRin;

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
    
    VRin = Vref - Vmeas - Vf;

    // VR equation
    if(VR_at_min) {
      f[1] = VR - VRmin;
    } else if(VR_at_max) {
      f[1] = VR - VRmax;
    } else {
      if(TA != 0) f[1] = (-VR + KA*VRin)/TA - dVR;
      else f[1] = -VR + KA*VRin;
    }

    if(Efd > Efdthresh) SE = satB*(Efd - satA)*(Efd - satA)/Efd;
    f[2] = (VR - (SE + KE)*Efd)/TE - dEfd;

    // xf equation
    f[3] = (-xf - KF/TF*Efd)/TF - dxf;    
  }
}

/**
   Non-zero pattern of the Jacobian (x denotes non-zero entry)
         Vmeas    VR    Efd    xf    delta    va    vb    vc
 eq.0 |    x                            x      x     x     x
 eq.1 |    x      x      x     x
 eq.2 |           x      x
 eq.3 |                  x     x

 Number of non-zeros = 5 + 4 + 2 + 2 = 13 
 * Get number of matrix values contributed by exciter
 * @return number of matrix values
 */
int Ieeet1::matrixNumValues()
{
  int nmat = 0;
  if(integrationtype == IMPLICIT) nmat = 13 + 5;
  return nmat;
}

  /**
   * Return values from Jacobian matrix
   * @param nvals: number of values to be inserted
   * @param values: pointer to matrix block values
   * @param rows: pointer to matrix block rows
   * @param cols: pointer to matrix block cols
   */
void Ieeet1::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;

  if(integrationtype != IMPLICIT) {
    *nvals = ctr;
    return;
  }
  
  int Vmeas_idx = p_gloc;
  int VR_idx    = p_gloc + 1;
  int Efd_idx   = p_gloc + 2;
  int xf_idx    = p_gloc + 3;
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

  rows[ctr]   = VR_idx; cols[ctr] = Vmeas_idx;
  rows[ctr+1] = VR_idx; cols[ctr+1] = VR_idx;
  rows[ctr+2] = VR_idx; cols[ctr+2] = Efd_idx;
  rows[ctr+3] = VR_idx; cols[ctr+3] = xf_idx;

  values[ctr] = values[ctr+1] = values[ctr+2] = values[ctr+3] = 0.0;

  if(VR_at_min || VR_at_max) {
    values[ctr+1] = 1.0;
  } else {
    if(TA != 0) {
      values[ctr] =   -KA/TA;
      values[ctr+1] = -1.0/TA - shift;
      values[ctr+2] = -(KA*KF/TF)/TA;
      values[ctr+3] = -KA/TA;
    } else {
      values[ctr] =   -KA/TA;
      values[ctr+1] = -1.0/TA;
      values[ctr+2] = -(KA*KF/TF)/TA;
      values[ctr+3] = -KA/TA;
    }
  }
      
  ctr += 4;

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
  
  *nvals = ctr;
		     
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double Ieeet1::getFieldVoltage()
{
  return Efd;
}

/** 
 * Get the value of the field voltage parameter
 * and its global location
 * @return value of field voltage
 */
double Ieeet1::getFieldVoltage(int *Efd_gloc)
{
  if(integrationtype == IMPLICIT) {
    *Efd_gloc = p_gloc + 2;
  } else {
    *Efd_gloc = -1;
  }
  return Efd;
}


bool Ieeet1::getFieldVoltagePartialDerivatives(int *xexc_loc,double *dEfd_dxexc,double *dEfd_dxgen)
{
  return false;
}

/**
 * Update the event function values
 */
void Ieeet1::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues)
{
  int offset    = getLocalOffset();
  int Vmeas_idx = offset;
  int VR_idx    = offset+1;
  int Efd_idx   = offset+2;
  int xf_idx    = offset+3;

  Vmeas = state[Vmeas_idx];
  VR    = state[VR_idx];
  Efd   = state[Efd_idx];
  xf    = state[xf_idx];

  /* Only considering limits on VR */
  double Vf,yLL,dVR_dt, VRin;
  Vf = xf + KF/TF*Efd;
  VRin = Vref - Vmeas - Vf;

  dVR_dt = (-VR + KA*VRin)/TA;

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
 * Event handler
 */
void Ieeet1::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  int offset    = getLocalOffset();
  int Vmeas_idx = offset;
  int VR_idx    = offset+1;
  int Efd_idx   = offset+2;
  int xf_idx    = offset+3;

  Vmeas = state[Vmeas_idx];
  VR    = state[VR_idx];
  Efd   = state[Efd_idx];
  xf    = state[xf_idx];

  /* Only considering limits on VR */
  double Vf,yLL,dVR_dt,VRin;
  Vf = xf + KF/TF*Efd;
  VRin = Vref - Vmeas - Vf;
  dVR_dt = (-VR + KA*VRin)/TA;

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
void Ieeet1::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  if(integrationtype == IMPLICIT) {
    gridpack::math::RealDAESolver::EventPtr e(new Ieeet1Event(this));

    eman->add(e);
  }
}

void Ieeet1Event::p_update(const double& t,gridpack::RealType *state)
{
  p_exc->eventFunction(t,state,p_current);
}

void Ieeet1Event::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  p_exc->eventHandlerFunction(triggered,t,state);
}
