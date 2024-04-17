/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   sexs.cpp
 * 
 * @brief  
 * 
 * @ Updated by Shuangshuang Jin on April 17, 2024.
 */

#include <sexs.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Sexs::Sexs(void)
{
  Vmeas = 0.0;
  xLL = 0.0;
  dxLL = 0.0;
  Efd = 0.0;
  dEfd = 0.0;
  TA_OVER_TB = 0.0;
  TB = 0.0;
  K = 0.0;
  TE = 0.0;
  EMIN = 0.0;
  EMAX = 0.0;

  zero_TE = false; 

  Efd_at_min = Efd_at_max = false;

  nxexc = 3;
}

Sexs::~Sexs(void)
{
}

void Sexs::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxexc = 0;
  *nvar = nxexc;
}

void Sexs::preStep(double time, double timestep)
{
  if(integrationtype != EXPLICIT) return;

  double vabc[3],vdq0[3];

  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;

  double delta = getGenerator()->getAngle();
  
  abc2dq0(vabc,p_time,delta,vdq0);
  double Vd, Vq;
  Vd = vdq0[0]; Vq = vdq0[1];
  
  Ec = sqrt(Vd*Vd + Vq*Vq);

  Vmeas = Ec;
  double Verr = Vref - Vmeas + Vs;

  double filter_blk_in;
  filter_blk_in = leadlagblock.getoutput(Verr,timestep,true);

  if(!zero_TE) {
    Efd = filterblock.getoutput(filter_blk_in,timestep,true);
  } else {
    Efd = gainblock.getoutput(filter_blk_in,EMIN,EMAX);
  }

}

void Sexs::postStep(double time)
{

}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to
 * SexsModel
 */
void Sexs::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  if (!data->getValue(EXCITER_TA_OVER_TB, &TA_OVER_TB, idx)) TA_OVER_TB = 0.0; // TA_OVER_TB
  if (!data->getValue(EXCITER_TB, &TB, idx)) TB = 0.0; // TB
  if (!data->getValue(EXCITER_K, &K, idx))   K  = 0.0; // K
  if (!data->getValue(EXCITER_EMAX, &EMAX, idx)) EMAX = 0.0; // EMAX
  if (!data->getValue(EXCITER_EMIN, &EMIN, idx)) EMIN = 0.0; // EMIN
  if (!data->getValue(EXCITER_TE, &TE, idx)) TE = 0.0; // TE
  //printf("%f, %f, %f, %f, %f, %f\n", TA_OVER_TB, TB, K, EMAX, EMIN, TE);

  TA = TA_OVER_TB*TB;

  if(integrationtype != IMPLICIT) {
    // Set up blocks
    // Set parameters for the first block
    leadlagblock.setparams(TA,TB);

    zero_TE = false;
    if(fabs(TE) < 1e-6) {
      zero_TE = true;
    }
    
    // Set parameters for the second block
    if(!zero_TE) {
      filterblock.setparams(K,TE,EMIN,EMAX,-1000.0,1000);
    } else {
      gainblock.setparams(K,EMIN,EMAX);
    }
  }

  Vs = 0.0; 
}


/**
 * Initialize exciter model before calculation
 * @para
 */
void Sexs::init(gridpack::RealType* xin)
{
  gridpack::RealType *x = xin+offsetb; // exciter array starts from this location
  double Ec = sqrt(VD*VD + VQ*VQ);

  // Get initial field voltage
  Efd = getInitialFieldVoltage();

  if(integrationtype != IMPLICIT) {
    // Initialization for explicit integration
  
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

    Vref = Ec + u1 - Vs;   // Voltage reference initial value
  } else {
      double yLL;

      Vmeas = Ec;

      // Efd is already set by the generator model
      yLL = Efd / K;

      Vref = yLL + Vmeas - Vs;

      if(TB != 0 && TA != 0) xLL = (1 - TA/TB)*(Vref - Vmeas + Vs);
      else xLL = Vref - Vmeas + Vs;

      x[0] = Vmeas;
      x[1] = xLL;
      x[2] = Efd;
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
bool Sexs::serialWrite(char *string, const int bufsize,const char *signal)
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
void Sexs::setValues(gridpack::RealType *val)
{
  gridpack::RealType *values = val+offsetb; // exciter array starts from this location

  if(integrationtype == EXPLICIT) return;

  if(p_mode == XVECTOBUS) {
    Vmeas = values[0];
    xLL = values[1];
    Efd = values[2];
  } else if(p_mode == XDOTVECTOBUS) {
    dxLL = values[1];
    dEfd = values[2];
  }
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Sexs::vectorGetValues(gridpack::RealType *values)
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
    f[0] = -Vmeas + Ec;

    // xLL equation
    if(TB != 0 && TA != 0) {
      f[1] = (-xLL + (1 - TA/TB)*(Vref - Vmeas + Vs))/TB - dxLL;
      yLL = xLL + TA/TB*(Vref - Vmeas + Vs);
    } else {
      f[1] = -xLL + Vref - Vmeas + Vs;
      yLL = xLL;
    }

    // Efd equation
    if(Efd_at_min) {
      f[2] = Efd - EMIN;
    } else if(Efd_at_max) {
      f[2] = Efd - EMAX;
    } else {
      if(TE != 0) f[2] = (-Efd + K*yLL)/TE - dEfd;
      else f[2] = -Efd + K*yLL;
    }

  }
}


/**
   Non-zero pattern of the Jacobian (x denotes non-zero entry)
         Vmeas   xLL     Efd     delta    va    vb    vc
 eq.0 |    x                       x       x     x     x  
 eq.1 |    x      x              
 eq.2 |    x      x       x      

 Number of non-zeros = 5 + 2 + 3 = 10 
 * Get number of matrix values contributed by exciter
 * @return number of matrix values
 */
int Sexs::matrixNumValues()
{
  int nmat = 0;
  if(integrationtype == IMPLICIT) nmat = 10;
  return nmat;
}


/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Sexs::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;

  if(integrationtype != IMPLICIT) {
    *nvals = ctr;
    return;
  }
    
  int xLL_idx = p_gloc;
  int Efd_idx = p_gloc + 1;
  int Vmeas_idx = p_gloc + 2;
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
  rows[ctr+1] = xLL_idx; cols[ctr+1] = delta_idx;
  rows[ctr+2] = xLL_idx; cols[ctr+2] = va_idx;
  rows[ctr+3] = xLL_idx; cols[ctr+3] = vb_idx;
  rows[ctr+4] = xLL_idx; cols[ctr+4] = vb_idx;

  values[ctr]   = -1.0;
  values[ctr+1] = (dEc_dVd*dVd_ddelta + dEc_dVq*dVq_ddelta);
  values[ctr+2] = (dEc_dVd*dVd_dvabc[0] + dEc_dVq*dVq_dvabc[0]);
  values[ctr+3] = (dEc_dVd*dVd_dvabc[1] + dEc_dVq*dVq_dvabc[1]);
  values[ctr+4] = (dEc_dVd*dVd_dvabc[2] + dEc_dVq*dVq_dvabc[2]);
  
  ctr += 5;

  rows[ctr]   = xLL_idx; cols[ctr] = Vmeas_idx;
  rows[ctr+1]   = xLL_idx; cols[ctr+1] = xLL_idx;

  values[ctr] = values[ctr+1] = 0.0;

  double param = (1 - TA/TB);
  double dyLL_dVmeas = 0.0, dyLL_dxLL = 0.0;
  if(TA != 0 && TB != 0) {
    values[ctr] = (param*-1)/TB;
    values[ctr+1] = -1.0/TB - shift;

    dyLL_dVmeas = TA/TB*-1.0;
    dyLL_dxLL = 1.0;
  } else {
    values[ctr] = -1.0;
    values[ctr+1] = -1.0;

    dyLL_dxLL = 1.0;
  }

  ctr += 2;

  rows[ctr] = Efd_idx;  cols[ctr] = Vmeas_idx;
  rows[ctr] = Efd_idx;  cols[ctr] = xLL_idx;
  rows[ctr+1] = Efd_idx; cols[ctr+1] = Efd_idx;

  values[ctr] = values[ctr+1] = values[ctr+2] = 0.0;

  if(Efd_at_min || Efd_at_max) {
    values[ctr+2] = 1.0;
  } else {
    if(TE != 0) {
      values[ctr] = (K*dyLL_dVmeas)/TE;
      values[ctr+1] = (K*dyLL_dxLL)/TE;
      values[ctr+2] = -1.0/TE - shift;
    } else {
      values[ctr] =   (K*dyLL_dVmeas);
      values[ctr+1] = (K*dyLL_dxLL);
      values[ctr+2] = -1.0/TE;
    }
  }
      
  ctr += 3;

  *nvals = ctr;

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
  if(integrationtype == IMPLICIT) {
    *Efd_gloc = p_gloc + 1;
  } else {
    *Efd_gloc = -1;
  }
  return Efd;
}

bool Sexs::getFieldVoltagePartialDerivatives(int *xexc_loc,double *dEfd_dxexc,double *dEfd_dxgen)
{
  return false;
}

/**
 * Update the event function values
 */
void Sexs::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType>& evalues)
{

}

/**
 * Event handler
 */
void Sexs::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{

}

/**
 * Set event
 */
void Sexs::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{

}

/**
 * Set the field voltage parameter inside the exciter
 * @param fldv value of the field voltage
 */
void Sexs::setFieldVoltage(double fldv)
{
  // This is the initial value of Efd using during initialization
  Efd = fldv;
}






