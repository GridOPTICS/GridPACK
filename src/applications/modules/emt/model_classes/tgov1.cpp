/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   tgov1.cpp
 *  
 * @brief TGOV1 governor model implementation 
 *
 */

#include <tgov1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Tgov1::Tgov1(void)
{
  x1 = 0.0; 
  x2 = 0.0;
  xout = 0.0;

  dx1 = 0.0;
  dx2 = 0.0;

  R = 0.0;
  T1 = 0.0;
  Vmax = 1000.0;
  Vmin = -1000.0;
  T2 = 0.0;
  T3 = 0.0;
  Dt = 0.0;
  x1_at_min = x1_at_max = false;

  nxgov = 3; // Number of variables
}

Tgov1::~Tgov1(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to BaseGovernorModel
 */
void Tgov1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGovModel::load(data,idx); // load parameters in base governor model

  if (!data->getValue(GOVERNOR_R, &R, idx)) R = 0.05; 
  if (!data->getValue(GOVERNOR_T1, &T1, idx)) T1 = 0.5; 
  if (!data->getValue(GOVERNOR_VMAX, &Vmax, idx)) Vmax = 1.0; 
  if (!data->getValue(GOVERNOR_VMIN, &Vmin, idx)) Vmin = 0.0; 
  if (!data->getValue(GOVERNOR_T2, &T2, idx)) T2 = 3.0; 
  if (!data->getValue(GOVERNOR_T3, &T3, idx)) T3 = 10.0; 
  if (!data->getValue(GOVERNOR_DT, &Dt, idx)) Dt = 0.0;

  // Set up transfer function blocks
  if(integrationtype != IMPLICIT) {
    leadlag_blk.setparams(T2,T3);
    delay_blk.setparams(1.0,T1,Vmin,Vmax,-1000.0,1000.0);
  }
}

/**
   Number of variables
*/ 
void Tgov1::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxgov = 0;

  *nvar = nxgov;
}
  
/**
   Prestep function
*/
void Tgov1::preStep(double time ,double timestep)
{
  if(integrationtype != EXPLICIT) return;

  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();
  
  // Delay block output and state update
  delay_blk_out = delay_blk.getoutput((Pref-dw)/R,timestep,true);

  // Leadlag block output and state update
  leadlag_blk_out = leadlag_blk.getoutput(delay_blk_out,timestep,true);

  // Output mechanical power
  Pmech = leadlag_blk_out - Dt*dw;

}

/**
   Poststep function
*/
void Tgov1::postStep(double time)
{
}


/**
 * Initialize governor model before calculation
 * @param [output] values - array where initialized governor variables should be set
 */
void Tgov1::init(gridpack::RealType* xin) 
{
  gridpack::RealType *x = xin+offsetb; // governor array starts from this location
  
  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();
  double Pmech = gen->getInitialMechanicalPower();

  if(integrationtype != IMPLICIT) {
    // Initialize leadlag block
    delay_blk_out = leadlag_blk.init_given_y(Pmech-Dt*dw);
    
    double delay_blk_in;
    
    // Initialize delay block
    delay_blk_in = delay_blk.init_given_y(delay_blk_out);
    
    // Reference power signal
    Pref = R*delay_blk_in + dw;

    return;
  }

  x2 = Pmech;
  x1 = Pmech;
  xout = Pmech;

  Pref = R*Pmech + dw;

  if(integrationtype == IMPLICIT) {
    x[0] = x1;
    x[1] = x2;
    x[2] = xout;
  }
}

/**
 * Write output from governors to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Tgov1::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out governor state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Tgov1::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto governors
 * @param values array containing governor state variables
*/
void Tgov1::setValues(gridpack::RealType *val)
{
  gridpack::RealType *values = val+offsetb; // governor array starts from this location

  if(integrationtype == EXPLICIT) return;

  if(p_mode == XVECTOBUS) {
    x1 = values[0];
    x2 = values[1];
    xout = values[2];
  } else if(p_mode == XDOTVECTOBUS) {
    dx1 = values[0];
    dx2 = values[1];
  }
}

/**
 * Return the values of the governor vector block
 * @param values: pointer to vector values
 */
void Tgov1::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // exciter array starts from this location

  if(integrationtype == EXPLICIT) return;
  
  int x1_idx = 0;
  int x2_idx = 1;
  int xout_idx = 2;

  double yLL;
  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();

  if(p_mode == RESIDUAL_EVAL) {
    // x1 equation
    f[x1_idx] = 1/T1*((1/R)*(Pref - dw) - x1) - dx1;
    f[x2_idx] = (-x2 + (1.0 - T2/T3)*x1)/T3 - dx2;

    Pmech = x2 + T2/T3*x1 - Dt*dw;

    f[xout_idx] = Pmech - xout;
  }
}

/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values
 */
int Tgov1::matrixNumValues()
{
  if(integrationtype == IMPLICIT) return 8;
  else return 0;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Tgov1::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;

  if(integrationtype != IMPLICIT) {
    *nvals = 0;
    return;
  }
  
  int x1_gloc = p_gloc;
  int x2_gloc = p_gloc+1;
  int xout_gloc = p_gloc+2;
  int    dw_gloc;
  double dw;

  dw = getGenerator()->getSpeedDeviation(&dw_gloc);

  // Partial derivatives of 1 equation
  rows[ctr] = x1_gloc; cols[ctr] = x1_gloc;
  rows[ctr+1] = x1_gloc; cols[ctr+1] = dw_gloc;
  
  values[ctr] = values[ctr+1] = 0.0;
  
  values[ctr] = -1.0/T1 - shift;
  values[ctr+1] = -1/T1*(1/R);

  ctr += 2;

  rows[ctr]   = x2_gloc; cols[ctr] = x1_gloc;
  rows[ctr+1] = x2_gloc; cols[ctr+1] = x2_gloc;

  values[ctr] = values[ctr+1] = 0.0;

  values[ctr] = (1.0 - T2/T3)/T3;
  values[ctr+1] = -1.0/T3 - shift;

  ctr += 2;

  rows[ctr]   = xout_gloc; cols[ctr] = x1_gloc;
  rows[ctr+1] = xout_gloc; cols[ctr+1] = x2_gloc;
  rows[ctr+2] = xout_gloc; cols[ctr+2] = xout_gloc;
  rows[ctr+3] = xout_gloc; cols[ctr+3] = dw_gloc;

  values[ctr]   = T2/T3;
  values[ctr+1] = 1.0;
  values[ctr+2] = -1.0;
  values[ctr+3] = -Dt;
  
  ctr += 4;

  *nvals = ctr;
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power 
 */
void Tgov1::setInitialMechanicalPower(double Pmech0)
{
  Pmech = Pmech0;
}

/** 
 * Get the value of the mechanical power parameter
 * @return value of the mechanical power 
 */
double Tgov1::getMechanicalPower()
{
  return xout;
}

/** 
 * Get the value of the mechanical power and its global location
 * @return value of the mechanical power
 *
 * Note: Used in Jacobian calculation
 */
double Tgov1::getMechanicalPower(int *Pmech_gloc)
{
  if(integrationtype == IMPLICIT) *Pmech_gloc = p_gloc + 2;
  else *Pmech_gloc = -1;

  return getMechanicalPower();
}


/**
 * Partial derivatives of Mechanical Power Pmech w.r.t. governor variables
 * @param xgov_loc locations of governor variables
 * @param dPmech_dxgov partial derivatives of mechanical power Pmech w.r.t governor variables
*/
bool Tgov1::getMechanicalPowerPartialDerivatives(int *xgov_loc,double *dPmech_dxgov)
{
  return false;
}

void Tgov1::setVcomp(double Vcomp)
{
}

/**
 * Update the event function values
 */
void Tgov1::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues)
{
  int offset   = getLocalOffset();
  int x1_idx   = offset;

  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();

  x1  = state[x1_idx];

  double dx1_dt = ((1/R)*(Pref - dw) - x1);
;

  /* Limits on x1 */
  if(!x1_at_min) {
    evalues[0] = x1 - Vmin;
  } else {
    evalues[0] = -dx1_dt; /* Release when derivative reaches 0 */
  }

  if(!x1_at_max) {
    evalues[1] = Vmax - x1;
  } else {
    evalues[1] = dx1_dt; /* Release when derivative reaches 0 */
  }
} 

/**
 * Event handler
 */
void Tgov1::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  int offset    = getLocalOffset();
  int x1_idx = offset;

  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();

  x1  = state[x1_idx];

  double dx1_dt = ((1/R)*(Pref - dw) - x1);

  if(triggered[0]) {
    if(!x1_at_min && dx1_dt < 0) {
      /* Hold x1 at Vmin */
      x1_at_min = true;
    } else {
      /* Release */
      x1_at_min = false;
    }
  }

  if(triggered[1]) {
    if(!x1_at_max && dx1_dt > 0) {
      /* Hold x1 at Pmax */
      x1_at_max = true;
    } else {
      /* Release */
      x1_at_max = false;
    }
  }
}

/**
 * Set event
 */
void Tgov1::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  if(integrationtype == IMPLICIT) {
    gridpack::math::RealDAESolver::EventPtr e(new Tgov1Event(this));

    eman->add(e);
  }
}

void Tgov1Event::p_update(const double& t,gridpack::RealType *state)
{
  p_gov->eventFunction(t,state,p_current);
}

void Tgov1Event::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  p_gov->eventHandlerFunction(triggered,t,state);
}
