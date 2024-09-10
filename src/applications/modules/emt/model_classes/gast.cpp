/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   gast.cpp
 *  
 * @brief GAST governor model implementation 
 *
 */

#include <gast.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>
#include <algorithm>

Gast::Gast(void)
{
  x1 = 0.0; 
  x2 = 0.0;
  x3 = 0.0;
  xout = 0.0;

  dx1 = 0.0;
  dx2 = 0.0;
  dx3 = 0.0;

  R = 0.0;
  T1 = 0.0;
  T2 = 0.0;
  T3 = 0.0;
  AT = 0.0;
  KT = 0.0;
  Vmax = 1000.0;
  Vmin = -1000.0;
  Dt = 0.0;
  
  Loadref = 0.0;
  x1_at_min = x1_at_max = false;

  nxgov = 4; // Number of variables
}

Gast::~Gast(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to BaseGovernorModel
 */
void Gast::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGovModel::load(data,idx); // load parameters in base governor model

  // Yuan: will need to add the new parameters into data parser
  if (!data->getValue(GOVERNOR_R, &R, idx)) R = 0.05; 
  if (!data->getValue(GOVERNOR_T1, &T1, idx)) T1 = 0.5; 
  if (!data->getValue(GOVERNOR_T2, &T2, idx)) T2 = 3.0; 
  if (!data->getValue(GOVERNOR_T3, &T3, idx)) T3 = 10.0; 
  if (!data->getValue(GOVERNOR_AT, &AT, idx)) AT = 10.0; 
  if (!data->getValue(GOVERNOR_KT, &KT, idx)) KT = 10.0; 
  if (!data->getValue(GOVERNOR_VMAX, &Vmax, idx)) Vmax = 1.0; 
  if (!data->getValue(GOVERNOR_VMIN, &Vmin, idx)) Vmin = 0.0; 
  if (!data->getValue(GOVERNOR_DT, &Dt, idx)) Dt = 0.0;
  
  // Set up transfer function blocks
  if(integrationtype != IMPLICIT) {
    /* Create string for setting name */
    std::string blkhead = std::to_string(busnum) + "_" + id + "GAST_";

    std::string delay_block_T1_name = blkhead + "delay_blk_T1";
    delay_blk_T1.setname(delay_block_T1_name.c_str());
    delay_blk_T1.setparams(1.0,T1,Vmin,Vmax,-1000.0,1000.0);
    
    delay_blk_T2.setparams(1.0,T2);
    delay_blk_T3.setparams(1.0,T3);
  }
}

/**
   Number of variables
*/ 
void Gast::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxgov = 0;

  *nvar = nxgov;
}
  
/**
   Prestep function
*/
void Gast::preStep(double time ,double timestep)
{
  if(integrationtype != EXPLICIT) return;

  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();
  
  // Low value gate
  double LVG_in_1, LVG_in_2, LVG_out;
  delay_blk_T3_out = delay_blk_T3.getoutput(delay_blk_T2_out);
  LVG_in_2 = KT * (AT - delay_blk_T3_out) + AT;
  LVG_in_1 = Loadref - dw/R;
  LVG_out = std::min(LVG_in_1, LVG_in_2);
  
  // Delay block T1 output and state update
  delay_blk_T1_out = delay_blk_T1.getoutput(LVG_out,timestep,true);

  // Leadlag block output and state update
  delay_blk_T2_out = delay_blk_T2.getoutput(delay_blk_T1_out,timestep,true);

  // Output mechanical power
  Pmech = delay_blk_T2_out - Dt*dw;

  xout = Pmech;

}

/**
   Poststep function
*/
void Gast::postStep(double time)
{
}


/**
 * Initialize governor model before calculation
 * @param [output] values - array where initialized governor variables should be set
 */
void Gast::init(gridpack::RealType* xin) 
{
  gridpack::RealType *x = xin+offsetb; // governor array starts from this location
  
  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();
  double Pmech = gen->getInitialMechanicalPower();

  if(integrationtype != IMPLICIT) {

    // Initialize delay block T2
    delay_blk_T2_out = Pmech+Dt*dw;
    delay_blk_T1_out = delay_blk_T2.init_given_y(delay_blk_T2_out);
    
    // Initialize delay block T1
    double delay_blk_T1_in;
    delay_blk_T1_in = delay_blk_T1.init_given_y(delay_blk_T1_out);
    
    // Initialize delay block T3
    delay_blk_T3_out = delay_blk_T3.init_given_u(delay_blk_T2_out);
    
    // Load reference signal
    Loadref = delay_blk_T1_in + dw/R;

    return;
  }

  x1 = Pmech;
  x2 = Pmech;
  x3 = Pmech;
  xout = Pmech;

  Loadref = Pmech + dw/R;
  
  x[0] = x1;
  x[1] = x2;
  x[2] = x3;
  x[3] = xout;
}

/**
 * Write output from governors to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Gast::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out governor state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Gast::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto governors
 * @param values array containing governor state variables
*/
void Gast::setValues(gridpack::RealType *val)
{
  gridpack::RealType *values = val+offsetb; // governor array starts from this location

  if(integrationtype == EXPLICIT) return;
  if(p_mode == XVECTOBUS) {
    x1 = values[0];
    x2 = values[1];
    x3 = values[2];
    xout = values[3];
  } else if(p_mode == XDOTVECTOBUS) {
    dx1 = values[0];
    dx2 = values[1];
    dx3 = values[2];
  }
}

/**
 * Return the values of the governor vector block
 * @param values: pointer to vector values
 */
void Gast::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // exciter array starts from this location

  if(integrationtype == EXPLICIT) return;
  
  int x1_idx = 0;
  int x2_idx = 1;
  int x3_idx = 2;
  int xout_idx = 3;

  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();

  if(p_mode == RESIDUAL_EVAL) {
    // x1 equation
    if(Loadref - 1/R*dw <= KT*(AT-x3) + AT) {
        if (x1_at_min) {
            f[x1_idx] = x1 - Vmin;
        } else if (x1_at_max) {
            f[x1_idx] = x1 - Vmax;
        } else {
            f[x1_idx] = 1 / T1 * (Loadref - 1 / R * dw) - 1 / T1 * x1 - dx1;
        }
    } else {
        if (x1_at_min) {
            f[x1_idx] = x1 - Vmin;
        } else if (x1_at_max) {
            f[x1_idx] = x1 - Vmax;
        } else {
            f[x1_idx] = 1/T1*(KT * (AT-x3) + AT) - 1/T1*x1 - dx1;
        }
    }
    // x2 equation
    f[x2_idx] = 1/T2*x1 - 1/T2*x2 - dx2;
    // x3 equation
    f[x3_idx] = 1/T3*x2 - 1/T3*x3 - dx3;
    // xout equation
    f[xout_idx] = x2 - Dt * dw - xout;
  }
}

/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values
 */
int Gast::matrixNumValues()
{ // Yuan: TODO: need to figure it out. # of non-zero values in matrixGetValues?
  // Yuan: return the # of assigned value positions? 8 for TGOV1
  if(integrationtype == IMPLICIT) return 10; 
  else return 0;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Gast::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;

  if(integrationtype != IMPLICIT) {
    *nvals = 0;
    return;
  }
  
  int x1_gloc = p_gloc;
  int x2_gloc = p_gloc+1;
  int x3_gloc = p_gloc+2;
  int xout_gloc = p_gloc+3;
  int    dw_gloc;
  double dw;

  dw = getGenerator()->getSpeedDeviation(&dw_gloc);

  // Row 1 of Jacobian matrix -> residual function f_x1
  rows[ctr] = x1_gloc; cols[ctr] = x1_gloc; 
  rows[ctr+1] = x1_gloc; cols[ctr+1] = x3_gloc;
  rows[ctr+2] = x1_gloc; cols[ctr+2] = dw_gloc;
  
  values[ctr] = values[ctr+1] = values[ctr+2] = 0.0;
  
  values[ctr] = -1.0/T1 - shift; 
  if(Loadref - 1/R*dw <= KT*(AT-x3) + AT) {
      values[ctr+1] = 0.0;
      values[ctr+2] = -1/T1*(1/R);
  }
  else {
      values[ctr+1] = -KT/T1;
      values[ctr+2] = 0.0;
  }
  
  ctr += 3;

  // Row 2 of Jacobian matrix -> f_x2
  rows[ctr]   = x2_gloc; cols[ctr] = x1_gloc;
  rows[ctr+1] = x2_gloc; cols[ctr+1] = x2_gloc;

  values[ctr] = values[ctr+1] = 0.0;

  values[ctr] = 1.0/T2;
  values[ctr+1] = -1.0/T2 - shift;

  ctr += 2;
  
  // Row 3 of Jacobian matrix -> f_x3
  rows[ctr] = x3_gloc; cols[ctr] = x2_gloc;
  rows[ctr+1] = x3_gloc; cols[ctr+1] = x3_gloc;
  
  values[ctr] = values[ctr+1] = 0.0;
  
  values[ctr] = 1.0/T3;
  values[ctr+1] = -1.0/T3 - shift;
  
  ctr += 2;

  // Row 4 of Jacobian matrix  -> f_xout
  rows[ctr] = xout_gloc; cols[ctr] = x2_gloc;
  rows[ctr+1] = xout_gloc; cols[ctr+1] = xout_gloc;
  rows[ctr+2] = xout_gloc; cols[ctr+2] = dw_gloc;
  
  values[ctr] = values[ctr+1] = values[ctr+2] = 0.0;

  values[ctr] = 1.0;
  values[ctr+1] = -1.0;
  values[ctr+2] = Dt;
  
  ctr += 3;

  *nvals = ctr;
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power 
 */
void Gast::setInitialMechanicalPower(double Pmech0)
{
  Pmech = Pmech0;
}

/** 
 * Get the value of the mechanical power parameter
 * @return value of the mechanical power 
 */
double Gast::getMechanicalPower()
{
  return xout;
}

/** 
 * Get the value of the mechanical power and its global location
 * @return value of the mechanical power
 *
 * Note: Used in Jacobian calculation
 */
double Gast::getMechanicalPower(int *Pmech_gloc)
{
  if(integrationtype == IMPLICIT) *Pmech_gloc = p_gloc + 3; // Yuan: why plus 2 here?
  else *Pmech_gloc = -1;

  return getMechanicalPower();
}


/**
 * Partial derivatives of Mechanical Power Pmech w.r.t. governor variables
 * @param xgov_loc locations of governor variables
 * @param dPmech_dxgov partial derivatives of mechanical power Pmech w.r.t governor variables
*/
bool Gast::getMechanicalPowerPartialDerivatives(int *xgov_loc,double *dPmech_dxgov)
{
  return false;
}

void Gast::setVcomp(double Vcomp)
{
}

/**
 * Update the event function values
 */
void Gast::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues)
{
  int offset   = getLocalOffset();
  int x1_idx   = offset;

  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();

  x1  = state[x1_idx];

  double dx1_dt;
  
  if(Loadref - 1/R*dw <= KT*(AT-x3) + AT)
      dx1_dt = (Loadref - dw/R - x1)/T1;
  else
      dx1_dt = (KT*(AT-x3)+AT - x1)/T1;

  /* Limits on x1 */
  // Yuan: don't quite understand this part. What does evalues do? 
  if(!x1_at_min) {
    // event function: detect if x1 has crossed Vmin, if x1>Vmin, no action, only if evalues[idx]=x1-Vmin changes from positive to negative, event handler triggered.
    evalues[0] = x1 - Vmin; 
  } else {
    // when x1 is at min, we detect if it is time to release it based on the derivative if x1 is at min and dx1_dt<0, that means it won't deviate from min
    // but when x1 is at min and dx1_dt change from <0 to >0, then evalues[idx] change from >0 to <0. Since we are using "CrossZeroNegative" to triggered
    // events as defined in the gast.hpp file. This "exiting min" event will be detected and flagged.
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
void Gast::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  int offset    = getLocalOffset();
  int x1_idx = offset;

  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();

  x1  = state[x1_idx];

  double dx1_dt;
  
  if(Loadref - 1/R*dw <= KT*(AT-x3) + AT)
      dx1_dt = (Loadref - dw/R - x1)/T1;
  else
      dx1_dt = (KT*(AT-x3)+AT - x1)/T1;

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
void Gast::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  if(integrationtype == IMPLICIT) {
    gridpack::math::RealDAESolver::EventPtr e(new GastEvent(this));

    eman->add(e);
  }
}

void GastEvent::p_update(const double& t,gridpack::RealType *state)
{
  p_gov->eventFunction(t,state,p_current);
}

void GastEvent::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  p_gov->eventHandlerFunction(triggered,t,state);
}
