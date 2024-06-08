/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hygov.cpp
 *  
 * @brief TGOV1 governor model implementation 
 *
 */

#include <hygov.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Hygov::Hygov(void)
{
  //x1 = 0.0; 
  //x2 = 0.0;
  //xout = 0.0;

  //dx1 = 0.0;
  //dx2 = 0.0;

  //R = 0.0;
  //T1 = 0.0;
  //Vmax = 1000.0;
  //Vmin = -1000.0;
  //T2 = 0.0;
  //T3 = 0.0;
  //Dt = 0.0;
  //x1_at_min = x1_at_max = false;

  //nxgov = 3; // Number of variables

  Pmech = 1.0;
  nref = 1.0;
  delta_w = 0.0;
}

Hygov::~Hygov(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to BaseGovernorModel
 */
void Hygov::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  //BaseEMTGovModel::load(data,idx); // load parameters in base governor model

  //if (!data->getValue(GOVERNOR_R, &R, idx)) R = 0.05; 
  //if (!data->getValue(GOVERNOR_T1, &T1, idx)) T1 = 0.5; 
  //if (!data->getValue(GOVERNOR_VMAX, &Vmax, idx)) Vmax = 1.0; 
  //if (!data->getValue(GOVERNOR_VMIN, &Vmin, idx)) Vmin = 0.0; 
  //if (!data->getValue(GOVERNOR_T2, &T2, idx)) T2 = 3.0; 
  //if (!data->getValue(GOVERNOR_T3, &T3, idx)) T3 = 10.0; 
  //if (!data->getValue(GOVERNOR_DT, &Dt, idx)) Dt = 0.0;

  //// Set up transfer function blocks
  //if(integrationtype != IMPLICIT) {
  //  leadlag_blk.setparams(T2,T3);
  //  delay_blk.setparams(1.0,T1,Vmin,Vmax,-1000.0,1000.0);
  //}

  BaseEMTGovModel::load(data, idx); // load parameters in base governor model

  if (!data->getValue(GOVERNOR_R, &R, idx)) R = 0.05;
  if (!data->getValue(GOVERNOR_r, &r, idx)) R = 0.05;
  if (!data->getValue(GOVERNOR_TR, &TR, idx)) TR = 0.5;
  if (!data->getValue(GOVERNOR_TF, &TF, idx)) TF = 3.0;
  if (!data->getValue(GOVERNOR_TG, &TG, idx)) TG = 10.0;
  if (!data->getValue(GOVERNOR_VELM, &VELM, idx)) VELM = 1000.0;
  if (!data->getValue(GOVERNOR_GMAX, &GMAX, idx)) GMAX = 1.0;
  if (!data->getValue(GOVERNOR_GMIN, &GMIN, idx)) GMIN = 0.0;
  if (!data->getValue(GOVERNOR_TW, &TW, idx)) TW = 10.0;
  if (!data->getValue(GOVERNOR_AT, &AT, idx)) AT = 1.0;
  if (!data->getValue(GOVERNOR_DT, &Dt, idx)) Dt = 0.0;
  if (!data->getValue(GOVERNOR_QNL, &qNL, idx))  qNL = 0.0;

  filter_block.setparams(1.0, TF);
  gate_block.setparams(1 / r, 1 / (r * TR), GMIN, GMAX, -10000, 10000);
  opening_block.setparams(1.0, TG);
  turbine_flow_block.setparams(TW);
}

/**
   Number of variables
*/ 
void Hygov::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxgov = 0;

  *nvar = nxgov;
}
  
/**
   Prestep function
*/
void Hygov::preStep(double time ,double timestep)
{
  //if(integrationtype != EXPLICIT) return;

  //BaseEMTGenModel* gen=getGenerator();
  //double dw = gen->getSpeedDeviation();
  //
  //// Delay block output and state update
  //delay_blk_out = delay_blk.getoutput((Pref-dw)/R,timestep,true);

  //// Leadlag block output and state update
  //leadlag_blk_out = leadlag_blk.getoutput(delay_blk_out,timestep,true);

  //// Output mechanical power
  //Pmech = leadlag_blk_out - Dt*dw;

  //xout = Pmech;


  if(integrationtype != EXPLICIT) return;

  BaseEMTGenModel* gen=getGenerator();
  delta_w = gen->getSpeedDeviation();

  filter_block_in = nref - (delta_w + R * gate_block_out);

  double filter_block_out = filter_block.getoutput(filter_block_in, timestep, true);

  gate_block_out = gate_block.getoutput(filter_block_out, timestep, true);

  opening_block_out = opening_block.getoutput(gate_block_out, timestep, true);

  double turbine_flow_block_in, h;

  h = opening_block_out / turbine_flow_block_out;
  h = h * h;

  turbine_flow_block_in = 1 - h;

  turbine_flow_block_out = turbine_flow_block.getoutput(turbine_flow_block_in, timestep, true);

  Pmech = AT * (turbine_flow_block_out - qNL) * h - opening_block_out * Dt * delta_w;

}

/**
   Poststep function
*/
void Hygov::postStep(double time)
{
}


/**
 * Initialize governor model before calculation
 * @param [output] values - array where initialized governor variables should be set
 */
void Hygov::init(gridpack::RealType* xin) 
{
  //gridpack::RealType *x = xin+offsetb; // governor array starts from this location
  //
  //BaseEMTGenModel* gen=getGenerator();
  //double dw = gen->getSpeedDeviation();
  //double Pmech = gen->getInitialMechanicalPower();

  //if(integrationtype != IMPLICIT) {
  //  // Initialize leadlag block
  //  delay_blk_out = leadlag_blk.init_given_y(Pmech-Dt*dw);
  //  
  //  double delay_blk_in;
  //  
  //  // Initialize delay block
  //  delay_blk_in = delay_blk.init_given_y(delay_blk_out);
  //  
  //  // Reference power signal
  //  Pref = R*delay_blk_in + dw;

  //  return;
  //}

  //x1 = Pmech;
  //x2 = (1 - T2/T3)*x1;
  //xout = Pmech;

  //Pref = R*Pmech + dw;

  //if(integrationtype == IMPLICIT) {
  //  x[0] = x1;
  //  x[1] = x2;
  //  x[2] = xout;
  //}

  gridpack::RealType* x = xin + offsetb; // governor array starts from this location

  BaseEMTGenModel* gen = getGenerator();
  delta_w = gen->getSpeedDeviation();
  // delta_w = 0.0;
  Pmech = gen->getInitialMechanicalPower();

  if (integrationtype != IMPLICIT) {
      turbine_flow_block_out = Pmech / AT + qNL;

      double turbine_flow_block_in;

      turbine_flow_block_in = turbine_flow_block.init_given_y(turbine_flow_block_out);

      opening_block_out = turbine_flow_block_out;

      gate_block_out = opening_block.init_given_y(opening_block_out);

      double gate_block_in;
      gate_block_in = gate_block.init_given_y(gate_block_out);

      filter_block_in = filter_block.init_given_y(gate_block_in);

      nref = filter_block_in + R * gate_block_out + delta_w;

      return;
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
bool Hygov::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out governor state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Hygov::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto governors
 * @param values array containing governor state variables
*/
void Hygov::setValues(gridpack::RealType *val)
{
  //gridpack::RealType *values = val+offsetb; // governor array starts from this location

  //if(integrationtype == EXPLICIT) return;

  //if(p_mode == XVECTOBUS) {
  //  x1 = values[0];
  //  x2 = values[1];
  //  xout = values[2];
  //} else if(p_mode == XDOTVECTOBUS) {
  //  dx1 = values[0];
  //  dx2 = values[1];
  //}

  if (integrationtype != EXPLICIT) return;
}

/**
 * Return the values of the governor vector block
 * @param values: pointer to vector values
 */
void Hygov::vectorGetValues(gridpack::RealType *values)
{
  //gridpack::RealType *f = values+offsetb; // exciter array starts from this location

  //if(integrationtype == EXPLICIT) return;
  //
  //int x1_idx = 0;
  //int x2_idx = 1;
  //int xout_idx = 2;

  //double yLL;
  //BaseEMTGenModel* gen=getGenerator();
  //double dw = gen->getSpeedDeviation();

  //if(p_mode == RESIDUAL_EVAL) {
  //  // x1 equation
  //  f[x1_idx] = 1/T1*((1/R)*(Pref - dw) - x1) - dx1;
  //  f[x2_idx] = (-x2 + (1.0 - T2/T3)*x1)/T3 - dx2;

  //  Pmech = x2 + T2/T3*x1 - Dt*dw;

  //  f[xout_idx] = Pmech - xout;
  //}

    if (integrationtype != EXPLICIT) return;
}

/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values
 */
int Hygov::matrixNumValues()
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
void Hygov::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  //int ctr = 0;

  //if(integrationtype != IMPLICIT) {
  //  *nvals = 0;
  //  return;
  //}
  //
  //int x1_gloc = p_gloc;
  //int x2_gloc = p_gloc+1;
  //int xout_gloc = p_gloc+2;
  //int    dw_gloc;
  //double dw;

  //dw = getGenerator()->getSpeedDeviation(&dw_gloc);

  //// Partial derivatives of 1 equation
  //rows[ctr] = x1_gloc; cols[ctr] = x1_gloc;
  //rows[ctr+1] = x1_gloc; cols[ctr+1] = dw_gloc;
  //
  //values[ctr] = values[ctr+1] = 0.0;
  //
  //values[ctr] = -1.0/T1 - shift;
  //values[ctr+1] = -1/T1*(1/R);

  //ctr += 2;

  //rows[ctr]   = x2_gloc; cols[ctr] = x1_gloc;
  //rows[ctr+1] = x2_gloc; cols[ctr+1] = x2_gloc;

  //values[ctr] = values[ctr+1] = 0.0;

  //values[ctr] = (1.0 - T2/T3)/T3;
  //values[ctr+1] = -1.0/T3 - shift;

  //ctr += 2;

  //rows[ctr]   = xout_gloc; cols[ctr] = x1_gloc;
  //rows[ctr+1] = xout_gloc; cols[ctr+1] = x2_gloc;
  //rows[ctr+2] = xout_gloc; cols[ctr+2] = xout_gloc;
  //rows[ctr+3] = xout_gloc; cols[ctr+3] = dw_gloc;

  //values[ctr]   = T2/T3;
  //values[ctr+1] = 1.0;
  //values[ctr+2] = -1.0;
  //values[ctr+3] = -Dt;
  //
  //ctr += 4;

  //*nvals = ctr;


    if (integrationtype != EXPLICIT) return;
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power 
 */
void Hygov::setInitialMechanicalPower(double Pmech0)
{
  Pmech = Pmech0;
}

/** 
 * Get the value of the mechanical power parameter
 * @return value of the mechanical power 
 */
double Hygov::getMechanicalPower()
{
  return xout;
}

/** 
 * Get the value of the mechanical power and its global location
 * @return value of the mechanical power
 *
 * Note: Used in Jacobian calculation
 */
double Hygov::getMechanicalPower(int *Pmech_gloc)
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
bool Hygov::getMechanicalPowerPartialDerivatives(int *xgov_loc,double *dPmech_dxgov)
{
  return false;
}

void Hygov::setVcomp(double Vcomp)
{
}

/**
 * Update the event function values
 */
void Hygov::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues)
{
//  int offset   = getLocalOffset();
//  int x1_idx   = offset;
//
//  BaseEMTGenModel* gen=getGenerator();
//  double dw = gen->getSpeedDeviation();
//
//  x1  = state[x1_idx];
//
//  double dx1_dt = ((1/R)*(Pref - dw) - x1);
//;
//
//  /* Limits on x1 */
//  if(!x1_at_min) {
//    evalues[0] = x1 - Vmin;
//  } else {
//    evalues[0] = -dx1_dt; /* Release when derivative reaches 0 */
//  }
//
//  if(!x1_at_max) {
//    evalues[1] = Vmax - x1;
//  } else {
//    evalues[1] = dx1_dt; /* Release when derivative reaches 0 */
//  }

    if (integrationtype != EXPLICIT) return;
} 

/**
 * Event handler
 */
void Hygov::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  //int offset    = getLocalOffset();
  //int x1_idx = offset;

  //BaseEMTGenModel* gen=getGenerator();
  //double dw = gen->getSpeedDeviation();

  //x1  = state[x1_idx];

  //double dx1_dt = ((1/R)*(Pref - dw) - x1);

  //if(triggered[0]) {
  //  if(!x1_at_min && dx1_dt < 0) {
  //    /* Hold x1 at Vmin */
  //    x1_at_min = true;
  //  } else {
  //    /* Release */
  //    x1_at_min = false;
  //  }
  //}

  //if(triggered[1]) {
  //  if(!x1_at_max && dx1_dt > 0) {
  //    /* Hold x1 at Pmax */
  //    x1_at_max = true;
  //  } else {
  //    /* Release */
  //    x1_at_max = false;
  //  }
  //}

    if (integrationtype != EXPLICIT) return;
}

/**
 * Set event
 */
void Hygov::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  //if(integrationtype == IMPLICIT) {
  //  gridpack::math::RealDAESolver::EventPtr e(new HygovEvent(this));

  //  eman->add(e);
  //}

    if (integrationtype != EXPLICIT) return;
}

void HygovEvent::p_update(const double& t,gridpack::RealType *state)
{
  //p_gov->eventFunction(t,state,p_current);
}

void HygovEvent::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  //p_gov->eventHandlerFunction(triggered,t,state);
}
