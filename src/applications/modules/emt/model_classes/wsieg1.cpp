/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wsieg1.cpp
 *  
 * @brief WSIEG1 governor model implementation 
 *
 * Model assumes there is no second generator
 * So only Pmech1 is active
 *
 */

#include <wsieg1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Wsieg1::Wsieg1(void)
{
  xLL = 0.0; 
  xGV = 0.0; 
  xT1 = 0.0; 
  xT2 = 0.0; 
  xT3 = 0.0; 
  xT4 = 0.0;
  xout = 0.0;

  dxLL = 0.0;
  dxGV = 0.0;
  dxT1 = 0.0;
  dxT2 = 0.0;
  dxT3 = 0.0;
  dxT4 = 0.0;

  K = 0.0;
  T1 = 0.0;
  T2 = 0.0;
  T3 = 0.0;
  Uo = 0.0;
  Uc = 0.0;
  Pmax = 0.0;
  Pmin = 0.0;
  T4 = 0.0;
  K1 = 0.0;
  K2 = 0.0;
  T5 = 0.0;
  K3 = 0.0;
  K4 = 0.0;
  T6 = 0.0;
  K5 = 0.0;
  K6 = 0.0;
  T7 = 0.0;
  K7 = 0.0;
  K8 = 0.0;
  SecondGenExists = false;
  uGV_at_min = uGV_at_max = false;
  xGV_at_min = xGV_at_max = false;

  nxgov = 7; // Number of variables
}

Wsieg1::~Wsieg1(void)
{
}

/**
 * Load parameters from DataCollection object into governor model
 * @param data collection of governor parameters from input files
 * @param index of governor on bus
 * TODO: might want to move this functionality to BaseGovernorModel
 */
void Wsieg1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGovModel::load(data,idx); // load parameters in base governor model
  
  // load parameters for the model type
  if (!data->getValue(GOVERNOR_K, &K, idx)) K = 0.0; // K
  if (!data->getValue(GOVERNOR_T1, &T1, idx)) T1 = 0.0; // T1
  if (!data->getValue(GOVERNOR_T2, &T2, idx)) T2 = 0.0; // T2
  if (!data->getValue(GOVERNOR_T3, &T3, idx)) T3 = 0.0; // T3
  if (!data->getValue(GOVERNOR_UO, &Uo, idx)) Uo = 0.0; // Uo
  if (!data->getValue(GOVERNOR_UC, &Uc, idx)) Uc = 0.0; // Uc
  if (!data->getValue(GOVERNOR_PMAX, &Pmax, idx)) Pmax = 0.0; // Pmax
  if (!data->getValue(GOVERNOR_PMIN, &Pmin, idx)) Pmin = 0.0; // Pmin
  if (!data->getValue(GOVERNOR_T4, &T4, idx)) T4 = 0.0; // T4
  if (!data->getValue(GOVERNOR_K1, &K1, idx)) K1 = 0.0; // K1
  if (!data->getValue(GOVERNOR_K2, &K2, idx)) K2 = 0.0; // K2
  if (!data->getValue(GOVERNOR_T5, &T5, idx)) T5 = 0.0; // T5
  if (!data->getValue(GOVERNOR_K3, &K3, idx)) K3 = 0.0; // K3
  if (!data->getValue(GOVERNOR_K4, &K4, idx)) K4 = 0.0; // K4
  if (!data->getValue(GOVERNOR_T6, &T6, idx)) T6 = 0.0; // T6
  if (!data->getValue(GOVERNOR_K5, &K5, idx)) K5 = 0.0; // K5
  if (!data->getValue(GOVERNOR_K6, &K6, idx)) K6 = 0.0; // K6
  if (!data->getValue(GOVERNOR_T7, &T7, idx)) T7 = 0.0; // T7
  if (!data->getValue(GOVERNOR_K7, &K7, idx)) K7 = 0.0; // K7
  if (!data->getValue(GOVERNOR_K8, &K8, idx)) K8 = 0.0; // K8
  if (!data->getValue(GOVERNOR_DB1, &Db1, idx)) Db1 = 0.0; // Db1
  if (!data->getValue(GOVERNOR_ERR, &Err, idx)) Err = 0.0; // Err
  if (!data->getValue(GOVERNOR_DB2, &Db2, idx)) Db2 = 0.0; // Db2
  if (!data->getValue(GOVERNOR_IBLOCK, &Iblock, idx)) Iblock = 0.0; // Iblock

  //Set flags for differential or algebraic equations
  iseq_diff[0] = (T1 == 0 || T2 == 0)?0:1;
  iseq_diff[1] = 1;
  iseq_diff[2] = (T4 == 0)?0:1;
  iseq_diff[3] = (T5 == 0)?0:1;
  iseq_diff[4] = (T6 == 0)?0:1;
  iseq_diff[5] = (T7 == 0)?0:1;
}

/**
 * Initialize governor model before calculation
 * @param [output] values - array where initialized governor variables should be set
 */
void Wsieg1::init(gridpack::RealType* xin) 
{
  gridpack::RealType *x = xin+offsetb; // governor array starts from this location
  
  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();
  double PGV;
  if (K1 + K3 + K5 + K7 > 0) 
    PGV = Pmech1 / (K1 + K3 + K5 + K7);
  else if (K2 + K4 + K6 + K8 > 0) 
    PGV = Pmech2 / (K2 + K4 + K6 + K8);
  else 
    PGV = 0;
  
  if (SecondGenExists && (Pmech2 != 0) && (K2 + K4 + K6 + K8 > 0) && (PGV != 0)) {
    double temp = Pmech2 / PGV * (K2 + K4 + K6 + K8);
    K2 = temp * K2;
    K4 = temp * K4;
    K6 = temp * K6;
    K8 = temp * K8;
  }
  xT4 = PGV;
  xT3 = PGV;
  xT2 = PGV;
  xT1 = PGV;
  xGV = PGV;

  Pref = PGV;

  // Note: (GV > Pmax) or (GV < Pmin) is an initial state violation
  if (Iblock == 1 && Pmin == 0) Pmin = PGV;
  if (Iblock == 2 && Pmax == 0) Pmax = PGV;
  if (Iblock == 3 && Pmin == 0) Pmin = PGV;
  if (Iblock == 3 && Pmax == 0) Pmax = PGV;
  if (iseq_diff[0]) xLL = K*dw * (1.0 - T2 / T1);
  else xLL = K*dw;

  xout = xT1 * K1 + xT2 * K3 + xT3 * K5 + xT4 * K7;

  x[0] = xLL;
  x[1] = xGV;
  x[2] = xT1;
  x[3] = xT2;
  x[4] = xT3;
  x[5] = xT4;
  x[6] = xout;
  
}

/**
 * Write output from governors to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Wsieg1::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out governor state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Wsieg1::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto governors
 * @param values array containing governor state variables
*/
void Wsieg1::setValues(gridpack::RealType *val)
{
  gridpack::RealType *values = val+offsetb; // governor array starts from this location

  if(p_mode == XVECTOBUS) {
    xLL = values[0];
    xGV = values[1];
    xT1 = values[2];
    xT2 = values[3];
    xT3 = values[4];
    xT4 = values[5];
    xout = values[6];
  } else if(p_mode == XDOTVECTOBUS) {
    dxLL = values[0];
    dxGV = values[1];
    dxT1 = values[2];
    dxT2 = values[3];
    dxT3 = values[4];
    dxT4 = values[5];
  }
}

/**
 * Return the values of the governor vector block
 * @param values: pointer to vector values
 */
void Wsieg1::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // exciter array starts from this location

  int x1_idx = 0;
  int x2_idx = 1;
  int x3_idx = 2;
  int x4_idx = 3;
  int x5_idx = 4;
  int x6_idx = 5;
  int xout_idx = 6;
  double yLL,uGV;
  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();

  if(p_mode == RESIDUAL_EVAL) {
    // xLL equation
    if(iseq_diff[x1_idx]) {
      f[x1_idx] = (-xLL + (1.0 - T2/T1)*K*dw)/T1 - dxLL;
      yLL = xLL + T2/T1*K*dw;
    } else {
      f[x1_idx] = -xLL + K*dw;
      yLL = xLL;
    }

    // xGV equation
    uGV = (Pref - xGV - yLL)/T3;
    if(uGV_at_min) uGV = Uc;
    else if(uGV_at_max) uGV = Uo;

    if(xGV_at_min) f[x2_idx] = xGV - Pmin;
    else if(xGV_at_max) f[x2_idx] = xGV - Pmax;
    else f[x2_idx] = uGV - dxGV;
    
    // xT1 equation
    if(iseq_diff[x3_idx]) f[x3_idx] = (xGV - xT1)/T4 - dxT1;
    else f[x3_idx] = xGV - xT1;

    // xT2 equation
    if(iseq_diff[x4_idx]) f[x4_idx] = (xT1 - xT2)/T5 - dxT2;
    else f[x4_idx] = xT1 - xT2;

    // xT3 equation
    if(iseq_diff[x5_idx]) f[x5_idx] = (xT2 - xT3)/T6 - dxT3;
    else f[x5_idx] = xT2 - xT3;

    // xT4 equation
    if(iseq_diff[x6_idx]) f[x6_idx] = (xT3 - xT4)/T7 - dxT4;
    else f[x6_idx] = xT3 - xT4;        
  }

  // xout equation
  f[xout_idx] = xT1 * K1 + xT2 * K3 + xT3 * K5 + xT4 * K7 - xout;

}

/**
 * Set Jacobian block
 * @param values a 2-d array of Jacobian block for the bus
 */
bool Wsieg1::setJacobian(gridpack::RealType **values)
{

  return true;
}

/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values
 */
int Wsieg1::matrixNumValues()
{
  return 24;
}

  /**
   * Return values from Jacobian matrix
   * @param nvals: number of values to be inserted
   * @param values: pointer to matrix block values
   * @param rows: pointer to matrix block rows
   * @param cols: pointer to matrix block cols
   */
void Wsieg1::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;

  int xLL_gloc = p_gloc;
  int xGV_gloc = p_gloc+1;
  int xT1_gloc = p_gloc+2;
  int xT2_gloc = p_gloc+3;
  int xT3_gloc = p_gloc+4;
  int xT4_gloc = p_gloc+5;
  int xout_gloc = p_gloc+6;
  double yLL;
  double dyLL_dxLL=0.0,dyLL_dxGV=0.0;
  double dyLL_dxT1=0.0,dyLL_dxT2=0.0;
  double dyLL_dxT3=0.0,dyLL_dxT4=0.0;
  double dyLL_ddw=0.0;
  int    dw_gloc;
  double dw;

  dw = getGenerator()->getSpeedDeviation(&dw_gloc);

  // Partial derivatives of xLL equation
  rows[ctr] = xLL_gloc; cols[ctr] = xLL_gloc;
  rows[ctr+1] = xLL_gloc; cols[ctr+1] = dw_gloc;
  
  values[ctr] = values[ctr+1] = 0.0;
  if(iseq_diff[0]) {
    values[ctr] = -1.0/T1 - shift;
    values[ctr+1] = (1.0 - T2/T1)*K/T1;
    dyLL_dxLL = 1.0;
    dyLL_ddw  = T2/T1*K;
  } else {
    values[ctr] = -1.0;
    values[ctr+1] = K;
    dyLL_dxLL = 1.0;
  }
  ctr += 2;

  rows[ctr]   = xGV_gloc; cols[ctr] = xLL_gloc;
  rows[ctr+1] = xGV_gloc; cols[ctr+1] = xGV_gloc;
  rows[ctr+2] = xGV_gloc; cols[ctr+2] = dw_gloc;

  values[ctr] = values[ctr+1] = values[ctr+2] = 0.0;
  values[ctr+3] = values[ctr+4] = values[ctr+5] = values[ctr+6] = 0.0;
  // Partial derivatives of xGV equation
  if(xGV_at_min || xGV_at_max) {
    values[ctr+1] = 1.0;
  } else {
    if(!uGV_at_min && !uGV_at_max) {
      values[ctr]   = -dyLL_dxLL/T3;
      values[ctr+1] = -1.0/T3 - shift;
      if(iseq_diff[0]) {
	values[ctr+2] = (-dyLL_ddw)/T2;
      }
    }
  }
  ctr += 3;

  // Partial derivatives of xT1 equation
  rows[ctr] = xT1_gloc; cols[ctr] = xGV_gloc;
  rows[ctr+1] = xT1_gloc; cols[ctr+1] = xT1_gloc;
  if(iseq_diff[2]) {
    values[ctr] = 1.0/T4;
    values[ctr+1] = -1.0/T4 - shift;
  } else {
    values[ctr] = 1.0;
    values[ctr+1] = -1.0;
  }
  ctr += 2;

  // Partial derivatives of xT2 equation
  rows[ctr] = xT2_gloc; cols[ctr] = xT1_gloc;
  rows[ctr+1] = xT2_gloc; cols[ctr+1] = xT2_gloc;
  if(iseq_diff[3]) {
    values[ctr] = 1.0/T5;
    values[ctr+1] = -1.0/T5 - shift;
  } else {
    values[ctr] = 1.0;
    values[ctr+1] = -1.0;
  }
  ctr += 2;

  // Partial derivatives of xT3 equation
  rows[ctr] = xT3_gloc; cols[ctr] = xT2_gloc;
  rows[ctr+1] = xT3_gloc; cols[ctr+1] = xT3_gloc;
  if(iseq_diff[4]) {
    values[ctr] = 1.0/T6;
    values[ctr+1] = -1.0/T6 - shift;
  } else {
    values[ctr] = 1.0;
    values[ctr+1] = -1.0;
  }
  ctr += 2;

  // Partial derivatives of xT4 equation
  rows[ctr] = xT4_gloc; cols[ctr] = xT3_gloc;
  rows[ctr+1] = xT4_gloc; cols[ctr+1] = xT4_gloc;
  if(iseq_diff[5]) {
    values[ctr] = 1.0/T7;
    values[ctr+1] = -1.0/T7 - shift;
  } else {
    values[ctr] = 1.0;
    values[ctr+1] = -1.0;
  }
  ctr += 2;

  rows[ctr] = xout_gloc;  cols[ctr] = xT1_gloc;
  rows[ctr+1] = xout_gloc; cols[ctr+1] = xT2_gloc;
  rows[ctr+2] = xout_gloc; cols[ctr+2] = xT3_gloc;
  rows[ctr+3] = xout_gloc; cols[ctr+3] = xT4_gloc;
  rows[ctr+4] = xout_gloc; cols[ctr+4] = xout_gloc;

  values[ctr]   = K1;
  values[ctr+1] = K2;
  values[ctr+2] = K3;
  values[ctr+3] = K4;
  values[ctr+4] = -1.0;

  ctr += 5;

  *nvals = ctr;
}

/**
 * Set the mechanical power parameter inside the governor
 * @param pmech value of the mechanical power 
 */
void Wsieg1::setInitialMechanicalPower(double Pmech0)
{
  Pmech1 = Pmech0;
  Pmech2 = Pmech0;
}

/** 
 * Get the value of the mechanical power parameter
 * @return value of the mechanical power 
 */
double Wsieg1::getMechanicalPower()
{
  Pmech1 = xout;

  return Pmech1;
}

/** 
 * Get the value of the mechanical power and its global location
 * @return value of the mechanical power
 *
 * Note: Used in Jacobian calculation
 */
double Wsieg1::getMechanicalPower(int *Pmech_gloc)
{
  *Pmech_gloc = p_gloc + 6;

  return getMechanicalPower();
}


/**
 * Partial derivatives of Mechanical Power Pmech w.r.t. governor variables
 * @param xgov_loc locations of governor variables
 * @param dPmech_dxgov partial derivatives of mechanical power Pmech w.r.t governor variables
*/
bool Wsieg1::getMechanicalPowerPartialDerivatives(int *xgov_loc,double *dPmech_dxgov)
{
  int i;

  for(i=0; i < nxgov; i++) xgov_loc[i] = offsetb + i;

  dPmech_dxgov[0] = dPmech_dxgov[1] = 0.0;
  dPmech_dxgov[2] = K1; 
  dPmech_dxgov[3] = K3; 
  dPmech_dxgov[4] = K5; 
  dPmech_dxgov[5] = K7; 

  return true;
}

void Wsieg1::setVcomp(double Vcomp)
{
}

/**
 * Update the event function values
 */
void Wsieg1::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues)
{
  int offset    = getLocalOffset();
  int xLL_idx = offset;
  int xGV_idx  = offset+1;
  int xT1_idx  = offset+2;
  int xT2_idx  = offset+3;
  int xT3_idx  = offset+4;
  int xT4_idx  = offset+5;

  double yLL;
  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();
  double uGV;

  xLL  = state[xLL_idx];
  xGV  = state[xGV_idx];
  xT1  = state[xT1_idx];
  xT2  = state[xT2_idx];
  xT3  = state[xT3_idx];
  xT4  = state[xT4_idx];

  if(iseq_diff[0]) {
    yLL = xLL + T2/T1*K*dw;
  } else {
    yLL = xLL;
  }

  uGV = (Pref - xGV - yLL)/T3;

  /* Limits on uGV */
  if(!uGV_at_min) {
    evalues[0] = uGV - Uc;
  } else {
    evalues[0] = Uc - uGV;
  }

  if(!uGV_at_max) {
    evalues[1] = Uo - uGV;
  } else {
    evalues[1] = uGV - Uo;
  }

  double dxGV_dt = uGV;

  /* Limits on xGV */
  if(!xGV_at_min) {
    evalues[2] = xGV - Pmin;
  } else {
    evalues[2] = -dxGV_dt; /* Release when derivative reaches 0 */
  }

  if(!xGV_at_max) {
    evalues[3] = Pmax - xGV;

  } else {
    evalues[3] = dxGV_dt; /* Release when derivative reaches 0 */
    //    printf("Va = %f, dVa_dt = %f\n",Va,dVa_dt);
  }
  
  //  printf("uGV = %f,xGV = %f,dxGV_dT = %f,Pmin = %f,Pmax = %f\n",uGV,xGV,dxGV_dt,Pmin,Pmax);
} 

/**
 * Reset limiter flags after a network resolve
 */
void Wsieg1::resetEventFlags()
{
  /* Note that the states are already pushed onto the network, so we can access these
     directly
  */
  double yLL;
  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();
  double uGV;

  if(iseq_diff[0]) {
    yLL = xLL + T2/T1*K*dw;
  } else {
    yLL = xLL;
  }

  uGV = (Pref - xGV - yLL)/T3;


  if(!uGV_at_min) {
    if(uGV - Uc < 0) uGV_at_min = true;
  } else {
    if(Uo - uGV < 0) uGV_at_min = false; /* Release */
  }

  if(!uGV_at_max) {
    if(Uo - uGV < 0) uGV_at_max = true;
  } else {
    if(uGV - Uo < 0) uGV_at_max = false; /* Release */
  }

  double dxGV_dt = uGV;

  if(!xGV_at_min) {
    if(xGV - Pmin < 0) xGV_at_min = true;
  } else {
    if(dxGV_dt > 0) xGV_at_min = false; /* Release */
  }

  if(!xGV_at_max) {
    if(Pmax - xGV < 0) xGV_at_max = true;
  } else {
    if(dxGV_dt < 0) xGV_at_max = false; /* Release */
  }
}

/**
 * Event handler
 */
void Wsieg1::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  int offset    = getLocalOffset();
  int xLL_idx = offset;
  int xGV_idx  = offset+1;
  int xT1_idx  = offset+2;
  int xT2_idx  = offset+3;
  int xT3_idx  = offset+4;
  int xT4_idx  = offset+5;

  double yLL;
  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();
  double uGV;

  xLL  = state[xLL_idx];
  xGV  = state[xGV_idx];
  xT1  = state[xT1_idx];
  xT2  = state[xT2_idx];
  xT3  = state[xT3_idx];
  xT4  = state[xT4_idx];

  if(iseq_diff[0]) {
    yLL = xLL + T2/T1*K*dw;
  } else {
    yLL = xLL;
  }

  uGV = (Pref - xGV - yLL)/T3;

  if(triggered[0]) {
    if(!uGV_at_min) {
      /* Hold uGV at Uc */
      uGV_at_min = true;
    } else {
      /* Release */
      uGV_at_max = false;
    }
  }

  if(triggered[1]) {
    if(!uGV_at_max) {
      /* Hold uGV at Uo */
      uGV_at_max = true;
    } else {
      /* Release */
      uGV_at_max = false;
    }
  }

  double dxGV_dt = uGV;
  if(triggered[2]) {
    if(!xGV_at_min && dxGV_dt < 0) {
      /* Hold xGV at Pmin */
      xGV_at_min = true;
    } else {
      /* Release */
      xGV_at_min = false;
    }
  }

  if(triggered[3]) {
    if(!xGV_at_max && dxGV_dt > 0) {
      /* Hold xGV at Pmax */
      xGV_at_max = true;
    } else {
      /* Release */
      xGV_at_max = false;
    }
  }
}

/**
 * Set event
 */
void Wsieg1::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  gridpack::math::RealDAESolver::EventPtr e(new Wsieg1Event(this));

  eman->add(e);
}

void Wsieg1Event::p_update(const double& t,gridpack::RealType *state)
{
  p_gov->eventFunction(t,state,p_current);
}

void Wsieg1Event::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  p_gov->eventHandlerFunction(triggered,t,state);
}
