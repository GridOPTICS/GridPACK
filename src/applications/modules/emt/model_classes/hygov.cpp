/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hygov.cpp
 *  
 * @brief HYGOV governor model implementation 
 *
 */

#include <hygov.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>
#include <algorithm>

Hygov::Hygov(void)
{
  Pmech = 1.0;
  nref = 1.0;

  nxgov = 0; // Number of variables
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
  BaseEMTGovModel::load(data,idx); // load parameters in base governor model

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
  if(integrationtype != EXPLICIT) return;

  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();

  filter_block_in = nref - (dw + R*gate_block_out);

  double filter_block_out = filter_block.getoutput(filter_block_in,timestep,true);

  gate_block_out = gate_block.getoutput(filter_block_out,timestep,true);

  opening_block_out = opening_block.getoutput(gate_block_out,timestep,true);

  double turbine_flow_block_in,h;

  h = opening_block_out/turbine_flow_block_out;
  h = h*h;

  turbine_flow_block_in = 1 - h;

  turbine_flow_block_out = turbine_flow_block.getoutput(turbine_flow_block_in,timestep,true);

  Pmech =  AT*(turbine_flow_block_out - qNL)*h - opening_block_out*Dt*dw;

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
  if(integrationtype == IMPLICIT) return;

  // Set up transfer function blocks
  /* Create string for setting name */
  std::string blkhead = std::to_string(busnum) + "_" + id + "HYGOV_";

  filter_block.setparams(1.0,TF);
  
  std::string gate_block_name = blkhead + "gate_blk";
  gate_block.setname(gate_block_name.c_str());

  gate_block.setparams(1/r,1/(r*TR),GMIN,GMAX,-10000,10000);
  opening_block.setparams(1.0,TG);
  turbine_flow_block.setparams(TW);

  // Initialization begins from here

  gridpack::RealType *x = xin+offsetb; // governor array starts from this location
  
  BaseEMTGenModel* gen=getGenerator();
  double dw = gen->getSpeedDeviation();
  double Pmech = gen->getInitialMechanicalPower();

  turbine_flow_block_out = Pmech/AT + qNL;

  double turbine_flow_block_in;

  turbine_flow_block_in = turbine_flow_block.init_given_y(turbine_flow_block_out);

  opening_block_out = turbine_flow_block_out;

  gate_block_out = opening_block.init_given_y(opening_block_out);

  double gate_block_in;
  gate_block_in = gate_block.init_given_y(gate_block_out);
  
  filter_block_in = filter_block.init_given_y(gate_block_in);

  nref = filter_block_in + R*gate_block_out + dw;

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
  gridpack::RealType *values = val+offsetb; // governor array starts from this location

  if(integrationtype == EXPLICIT) return;

}

/**
 * Return the values of the governor vector block
 * @param values: pointer to vector values
 */
void Hygov::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // exciter array starts from this location

  if(integrationtype == EXPLICIT) return;
  
}

/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values
 */
int Hygov::matrixNumValues()
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
void Hygov::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;

  *nvals = ctr;
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
  return Pmech;
}

/** 
 * Get the value of the mechanical power and its global location
 * @return value of the mechanical power
 *
 * Note: Used in Jacobian calculation
 */
double Hygov::getMechanicalPower(int *Pmech_gloc)
{
  if(integrationtype == IMPLICIT) *Pmech_gloc = -1;
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
} 

/**
 * Event handler
 */
void Hygov::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
}

/**
 * Set event
 */
void Hygov::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
}

void HygovEvent::p_update(const double& t,gridpack::RealType *state)
{
}

void HygovEvent::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
}
