/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wtpta1.cpp
 *  
 * @brief WTPTA1 pitch controller model implementation 
 *
 *
 */

#include <wtpta1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Wtpta1::Wtpta1(void)
{
  nxrmech = 0;
  Theta = 0.0;
}

Wtpta1::~Wtpta1(void)
{
}

void Wtpta1::getnvar(int *nvar)
{
  *nvar = nxrmech;
}

/**
 * Load parameters from DataCollection object into wtpta1 model
 * @param data collection of wtpta1 parameters from input files
 * @param index of wtpta1 on bus
 * TODO: might want to move this functionality to BaseExciterModel
 */
void Wtpta1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTRMechModel::load(data,idx); // load parameters in base wtpta1 model

  if (!data->getValue(WIND_PC_KIW, &Kiw, idx)) Kiw = 0.0; 
  if (!data->getValue(WIND_PC_KPW, &Kpw, idx)) Kpw = 0.0; 
  if (!data->getValue(WIND_PC_KIC, &Kic, idx)) Kic = 0.0; 
  if (!data->getValue(WIND_PC_KPC ,&Kpc,idx)) Kpc = 0.0;
  if (!data->getValue(WIND_PC_KCC, &Kcc, idx)) Kcc = 0.0;
  if (!data->getValue(WIND_PC_TP, &Tp, idx)) Tp = 0.0;
  if (!data->getValue(WIND_PC_THETAMAX,&Thetamax,idx))  Thetamax = 0.0; 
  if (!data->getValue(WIND_PC_THETAMIN,&Thetamin,idx)) Thetamin = 0.0;
  if (!data->getValue(WIND_PC_RTHETAMAX,&dThetamax,idx))  dThetamax = 0.0;
  if (!data->getValue(WIND_PC_RTHETAMIN,&dThetamin,idx)) dThetamin = 0.0;
  
  // Set parameters for blocks
  pitchcomp_blk.setparams(Kpc,Kic,Thetamin,Thetamax,-1000.0,1000.0);
  pitchctrl_blk.setparams(Kpw,Kiw,Thetamin,Thetamax,-1000.0,1000.0);
  lag_blk.setparams(1.0,Tp,Thetamin,Thetamax,dThetamin,dThetamax,-1000.0,1000.0);

}

/**
 * Initialize model before calculation
 * @param [xin] values - array where initialized wtpta1 variables should be set
 */
void Wtpta1::init(gridpack::RealType* xin) 
{
  gridpack::RealType *x = xin+offsetb; // wtpta1 array starts from this location

  double Thetacmd,u1,u2;
  /* Create string for setting name */
  std::string blkhead = std::to_string(busnum) + "_" + id + "WTPTA1_";

  std::string lag_block_name = blkhead + "lag_blk";
  lag_blk.setname(lag_block_name.c_str());

  Theta = getDriveTrainController()->getAeroDynamicController()->getTheta();

  Thetacmd = lag_blk.init_given_y(Theta);

  std::string pitchcomp_block_name = blkhead + "pitchcomp_blk";
  pitchcomp_blk.setname(pitchcomp_block_name.c_str());

  Pord = getElectricalController()->getPord();

  Pord0 = getPlantController()->getPref();

  domega_t = getDriveTrainController()->getTurbineSpeedDeviation();

  omega_ref = getTorqueController()->getOmegaref();
  
  u1 = pitchcomp_blk.init_given_y(0.0);

  std::string pitchctrl_block_name = blkhead + "pitchctrl_blk";
  pitchctrl_blk.setname(pitchctrl_block_name.c_str());

  u2 = pitchctrl_blk.init_given_y(Thetacmd);

}

/**
 * Write output from wtpta1s to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Wtpta1::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out wtpta1 state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Wtpta1::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto wtpta1s
 * @param values array containing wtpta1 state variables
*/
void Wtpta1::setValues(gridpack::RealType *val)
{
  gridpack::RealType *values = val+offsetb; // wtpta1 array starts from this location

  if(p_mode == XVECTOBUS) {
  } else if(p_mode == XDOTVECTOBUS) {
  }
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Wtpta1::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // wtpta1 array starts from this location
}

/**
   Prestep function
*/
void Wtpta1::preStep(double time ,double timestep)
{
  double y1,y2;
  
  if(integrationtype != EXPLICIT) return;

  Pord = getElectricalController()->getPord();

  Pord0 = getPlantController()->getPref();

  domega_t = getDriveTrainController()->getTurbineSpeedDeviation();

  omega_ref = getTorqueController()->getOmegaref();

  y1 = pitchcomp_blk.getoutput(Pord-Pord0,timestep,true);
  y2 = pitchctrl_blk.getoutput((1+domega_t)-omega_ref + Kcc*(Pord - Pord0),timestep,true);

  Theta = lag_blk.getoutput(y1+y2,timestep,true);

}

/**
   Poststep function
*/
void Wtpta1::postStep(double time)
{
}

/**
 * Get number of matrix values contributed by wtpta1
 * @return number of matrix values
 */
int Wtpta1::matrixNumValues()
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
void Wtpta1::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;
  *nvals = ctr;		     
}

/**
 * setTurbineSpeedDeviation - sets the turbine speed deviation
 * @param domega_turb : turbine speed deviation
 * From drive train model
 **/
void Wtpta1::setTurbineSpeedDeviation(double domega_turb)
{
  domega_t = domega_turb;
}

/**
 * setPord - sets Pord
 * @param Pord - electric Pord
 * From electrical controller model
 **/
void Wtpta1::setPord(double Pord_in)
{
  Pord = Pord_in;
}

/**
 *  setPord0 - Sets initial power order
 *  @param Pord0 - Initial value of reference power
 *
 **/
void Wtpta1::setPord0(double Pord0_in)
{
  Pord0 = Pord0_in;
}

/**
 *  setOmegaref - Output of torque controller
 * @param - omega_ref : reference speed
 **/
void Wtpta1::setOmegaref(double omega_ref_in)
{
  omega_ref = omega_ref_in;
}

/**
 * setTheta - Set initial value of  pitch controller
 * @param theta0 - initial value of pitch controller
 **/
void Wtpta1::setTheta(double Theta0)
{
  Theta = Theta0;
}

/**
 * getTheta - Get output of pitch controller
 * @output Theta - pitch angle (degrees) output of pitch controller
 **/
double Wtpta1::getTheta()
{
  return Theta;
}

/**
 * Update the event function values
 */
void Wtpta1::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues)
{
  int offset    = getLocalOffset();
} 

/**
 * Event handler
 */
void Wtpta1::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  int offset    = getLocalOffset();
}

/**
 * Set event
 */
void Wtpta1::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  gridpack::math::RealDAESolver::EventPtr e(new Wtpta1Event(this));

  eman->add(e);
}

void Wtpta1Event::p_update(const double& t,gridpack::RealType *state)
{
  p_pitchcon->eventFunction(t,state,p_current);
}

void Wtpta1Event::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  p_pitchcon->eventHandlerFunction(triggered,t,state);
}
