/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wtara1.cpp
 *  
 * @brief WTARA1 Aerodynamic model 
 *
 *
 */

#include <wtara1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Wtara1::Wtara1(void)
{
  nxrmech = 0;
  Theta = 0.0;
  Pmech0 = 0.0;
  domega_t = 0.0;
  Taero = 0.0;
}

Wtara1::~Wtara1(void)
{
}

void Wtara1::getnvar(int *nvar)
{
  *nvar = nxrmech;
}

/**
 * Load parameters from DataCollection object into wtara1 model
 * @param data collection of wtara1 parameters from input files
 * @param index of wtara1 on bus
 * TODO: might want to move this functionality to BaseExciterModel
 */
void Wtara1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTRMechModel::load(data,idx); // load parameters in base wtara1 model

  if (!data->getValue(WIND_AD_KA, &Ka, idx)) Ka = 0.007; // Ka
  if (!data->getValue(WIND_AD_THETA, &Theta0, idx)) Theta0 = 0.0; // theta
}

/**
 * Initialize model before calculation
 * @param [xin] values - array where initialized wtara1 variables should be set
 */
void Wtara1::init(gridpack::RealType* xin) 
{
  gridpack::RealType *x = xin+offsetb; // wtara1 array starts from this location

  double Pmech;
  Pmech = Pmech0 - Ka * Theta * (Theta - Theta0);

  Taero = Pmech/(1 + domega_t);

}

/**
 * Write output from wtara1s to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Wtara1::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out wtara1 state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Wtara1::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto wtara1s
 * @param values array containing wtara1 state variables
*/
void Wtara1::setValues(gridpack::RealType *val)
{
  gridpack::RealType *values = val+offsetb; // wtara1 array starts from this location

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
void Wtara1::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // wtara1 array starts from this location
}

/**
   Prestep function
*/
void Wtara1::preStep(double time ,double timestep)
{
  if(integrationtype != EXPLICIT) return;
  
  double Pmech;
  Pmech = Pmech0 - Ka * Theta * (Theta - Theta0);

  Taero = Pmech/(1 + domega_t);
}

/**
   Poststep function
*/
void Wtara1::postStep(double time)
{
}

/**
 * Get number of matrix values contributed by wtara1
 * @return number of matrix values
 */
int Wtara1::matrixNumValues()
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
void Wtara1::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;
  *nvals = ctr;		     
}

/**
 * Set the value of initial mechanical power
 * @param: Pmech0 => Initial mechanical power
 */
void  Wtara1::setPmech(double Pmech0_in)
{
  Pmech0 = Pmech0_in;
}

/**
 * setTurbineSpeedDeviation - sets the turbine speed deviation
 * @param domega_turb : turbine speed deviation
 * From drive train model
 **/
void Wtara1::setTurbineSpeedDeviation(double domega_turb)
{
  domega_t = domega_turb;
}

/**
 * setTheta - sets pitch angle
 * @param Theta - pitch angle
 **/
void Wtara1::setTheta(double Theta_in)
{
  Theta = Theta_in;
}

/**
 * getTaero - returns the aero-dynamic torque
 * @output Taero - aero dynamic torque, output of aerodynamic model
 **/
double Wtara1::getTaero()
{
  return Taero;
}

/**
 * getTheta - Get initial value of pitch controller
 * @output theta0 - initial value of pitch controller
 **/
double Wtara1::getTheta() {
  return Theta0;
}

/**
 * Update the event function values
 */
void Wtara1::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues)
{
  int offset    = getLocalOffset();
} 

/**
 * Event handler
 */
void Wtara1::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  int offset    = getLocalOffset();
}

/**
 * Set event
 */
void Wtara1::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  gridpack::math::RealDAESolver::EventPtr e(new Wtara1Event(this));

  eman->add(e);
}

void Wtara1Event::p_update(const double& t,gridpack::RealType *state)
{
  p_aerocon->eventFunction(t,state,p_current);
}

void Wtara1Event::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  p_aerocon->eventHandlerFunction(triggered,t,state);
}
