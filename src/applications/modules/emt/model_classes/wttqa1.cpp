/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wttqa1.cpp
 *  
 * @brief WTTQA1 pitch controller model implementation 
 *
 *
 */

#include <wttqa1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Wttqa1::Wttqa1(void)
{
  nxrmech = 0;
  domega_g = 0.0;
  omega_ref = 1.0;
  Tflag = 1;
}

Wttqa1::~Wttqa1(void)
{
}

void Wttqa1::getnvar(int *nvar)
{
  *nvar = nxrmech;
}

/**
 * Load parameters from DataCollection object into wttqa1 model
 * @param data collection of wttqa1 parameters from input files
 * @param index of wttqa1 on bus
 * TODO: might want to move this functionality to BaseExciterModel
 */
void Wttqa1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTRMechModel::load(data,idx); // load parameters in base wttqa1 model

    if (!data->getValue(WIND_TC_TFLAG, &Tflag, idx)) Tflag = 1;  
  if (!data->getValue(WIND_TC_TP,&Tp,idx)) Tp = 0.0;
  if (!data->getValue(WIND_TC_TWREF,&Twref,idx)) Twref = 0.0; //Twref
  if (!data->getValue(WIND_TC_KPP,&Kpp,idx)) Kpp = 0.0; // Kpp
  if (!data->getValue(WIND_TC_KIP,&Kip,idx)) Kip = 0.0; // Kip
  if (!data->getValue(WIND_TC_TEMAX,&Temax,idx)) Temax = 0.0; // Temax
  if (!data->getValue(WIND_TC_TEMIN,&Temin,idx)) Temin = 0.0; //Temin
  if (!data->getValue(WIND_TC_P1,&p1,idx)) p1 = 0.0; // p1
  if (!data->getValue(WIND_TC_SPD1,&spd1,idx)) spd1 = 1.0; // spd1
  if (!data->getValue(WIND_TC_P2,&p2,idx)) p2 = 2.0; // p2
  if (!data->getValue(WIND_TC_SPD2,&spd2,idx)) spd2 = 1.0; // spd2
  if (!data->getValue(WIND_TC_P3,&p3,idx)) p3 = 0.0; // p3
  if (!data->getValue(WIND_TC_SPD3,&spd3,idx)) spd3 = 0.0; //spd3
  if (!data->getValue(WIND_TC_P4,&p4,idx)) p4 = 0.0; // p4
  if (!data->getValue(WIND_TC_SPD4,&spd4,idx)) spd4 = 0.0; // spd4
  if(!data->getValue(GENERATOR_MBASE, &MBase, idx)) MBase = 100.0;

  if(!data->getValue(WIND_TC_TRATE,&Trate,idx)) Trate = MBase;
  if(fabs(Trate) < 1e-6) Trate = MBase;

  // Set parameters in the blocks
  Pelec_filter_blk.setparams(1.0,Tp);
  wref_filter_blk.setparams(1.0,Twref);
  Tref_pi_blk.setparams(Kpp,Kip,Temin,Temax,-1000.0,1000.0);

  double x[4],y[4];
  x[0] = p1; y[0] = spd1;
  x[1] = p2; y[1] = spd2;
  x[2] = p3; y[2] = spd3;
  x[3] = p4; y[3] = spd4;

  Pomega_blk.setparams(4,x,y);
}

/**
 * Initialize model before calculation
 * @param [xin] values - array where initialized wttqa1 variables should be set
 */
void Wttqa1::init(gridpack::RealType* xin) 
{
  gridpack::RealType *x = xin+offsetb; // wttqa1 array starts from this location

  double Qg;
  getGenerator()->getPower(p_time,&Pelec,&Qg);

  Pref = getPlantController()->getPref();

  
    /* Create string for setting name */
  std::string blkhead = std::to_string(busnum) + "_" + id + "WTTQA1_";

  std::string Pelec_filter_block_name = blkhead + "Pelec_filter_blk";
  Pelec_filter_blk.setname(Pelec_filter_block_name.c_str());
  // Initialize Pelec filter block
  Pelec_filter_blk_out = Pelec_filter_blk.init_given_u(Pelec);

  Pomega_blk_out = Pomega_blk.getoutput(Pelec_filter_blk_out);

  std::string wref_filter_name = blkhead + "wref_filter_name";
  wref_filter_blk.setname(wref_filter_name.c_str());

  omega_ref = wref_filter_blk.init_given_u(Pomega_blk_out);
  
  // omega_ref is the same as the generator speed that we want to
  // keep the turbine at. This is the generator speed we want
  // to drive at depending on the electric power input

  // speed deviation
  domega_g = omega_ref - 1.0; // note speed deviation is negative

  Pref = Pelec;

  std::string Tref_block_name = blkhead + "Tref_blk";
  Tref_pi_blk.setname(Tref_block_name.c_str());
  Tref_pi_blk.init_given_y(Pelec/(1+domega_g));
}

/**
 * Write output from wttqa1s to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Wttqa1::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out wttqa1 state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Wttqa1::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto wttqa1s
 * @param values array containing wttqa1 state variables
*/
void Wttqa1::setValues(gridpack::RealType *val)
{
  gridpack::RealType *values = val+offsetb; // wttqa1 array starts from this location

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
void Wttqa1::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // wttqa1 array starts from this location
}

/**
   Prestep function
*/
void Wttqa1::preStep(double time ,double timestep)
{
  double dtorque;
  bool   updatestate = !Vdip;
  
  Pelec_filter_blk_out = Pelec_filter_blk.getoutput(Pelec,timestep,true);

  Pomega_blk_out = Pomega_blk.getoutput(Pelec_filter_blk_out);
  omega_ref = wref_filter_blk.getoutput(Pomega_blk_out,timestep,true);

  if(Tflag == 1) {
    dtorque = (Pref0 - Pelec_filter_blk_out)/(1 + domega_g);
    Tref_pi_blk_out = Tref_pi_blk.getoutput(dtorque,timestep,updatestate);
  } else {
    Tref_pi_blk_out = Tref_pi_blk.getoutput(1+domega_g-omega_ref,timestep,updatestate);
  }

  Pref = Tref_pi_blk_out*(1 + domega_g);
}

/**
   Poststep function
*/
void Wttqa1::postStep(double time)
{
}

/**
 * Get number of matrix values contributed by wttqa1
 * @return number of matrix values
 */
int Wttqa1::matrixNumValues()
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
void Wttqa1::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;
  *nvals = ctr;		     
}

/**
 *  Set Pref0 - Sets reference power
 *  @param Pref0 - reference power
 *  Set by plant controller model
 **/
void Wttqa1::setPref0(double Pref0_in)
{
  Pref0 = Pref0_in;
}

/**
 *  setPelec - Electrical power Pelec
 *  @param Pelec - Electrical power input
 **/
void Wttqa1::setPelec(double Pg)
{
  Pelec = Pg;
}

/**
 *  setGeneratorSpeedDeviation - Set the speed deviation
 *  @param - domega : speed deviation
 *  From drive train model
 **/
void Wttqa1::setGeneratorSpeedDeviation(double domega_in)
{
  domega_g = domega_in;
}

/**
 * setVdip - Voltage dip flag
 * @param vdip - flag to indicate voltage dip
 * From elecrical controller
 **/
void Wttqa1::setVdip(bool Vdip_in)
{
  Vdip = Vdip_in;
}

/**
 * getPref - Output of torque controller
 * @param  - Pref : reference power
 **/
double Wttqa1::getPref()
{
  return Pref;
}

/**
 *  getOmegaref - Output of torque controller
 * @param - omega_ref : reference speed
 **/
double Wttqa1::getOmegaref()
{
  return omega_ref;
}

/**
 * Update the event function values
 */
void Wttqa1::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues)
{
  int offset    = getLocalOffset();
} 

/**
 * Event handler
 */
void Wttqa1::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  int offset    = getLocalOffset();
}

/**
 * Set event
 */
void Wttqa1::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  gridpack::math::RealDAESolver::EventPtr e(new Wttqa1Event(this));

  eman->add(e);
}

void Wttqa1Event::p_update(const double& t,gridpack::RealType *state)
{
  p_torquecon->eventFunction(t,state,p_current);
}

void Wttqa1Event::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  p_torquecon->eventHandlerFunction(triggered,t,state);
}
