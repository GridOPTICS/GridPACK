/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wtdta1.cpp
 *  
 * @brief WTDTA1 model implementation 
 *
 *
 */

#include <wtdta1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Wtdta1::Wtdta1(void)
{
  nxrmech = 0;
  domega_g = 0.0;
  domega_t = 0.0;
  double_mass = false;

}

Wtdta1::~Wtdta1(void)
{
}

void Wtdta1::getnvar(int *nvar)
{
  *nvar = nxrmech;
}

/**
 * Load parameters from DataCollection object into wtdta1 model
 * @param data collection of wtdta1 parameters from input files
 * @param index of wtdta1 on bus
 * TODO: might want to move this functionality to BaseExciterModel
 */
void Wtdta1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTRMechModel::load(data,idx); // load parameters in base wtdta1 model

  if (!data->getValue(WIND_DT_H, &H, idx)) H = 0.0;
  if (!data->getValue(WIND_DT_DAMP,&D,idx)) D = 0.0;
  if (!data->getValue(WIND_DT_HFRAC,&Hfrac,idx)) Hfrac = 0.0;
  if (!data->getValue(WIND_DT_FREQ1,&Freq1,idx)) Freq1 = 0.0;
  if (!data->getValue(WIND_DT_DSHAFT,&Dshaft,idx)) Dshaft = 0.0;
  
  Ht = Hfrac*H;
  Hg = H - Ht;
  if(fabs(Hfrac) > 1e-6) double_mass = true;

  // Set parameters for blocks
  domegag_blk.setparams(2*Hg);
  dthetag_blk.setparams(1.0);

  if(double_mass) {
    domegat_blk.setparams(2*Ht);
    Kshaft = (2*Ht *Hg * 2*PI*Freq1)/(H*OMEGA_S);

    Tshaft_blk.setparams(1/Kshaft);
  }
}

/**
 * Initialize model before calculation
 * @param [xin] values - array where initialized wtdta1 variables should be set
 */
void Wtdta1::init(gridpack::RealType* xin) 
{
  gridpack::RealType *x = xin+offsetb; // wtdta1 array starts from this location
  double ang = getGenerator()->getAngle();

  /* Create string for setting name */
  std::string blkhead = std::to_string(busnum) + "_" + id + "WTDTA1_";

  double u;

  dtheta_g = ang;
  u = dthetag_blk.init_given_y(dtheta_g);

  s0 = 0.0;

  domega_g = s0; // speed deviation
  std::string domegag_block_name = blkhead + "domegag_blk";
  domegag_blk.setname(domegag_block_name.c_str());

  u = domegag_blk.init_given_y(domega_g);

  double pg,qg;
  getGenerator()->getPower(p_time,&pg,&qg);
  Te = pg;

  Tm = Te + D*domega_g + u;

  domega_t = domega_g;
  if(double_mass) {
    domegat_blk.init_given_y(domega_t);
    Tshaft_blk.init_given_y(0.0);
  }
}

/**
 * Write output from wtdta1s to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Wtdta1::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out wtdta1 state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Wtdta1::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto wtdta1s
 * @param values array containing wtdta1 state variables
*/
void Wtdta1::setValues(gridpack::RealType *val)
{
  gridpack::RealType *values = val+offsetb; // wtdta1 array starts from this location

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
void Wtdta1::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // wtdta1 array starts from this location
}

/**
   Prestep function
*/
void Wtdta1::preStep(double time ,double timestep)
{
  double u;
  
  if(integrationtype != EXPLICIT) return;

  double pg,qg;
  getGenerator()->getPower(time,&pg,&qg);
  Te = pg/(1+domega_g);

  u = Tm - Te - D*domega_g;
  domega_g = domegag_blk.getoutput(u,timestep,true);

  u = OMEGA_S*(domega_g - s0);

  dtheta_g = dthetag_blk.getoutput(u,timestep,true);
}

/**
   Poststep function
*/
void Wtdta1::postStep(double time)
{
}

/**
 * Get number of matrix values contributed by wtdta1
 * @return number of matrix values
 */
int Wtdta1::matrixNumValues()
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
void Wtdta1::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;
  *nvals = ctr;		     
}


/**
 * setTmech - sets mechanical torque
 * @param Tmech - mechanical torque
 * From aerodynamic model
 **/
void Wtdta1::setTmech(double Tmech_in)
{
  Tm = Tmech_in;
}
  
/**
 * setTelec - sets electrical torque
 * @param Telec - electrical torque
 * From generator
   **/
void Wtdta1::setTelec(double Telec_in)
{
  Te = Telec_in;
}

/**
 * getTmech - get initial mechanical torque
 * @return Tm - initial mechanical torque
   **/
double Wtdta1::getTmech()
{
  return Tm;
}


/**
 * getTurbineSpeedDeviation - gets the turbine speed deviation
 * @ouput turbine speed deviation
 **/
double Wtdta1::getTurbineSpeedDeviation()
{
  return domega_t;
}

/**
 * getGenSpeedDeviation - gets the generator speed deviation
 * @ouput generator speed deviation
 **/
double Wtdta1::getGeneratorSpeedDeviation()
{
  return domega_g;
}

/**
 * getGenRotorAngleDeviation - gets the rotor angle deviation
 * @ouput rotor angle deviation
 **/
double Wtdta1::getRotorAngleDeviation()
{
  return dtheta_g;
}

/**
 *  setOmegaref - Output of torque controller
 * @param - omega_ref : reference speed
 **/
void Wtdta1::setOmegaref(double omega_ref)
{
  s0 = omega_ref - 1.0;

}


/**
 * Update the event function values
 */
void Wtdta1::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues)
{
  int offset    = getLocalOffset();
} 

/**
 * Event handler
 */
void Wtdta1::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  int offset    = getLocalOffset();
}

/**
 * Set event
 */
void Wtdta1::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  gridpack::math::RealDAESolver::EventPtr e(new Wtdta1Event(this));

  eman->add(e);
}

void Wtdta1Event::p_update(const double& t,gridpack::RealType *state)
{
  p_dtrain->eventFunction(t,state,p_current);
}

void Wtdta1Event::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  p_dtrain->eventHandlerFunction(triggered,t,state);
}
