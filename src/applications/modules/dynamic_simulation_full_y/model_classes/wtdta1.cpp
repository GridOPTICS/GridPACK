/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   wtdta1.cpp
 * @author Shuangshuang Jin 
 * @Created on:    Nov 15, 2022
 * 
 * @Last updated
 * Shrirang Abhyankar
 * Added all pieces
 *
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <string>
#include <cstdio>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "wtdta1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Wtdta1Model::Wtdta1Model(void)
{
  domega_g = 0.0;
  domega_t = 0.0;
  double_mass = false;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Wtdta1Model::~Wtdta1Model(void)
{
}

/**
 * Load parameters from DataCollection object into mechanical model
 * @param data collection of mechanical parameters from input files
 * @param index of mechanical on bus
 * TODO: might want to move this functionality to
 * Wtdta1Model
 */
void gridpack::dynamic_simulation::Wtdta1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  double pi = 4.0*atan(1.0);
  double omega0 = 2*pi*60.0;
  
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
    Kshaft = (2*Ht *Hg * 2*pi*Freq1)/(H*omega0);

    Tshaft_blk.setparams(1/Kshaft);
  }
}

/**
 * Initialize  model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Wtdta1Model::init(double mag, double ang, double ts)
{
  double u;

  dtheta_g = ang;
  u = dthetag_blk.init_given_y(dtheta_g);

  domega_g = s0; // speed deviation
  u = domegag_blk.init_given_y(domega_g);

  Tm = Te + D*domega_g + u;

  domega_t = domega_g;
  if(double_mass) {
    domegat_blk.init_given_y(domega_t);
    Tshaft_blk.init_given_y(0.0);
  }
    
}

void gridpack::dynamic_simulation::Wtdta1Model::computeModel(double t_inc, IntegrationStage int_flag)
{
  double u;
  double pi = 4.0*atan(1.0);
  double omega0 = 2*pi*60.0;

  u = Tm - Te - D*domega_g;
  domega_g = domegag_blk.getoutput(u,t_inc,int_flag,true);

  u = omega0*(domega_g - s0);

  dtheta_g = dthetag_blk.getoutput(u,t_inc,int_flag,true);

}
/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wtdta1Model::predictor(double t_inc, bool flag)
{
  computeModel(t_inc,PREDICTOR);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Wtdta1Model::corrector(double t_inc, bool flag)
{
  computeModel(t_inc,CORRECTOR);
}

/**
 * setTmech - sets mechanical torque
 * @param Tmech - mechanical torque
 * From aerodynamic model
 **/
void gridpack::dynamic_simulation::Wtdta1Model::setTmech(double Tmech_in)
{
  Tm = Tmech_in;
}
  
/**
 * setTelec - sets electrical torque
 * @param Telec - electrical torque
 * From generator
   **/
void gridpack::dynamic_simulation::Wtdta1Model::setTelec(double Telec_in)
{
  Te = Telec_in;
}

/**
 * getTmech - get initial mechanical torque
 * @return Tm - initial mechanical torque
   **/
double gridpack::dynamic_simulation::Wtdta1Model::getTmech()
{
  return Tm;
}


/**
 * getTurbineSpeedDeviation - gets the turbine speed deviation
 * @ouput turbine speed deviation
 **/
double gridpack::dynamic_simulation::Wtdta1Model::getTurbineSpeedDeviation()
{
  return domega_t;
}

/**
 * getGenSpeedDeviation - gets the generator speed deviation
 * @ouput generator speed deviation
 **/
double gridpack::dynamic_simulation::Wtdta1Model::getGeneratorSpeedDeviation()
{
  return domega_g;
}

/**
 * getGenRotorAngleDeviation - gets the rotor angle deviation
 * @ouput rotor angle deviation
 **/
double gridpack::dynamic_simulation::Wtdta1Model::getRotorAngleDeviation()
{
  return dtheta_g;
}

/**
 *  setOmegaref - Output of torque controller
 * @param - omega_ref : reference speed
 **/
void gridpack::dynamic_simulation::Wtdta1Model::setOmegaref(double omega_ref)
{
  s0 = omega_ref - 1.0;
}
