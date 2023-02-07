/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   epria1.cpp
 * @author Shrirang Abhyankar
 * @Added: February 7, 2023
 * 
 * @brief  EPRI IBR generator/converter model
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_generator_model.hpp"
#include "epria1.hpp"


/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Epria1Generator::Epria1Generator(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Epria1Generator::~Epria1Generator(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::Epria1Generator::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{

  Model_PrintInfo();
  
  data->getValue(BUS_NUMBER,&p_bus_num);
  data->getValue(CASE_SBASE,&p_sbase);
  data->getValue(GENERATOR_ID,&p_gen_id,idx);
  if (!data->getValue(GENERATOR_PG, &p_pg,idx)) p_pg = 0.0;
  if (!data->getValue(GENERATOR_QG, &p_qg,idx)) p_qg = 0.0;
  if (!data->getValue(GENERATOR_STAT, &p_status,idx)) p_status = 0;
  if (!data->getValue(GENERATOR_MBASE, &p_mbase, idx))  p_mbase = 100.0; // MBase
  if(fabs(p_mbase) < 1e-6) p_mbase = p_pg; // Set mbase = p_pg if it is not given in file

  if(!data->getValue(EPRIA1_PARAM1, &modelparams.Vbase, idx)) modelparams.Vbase = 0.0;
  if(!data->getValue(EPRIA1_PARAM2, &modelparams.Sbase, idx)) modelparams.Sbase = 0.0;
  if(!data->getValue(EPRIA1_PARAM3, &modelparams.Vdcbase, idx)) modelparams.Vdcbase = 0.0;
  if(!data->getValue(EPRIA1_PARAM4, &modelparams.KpI, idx)) modelparams.KpI = 0.0;
  if(!data->getValue(EPRIA1_PARAM5, &modelparams.KiI, idx)) modelparams.KiI = 0.0;
  if(!data->getValue(EPRIA1_PARAM6, &modelparams.KpPLL, idx)) modelparams.KpPLL = 0.0;
  if(!data->getValue(EPRIA1_PARAM7, &modelparams.KiPLL, idx)) modelparams.KiPLL = 0.0;
  if(!data->getValue(EPRIA1_PARAM8, &modelparams.KpP, idx)) modelparams.KpP = 0.0;
  if(!data->getValue(EPRIA1_PARAM9, &modelparams.KiP, idx)) modelparams.KiP = 0.0;
  if(!data->getValue(EPRIA1_PARAM10, &modelparams.KpQ, idx)) modelparams.KpQ = 0.0;
  if(!data->getValue(EPRIA1_PARAM11, &modelparams.KiQ, idx)) modelparams.KiQ = 0.0;
  if(!data->getValue(EPRIA1_PARAM12, &modelparams.Imax, idx)) modelparams.Imax = 0.0;
  if(!data->getValue(EPRIA1_PARAM13, &modelparams.PQflag, idx)) modelparams.PQflag = 0;
  if(!data->getValue(EPRIA1_PARAM14, &modelparams.Vdip, idx)) modelparams.Vdip = 0.0;
  if(!data->getValue(EPRIA1_PARAM15, &modelparams.Vup, idx)) modelparams.Vup = 0.0;
  if(!data->getValue(EPRIA1_PARAM16, &modelparams.Rchoke, idx)) modelparams.Rchoke = 0.0;
  if(!data->getValue(EPRIA1_PARAM17, &modelparams.Lchoke, idx)) modelparams.Lchoke = 0.0;
  if(!data->getValue(EPRIA1_PARAM18, &modelparams.Cfilt, idx)) modelparams.Cfilt = 0.0;
  if(!data->getValue(EPRIA1_PARAM19, &modelparams.Rdamp, idx)) modelparams.Rdamp = 0.0;

  model.Parameters      = (void*)&modelparams;
  model.ExternalInputs  = (void*)&modelinputs;
  model.ExternalOutputs = (void*)&modeloutputs;
  model.DoubleStates    = (double*)modelstates;

}

/**
 * Initialize generator model before calculation
 * @param Vm voltage magnitude
 * @param Va voltage angle
 */
void gridpack::dynamic_simulation::Epria1Generator::init(double Vm,
    double Va, double ts)
{
  theta = Va; // save to local variable for later use

  // Real and imaginary components of voltage
  VR = Vm*cos(Va);
  VI = Vm*sin(Va);

  // Conversion to MW and MVAr
  p_pg *= p_sbase;
  p_qg *= p_sbase;

  // Prepare model inputs
  modelinputs.Va = VR*modelparams.Vbase;
  modelinputs.Vb = VI*modelparams.Vbase;

  // Real and reactive power reference
  modelinputs.Pref = p_pg;
  modelinputs.Qref = p_qg;

  // Voltage magnitude reference
  modelinputs.Vref = Vm;

  int ok = Model_Initialize(&model);

}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType gridpack::dynamic_simulation::Epria1Generator::INorton()
{
  return p_INorton;
}


/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType gridpack::dynamic_simulation::Epria1Generator::NortonImpedence()
{
  return Y_a;
}


/**
 * Predict part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Epria1Generator::predictor_currentInjection(bool flag)
{
  p_INorton = INorton();
}


/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Epria1Generator::predictor(
    double t_inc, bool flag)
{

}

/**
 * Corrector part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Epria1Generator::corrector_currentInjection(bool flag)
{
  p_INorton = INorton();
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Epria1Generator::corrector(
    double t_inc, bool flag)
{

}

bool gridpack::dynamic_simulation::Epria1Generator::tripGenerator()
{
	return false;
}

/**
 * Set voltage on each generator
 */
void gridpack::dynamic_simulation::Epria1Generator::setVoltage(
    gridpack::ComplexType voltage)
{
  Vt    = abs(voltage);
  theta = atan2(imag(voltage), real(voltage));
  VR    = real(voltage);
  VI    = imag(voltage);
}

/**
 * Set frequency on each generator, frequency is perunit
 */
void gridpack::dynamic_simulation::Epria1Generator::setFreq(double dFreq)
{
   busfreq = dFreq;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
bool gridpack::dynamic_simulation::Epria1Generator::serialWrite(
    char* string, const int bufsize, const char *signal)
{
  bool ret = false;
  if (!strcmp(signal,"watch")) {
    if(getWatch()) {
      ret = true; 
    }
  } else if(!strcmp(signal,"watch_header")) {
    if(getWatch()) {
      char buf[128];
	ret = true;
    }
  }
  return ret;
}

/**
 * return a vector containing any generator values that are being
 * watched
 * @param vals vector of watched values
 */
void gridpack::dynamic_simulation::Epria1Generator::getWatchValues(
    std::vector<double> &vals)
{
}
