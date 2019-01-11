/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   base_generator_model.hpp
 * @author Bruce Palmer
 * @Last modified:   May 18, 2015
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_generator_model.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "classical.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::ClassicalGenerator::ClassicalGenerator(void)
{
  p_PI = 4.0*atan(1.0);
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::ClassicalGenerator::~ClassicalGenerator(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::ClassicalGenerator::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  p_sbase = 100.0;

  int gstatus;
  double pg, qg, mva, r, dstr, dtr;
  double h, d0;

  data->getValue(BUS_NUMBER,&p_bus_id);
  data->getValue(GENERATOR_ID,&p_ckt,idx);
  if (!data->getValue(GENERATOR_PG, &p_pg,idx)) p_pg = 0.0;
  if (!data->getValue(GENERATOR_QG, &p_qg,idx)) p_qg = 0.0;
  if (!data->getValue(GENERATOR_STAT, &p_status,idx)) p_status = 0;
  p_pg /= p_sbase;
  p_qg /= p_sbase;

  if (!data->getValue(GENERATOR_MBASE, &p_mva, idx)) p_mva = 0.0;
  if (!data->getValue(GENERATOR_RESISTANCE, &p_r, idx)) p_r=0.0; // r
  if (!data->getValue(GENERATOR_SUBTRANSIENT_REACTANCE, &p_dstr,idx)) p_dstr=0.0; // dstr
  gridpack::ComplexType zsrc;
  if (!data->getValue(GENERATOR_ZSOURCE, &zsrc, idx))
    zsrc=gridpack::ComplexType(0.0,0.0); // dtr
  p_dtr = imag(zsrc);
  if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H, &p_h, idx)) p_h = 0.0; // h
  if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0, &p_d0, idx)) p_d0 = 0.0; // d0
}

/**
 * Initialize generator model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::ClassicalGenerator::init(double mag,
    double ang, double ts) 
{ 
  p_pg *= p_sbase; // p_pg *= p_sbase if ds2 read network from powerflow solution!
  p_qg *= p_sbase; // p_qg *= p_sbase if ds2 read network from powerflow solution!
  p_mva = p_sbase / p_mva;
  p_d0 = p_d0 / p_mva;
  p_h = p_h / p_mva;
  p_dtr = p_dtr * p_mva;
  p_pelect = p_pg;
  double eterm = mag;
  double vi = ang;
  gridpack::ComplexType v(0.0, vi);
  v = eterm * exp(v);
  p_volt = v;
  double pelect = p_pg;
  double qelect = p_qg;
  double currr = sqrt(pelect * pelect + qelect * qelect) / eterm;
  double phi = atan2(qelect, pelect);
  double curri = ang - phi;
  gridpack::ComplexType curr(0.0, curri);
  curr = currr * exp(curr);
  gridpack::ComplexType jay(0.0, 1.0);
  gridpack::ComplexType temp = v + jay * p_dtr * curr;
  p_eprime_s0 = temp;
  temp = atan2(imag(p_eprime_s0), real(p_eprime_s0));
  p_mac_ang_s0 = temp;
  p_mac_spd_s0 = gridpack::ComplexType(1.0,0.0);
  p_eqprime = gridpack::ComplexType(abs(p_eprime_s0),0.0);
  p_pmech = gridpack::ComplexType(abs(p_pelect),0.0);
  // Initialize other variables 
  p_mac_ang_s1 = gridpack::ComplexType(0.0,0.0);
  p_mac_spd_s1 = gridpack::ComplexType(0.0,0.0);
  p_dmac_ang_s0 = gridpack::ComplexType(0.0,0.0);
  p_dmac_spd_s0 = gridpack::ComplexType(0.0,0.0);
  p_dmac_ang_s1 = gridpack::ComplexType(0.0,0.0);
  p_dmac_spd_s1 = gridpack::ComplexType(0.0,0.0);
  p_eprime_s1 = gridpack::ComplexType(0.0,0.0);
  p_INorton = gridpack::ComplexType(0.0,0.0);
}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType gridpack::dynamic_simulation::ClassicalGenerator::INorton()
{
  return p_INorton;
}


/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType gridpack::dynamic_simulation::ClassicalGenerator::NortonImpedence()
{
  double ra = p_r * p_sbase / p_mva;
  double xd;
  if (p_dstr == 0.0) {
    xd = p_dtr * p_sbase / p_mva;
  }
  gridpack::ComplexType Y_a(ra, xd);
  Y_a = 1.0 / Y_a;
  return Y_a;
}

/**
 * Predict part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::ClassicalGenerator::predictor_currentInjection(bool flag)
{
  if (!flag) {
    p_mac_ang_s0 = p_mac_ang_s1;
    p_mac_spd_s0 = p_mac_spd_s1;
    p_eprime_s0 = p_eprime_s1;
  }
  
  gridpack::ComplexType jay(0.0, 1.0);
  // Calculate INorton_full
  p_INorton = p_eprime_s0 / (p_dtr * jay);
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::ClassicalGenerator::predictor(
    double t_inc, bool flag)
{
  gridpack::ComplexType jay, curr;
  double sysFreq = 60.0;
  const double basrad = 2.0*p_PI*sysFreq;
  // Reset values of machine parameters after first time step
  if (!flag) {
    p_mac_ang_s0 = p_mac_ang_s1;
    p_mac_spd_s0 = p_mac_spd_s1;
    p_eprime_s0 = p_eprime_s1;
  }
  // Evaluate updated values of machine parameters for integration
  jay = gridpack::ComplexType(0.0,1.0);

  // terminal curr: curr = (eprime - volt) / GEN_dtr) -----------
  curr = p_eprime_s0;
  curr -= p_volt;
  curr = curr/(jay*p_dtr);
  curr = conj(curr);
  p_pelect = real(p_eprime_s0 * curr);
  // dmac_ang:
  p_dmac_ang_s0 = (p_mac_spd_s0 - 1.0) * basrad;
  // dmac_spd:
  p_dmac_spd_s0 = ((p_pmech - p_d0 * (p_mac_spd_s0 - 1.0)) / p_mac_spd_s0 - p_pelect) / (2.0 * p_h);
  p_mac_ang_s1 = p_mac_ang_s0 + p_dmac_ang_s0 * t_inc;
  p_mac_spd_s1 = p_mac_spd_s0 + p_dmac_spd_s0 * t_inc;

  p_eprime_s1 = exp(p_mac_ang_s1 * jay) * p_eqprime;
}

/**
 * Correct part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::ClassicalGenerator::corrector_currentInjection(bool flag)
{
  gridpack::ComplexType jay(0.0, 1.0);
  p_INorton = p_eprime_s1 / (p_dtr * jay);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::ClassicalGenerator::corrector(
    double t_inc, bool flag)
{
  gridpack::ComplexType jay, curr;
  double sysFreq = 60.0;
  const double basrad = 2.0*p_PI*sysFreq;
  // Evaluate updated values of machine parameters for integration
  jay = gridpack::ComplexType(0.0,1.0);
  p_INorton = gridpack::ComplexType(0.0,0.0);
  // terminal curr: curr = (eprime - volt) / GEN_dtr) -----------
  curr = p_eprime_s1;
  curr -= p_volt;
  curr = curr/(jay*p_dtr);
  curr = conj(curr);
  p_pelect = real(p_eprime_s1 * curr);
  // dmac_ang:
  p_dmac_ang_s1 = (p_mac_spd_s1 - 1.0) * basrad;
  // dmac_spd:
  p_dmac_spd_s1 = ((p_pmech - p_d0 * (p_mac_spd_s1 - 1.0)) / p_mac_spd_s1 - p_pelect) / (2.0 * p_h);

  p_mac_ang_s1 = p_mac_ang_s0 + (p_dmac_ang_s0 + p_dmac_ang_s1)
    / 2.0 * t_inc;
  p_mac_spd_s1 = p_mac_spd_s0 + (p_dmac_spd_s0 +
      p_dmac_spd_s1) / 2.0 * t_inc;

  p_eprime_s1 = exp(p_mac_ang_s1 * jay) * p_eqprime;
}

/**
 * Set voltage on each generator
 */
void gridpack::dynamic_simulation::ClassicalGenerator::setVoltage(
    gridpack::ComplexType voltage)
{
  p_volt = voltage;
}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool gridpack::dynamic_simulation::ClassicalGenerator::
serialWrite(char *string, const int bufsize, const char *signal)
{
  bool ret = false;
  if (!strcmp(signal,"watch_header")) {
    if (getWatch()) {
      char buf[128];
      std::string tag;
      if (p_ckt[0] != ' ') {
        tag = p_ckt;
      } else {
        tag = p_ckt[1];
      }
      sprintf(buf,", %d_%s_angle, %d_%s_speed",p_bus_id,tag.c_str(),
          p_bus_id,tag.c_str());
      if (strlen(buf) <= bufsize) {
        sprintf(string,"%s",buf);
        ret = true;
      } else {
        ret = false;
      }
    } else {
      ret = false;
    }
  } else if (!strcmp(signal,"watch")) {
    if (getWatch()) {
      char buf[128];
      sprintf(buf,", %f, %f",real(p_mac_ang_s1),real(p_mac_spd_s1));
      if (strlen(buf) <= bufsize) {
        sprintf(string,"%s",buf);
        ret = true;
      } else {
        ret = false;
      }
    } else {
      ret = false;
    }
  } else if (!strcmp(signal,"debug_initial")) {
    sprintf(string,"PG: %f QG: %f\n",p_pg,p_qg);
    ret = true;
  }
  return ret;
}

/**
 * Get the roter angle
 */
double gridpack::dynamic_simulation::ClassicalGenerator::getAngle()
{
  return real(p_mac_ang_s1);
}
  
/**
 * return a vector containing any generator values that are being
 * watched
 * @param vals vector of watched values
 */
void gridpack::dynamic_simulation::ClassicalGenerator::getWatchValues(
    std::vector<double> &vals)
{
  vals.clear();
  if (getWatch()) {
    vals.push_back(real(p_mac_ang_s1));
    vals.push_back(real(p_mac_spd_s1));
  }
}
