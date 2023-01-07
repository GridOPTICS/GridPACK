/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   gdform.cpp
 * @author Shrirang Abhyankar
 * @Last modified:   Jan. 6, 2023
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
#include "gridpack/parser/dictionary.hpp"
#include "base_generator_model.hpp"
#include "gdform.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::GridFormingGenerator::GridFormingGenerator(void)
{
	double pi = 4.0*atan(1.0);
	omega0 = 2*pi*60.0;
	zero_Tf = false;
	p_tripped = false;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::GridFormingGenerator::~GridFormingGenerator(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::GridFormingGenerator::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  p_sbase = 100.0;

  data->getValue(BUS_NUMBER,&p_bus_num);
  data->getValue(CASE_SBASE,&p_sbase);
  data->getValue(GENERATOR_ID,&p_gen_id,idx);

  if (!data->getValue(GENERATOR_PG, &p_pg,idx)) p_pg = 0.0;
  if (!data->getValue(GENERATOR_QG, &p_qg,idx)) p_qg = 0.0;
  if (!data->getValue(GENERATOR_STAT, &p_status,idx)) p_status = 0;
  
  if (!data->getValue(GENERATOR_MBASE, &p_mbase, idx)) p_mbase = 1000.0; // Mbase

  if(fabs(p_mbase) < 1e-6) p_mbase = p_pg; // Set mbase = p_pg if it is not given in file

  if (!data->getValue(GENERATOR_ZSOURCE, &Zsource, idx)) Zsource = gridpack::ComplexType(0.01,0.02);
  if (!data->getValue(GENERATOR_XL  , &XL, idx)) XL = 0.075; 
  if (!data->getValue(GENERATOR_TF, &Tf, idx)) Tf = 0.01666; //
  if(fabs(Tf) < 1e-6) zero_Tf = true;
  if (!data->getValue(GENERATOR_MQ, &mq, idx)) mq=0.05; // 
  if (!data->getValue(GENERATOR_KPV, &kpv, idx)) kpv=0.0; // 
  if (!data->getValue(GENERATOR_KIV, &kiv, idx)) kiv=5.86; // 
  if (!data->getValue(GENERATOR_EMAX, &Emax, idx)) Emax = 2.0; // 
  if (!data->getValue(GENERATOR_EMIN, &Emin, idx)) Emin = -2.0; //

  Edroop_min = Emin;
  Edroop_max = Emax;
  if (!data->getValue(GENERATOR_MP, &mp, idx)) mp=3.77; // 
  if (!data->getValue(GENERATOR_KPPMAX, &kppmax, idx)) kppmax=3.0; // 
  if (!data->getValue(GENERATOR_KIPMAX, &kipmax, idx)) kipmax=30.0; // 
  if (!data->getValue(GENERATOR_PSET, &Pset, idx)) Pset=0.44; // 
  if (!data->getValue(GENERATOR_PMAX, &Pmax, idx)) Pmax=1.0; // 
  if (!data->getValue(GENERATOR_PMIN, &Pmin, idx)) Pmin=0.0; // 
  if (!data->getValue(GENERATOR_IMAX, &Imax, idx)) Imax=2.5;
  if (!data->getValue(GENERATOR_WMAX, &wmax, idx)) wmax=0.4;
  if (!data->getValue(GENERATOR_WMIN, &wmin, idx)) wmin=-0.4;

  Ra = real(Zsource);
  XL = imag(Zsource);
  
  mp_org = mp;
  mq_org = mq;

  // Initialize blocks
  if(!zero_Tf) {
    P_filter_blk.setparams(1.0,Tf);
    Q_filter_blk.setparams(1.0,Tf);
    V_filter_blk.setparams(1.0,Tf);
  }

  Edroop_PI_blk.setparams(kpv,kiv);

  Pmax_PI_blk.setparams(kppmax,kipmax,-1000.0,0.0,-1000.0,0.0);
  Pmin_PI_blk.setparams(kppmax,kipmax,0.0,1000.0,0.0,1000.0);

  Delta_blk.setparams(1.0);
}

/**
 * Initialize generator model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 */
void gridpack::dynamic_simulation::GridFormingGenerator::init(double Vm,
    double Va, double ts)
{
  theta = Va;

  // Change of base from system MVA base to Machine base
  p_pg *= p_sbase/p_mbase;
  p_qg *= p_sbase/p_mbase;

  // Real and imaginary components of voltage
  VR = Vm*cos(Va);
  VI = Vm*sin(Va);

  V = gridpack::ComplexType(VR,VI);
  S = gridpack::ComplexType(p_pg,p_qg);
  I = conj(S/V);

  Im = abs(I);
  
  E = V + Zsource*I;

  Edroop = abs(E);
  delta  = arg(E);

  // Initialize blocks
  Edroop_PI_blk.init_given_y(Edroop);
  Delta_blk.init_given_y(delta);

  Pmax_PI_blk.init_given_y(0.0);
  Pmin_PI_blk.init_given_y(0.0);

  if(!zero_Tf) {
    P_filter_blk.init_given_y(p_pg);
    Q_filter_blk.init_given_y(p_qg);
    V_filter_blk.init_given_y(Vm);
  }
  
  Vset = Vm + p_qg*mq;
  Pset = p_pg;
  
  mp_org = mp;
  mq_org = mq;
  Vset_org = Vset;
  Pset_org = Pset;
  
  p_Norton_Ya = NortonImpedence();
}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType gridpack::dynamic_simulation::GridFormingGenerator::INorton()
{
  I = (E - V)/Zsource; // output current

  Im = abs(I);

  if(Im > Imax) {
    double Ir,Ii,Ia;
    Ia = arg(I);
    Ir = Imax*cos(Ia);
    Ii = Imax*sin(Ia);
    I = gridpack::ComplexType(Ir,Ii);
    
    E = V + Zsource*I;
    Edroop_max = abs(E);
    Edroop_min = 0.0;
  } else {
    Edroop_max = Emax;
    Edroop_min = Emin;
  }

  gridpack::ComplexType Ie;

  Ie = I + V/Zsource; // Current from internal voltage source
  
  p_INorton = Ie*p_mbase/p_sbase;
  
  return p_INorton;
}


/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType gridpack::dynamic_simulation::GridFormingGenerator::NortonImpedence()
{
  double den;

  den = (Ra*Ra + XL*XL);

  B = -XL/den;
  G = Ra/den;

  B *= p_mbase/p_sbase;
  G *= p_mbase/p_sbase;
  gridpack::ComplexType Y_a(G, B);
  return Y_a;
}


/**
 * Predict part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GridFormingGenerator::predictor_currentInjection(bool flag)
{
  p_INorton = INorton();
}

void gridpack::dynamic_simulation::GridFormingGenerator::computeModel(double t_inc, IntegrationStage int_flag, bool flag)
{
  S = V*conj(I);

  p_pg = real(S);
  p_qg = imag(S);
  Vt   = abs(V);

  if(!zero_Tf) {
    Pinv = P_filter_blk.getoutput(p_pg,t_inc,int_flag,true);
    Qinv = Q_filter_blk.getoutput(p_qg,t_inc,int_flag,true);
    Vmeas = V_filter_blk.getoutput(Vt, t_inc,int_flag,true);
  } else {
    Pinv = p_pg;
    Qinv = p_qg;
    Vmeas = Vt;
  }

  Edroop = Edroop_PI_blk.getoutput(Vset-mq*Qinv-Vmeas,t_inc,Edroop_min,Edroop_max,Edroop_min,Edroop_max,int_flag,true);
  
  Pmax_PI_blk_out = Pmax_PI_blk.getoutput(Pmax-Pinv,t_inc,int_flag,true);
  Pmin_PI_blk_out = Pmin_PI_blk.getoutput(Pmin-Pinv,t_inc,int_flag,true);

  double domega;

  domega = omega0*(mp*(Pset - Pinv) + Pmax_PI_blk_out + Pmin_PI_blk_out);

  omega = omega0 + domega;

  delta = Delta_blk.getoutput(domega,t_inc,int_flag,true);

  double Er,Ei;
  Er = Edroop*cos(delta);
  Ei = Edroop*sin(delta);

  E = gridpack::ComplexType(Er,Ei);
}
  
/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GridFormingGenerator::predictor(
    double t_inc, bool flag)
{
  computeModel(t_inc,PREDICTOR,flag);
}

/**
 * Corrector part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GridFormingGenerator::corrector_currentInjection(bool flag)
{
  p_INorton = INorton();
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GridFormingGenerator::corrector(
    double t_inc, bool flag)
{
  computeModel(t_inc,CORRECTOR,flag);
}

bool gridpack::dynamic_simulation::GridFormingGenerator::tripGenerator()
{
	p_tripped = true;
	
	return true;
}

/**
* return true if modify the generator parameters successfully
* input controlTyp: 0: GFI mp adjust; 1: GFI mq adjust; 2: GFI Pset adjust; 3: GFI Qset adjust; others: invalid 
* input newParValScaletoOrg:  GFI new parameter scale factor to the very initial parameter value at the begining of dynamic simulation
* 
*/
bool gridpack::dynamic_simulation::GridFormingGenerator::applyGeneratorParAdjustment(int controlType, double newParValScaletoOrg){
	
  if (controlType == 0){
    mp = mp_org * newParValScaletoOrg;
    return true;
  }else if(controlType == 1){
    mq = mq_org * newParValScaletoOrg;
    return true;
  }else if(controlType == 2) {
    Pset = Pset_org  * newParValScaletoOrg;
    return true;
  }else if(controlType == 3) {
    Vset = Vset_org  * newParValScaletoOrg;
    return true;
  }else{
    return false;
  }	
}

/**
 * Set voltage on each generator
 */
void gridpack::dynamic_simulation::GridFormingGenerator::setVoltage(
    gridpack::ComplexType voltage)
{
  V     = voltage;
  Vt    = abs(voltage);
  theta = atan2(imag(voltage), real(voltage));
  VR    = real(voltage);
  VI    = imag(voltage);
}

/**
 * Set frequency on each generator, frequency is perunit
 */
void gridpack::dynamic_simulation::GridFormingGenerator::setFreq(double dFreq)
{
   busfreq = dFreq;
}


/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
bool gridpack::dynamic_simulation::GridFormingGenerator::serialWrite(
    char* string, const int bufsize, const char *signal)
{
  double Pg,Qg;
  bool   ret=false;

  Pg = p_pg*p_mbase/p_sbase;
  Qg = p_qg*p_mbase/p_sbase;
  
  if (!strcmp(signal,"standard")) {
    sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f	%12.6f\n",
          p_bus_num, p_gen_id.c_str(), Edroop, delta, Pmax_PI_blk_out, Pmin_PI_blk_out, omega);
    ret = true;
  } else if (!strcmp(signal,"init_debug")) {
    sprintf(string," %8d  %2s Something\n",p_bus_num,p_gen_id.c_str());
    ret = true;
  } else if (!strcmp(signal,"debug_initial")) {
    sprintf(string,"Bus: %d PG: %f QG: %f\n",p_bus_num,Pg,Qg);
    ret = true;
  } else if (!strcmp(signal,"watch")) {
    if (getWatch()) {
      char buf[256];
      sprintf(string,",%12.6f,%12.6f, %12.6f, %12.6f,%12.6f",
	      Vt,Pg, Qg, busfreq,Im*p_mbase/p_sbase);
      ret = true;
    }
  } else if(!strcmp(signal,"watch_header")) {
    if(getWatch()) {
      char buf[128];
      std::string tag;
      if(p_gen_id[0] != ' ') {
	tag = p_gen_id;
      } else {
	tag = p_gen_id[1];
      }
      sprintf(buf,", %d_%s_V,%d_%s_Pg, %d_%s_Qg, %d_%s_freq, %d_%s_Igen",p_bus_num,tag.c_str(),p_bus_num,tag.c_str(),
	      p_bus_num,tag.c_str(),p_bus_num,tag.c_str(),	      p_bus_num,tag.c_str());
      
      if (strlen(buf) <= bufsize) {
	sprintf(string,"%s",buf);
	ret = true;
      } else {
      	ret = false;
      }
    }
  } else if (!strcmp(signal,"debug_initial")) {
    ret = false;
  }
  return ret;
}

/**
 * return a vector containing any generator values that are being
 * watched
 * @param vals vector of watched values
 */
void gridpack::dynamic_simulation::GridFormingGenerator::getWatchValues(
    std::vector<double> &vals)
{
  vals.clear();
  vals.push_back(delta);
  vals.push_back(omega);
  
  if (p_generatorObservationPowerSystemBase){
	vals.push_back(p_pg*p_mbase/p_sbase);  //output at system mva base
	vals.push_back(p_qg*p_mbase/p_sbase);  //output at system mva base
  }else{
	vals.push_back(p_pg);  //output at generator mva base
	vals.push_back(p_qg);  //output at generator mva base
  }
}
