/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   gensal.cpp
 * 
 * @brief: Salient pole generator model - no saturation  
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
#include "gensal.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::GensalGenerator::GensalGenerator(void)
{
  dx1d_0 = 0;
  dx2w_0 = 0;
  dx3Eqp_0 = 0;
  dx4Psidp_0 = 0;
  dx5Psiqpp_0 = 0;
  dx1d_1 = 0;
  dx2w_1 = 0;
  dx3Eqp_1 = 0;
  dx4Psidp_1 = 0;
  dx5Psiqpp_1 = 0;
  Vstab = 0.0;
  p_tripped = false;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::GensalGenerator::~GensalGenerator(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::GensalGenerator::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if(!data->getValue(CASE_SBASE,&p_sbase)) p_sbase = 100.0;

  data->getValue(BUS_NUMBER,&p_bus_id);
  data->getValue(GENERATOR_ID,&p_ckt,idx);
  if (!data->getValue(GENERATOR_PG, &p_pg,idx)) p_pg = 0.0;
  if (!data->getValue(GENERATOR_QG, &p_qg,idx)) p_qg = 0.0;
  if (!data->getValue(GENERATOR_STAT, &p_status,idx)) p_status = 0;

  p_pg *= p_sbase;
  p_qg *= p_sbase;

  if (!data->getValue(GENERATOR_MBASE, &MBase, idx)) MBase = 0.0; // MBase
  if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H, &H, idx)) H = 0.0; // H
  if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0, &D, idx)) D = 0.0; // D
  if (!data->getValue(GENERATOR_RESISTANCE, &Ra, idx)) Ra=0.0; // Ra
  if (!data->getValue(GENERATOR_XD, &Xd, idx)) Xd=0.0; // Xd
  if (!data->getValue(GENERATOR_XQ, &Xq, idx)) Xq=0.0; // Xq
  if (!data->getValue(GENERATOR_XDP, &Xdp, idx)) Xdp=0.0; // Xdp
  if (!data->getValue(GENERATOR_XDPP, &Xdpp, idx)) Xdpp=0.0; // Xdpp
  if (!data->getValue(GENERATOR_XL, &Xl, idx)) Xl=0.0; // Xl
  if (!data->getValue(GENERATOR_TDOP, &Tdop, idx)) Tdop=0.0; // Tdop
  if (!data->getValue(GENERATOR_TDOPP, &Tdopp, idx)) Tdopp=0.0; // Tdopp
  if (!data->getValue(GENERATOR_TQOPP, &Tqopp, idx)) Tqopp=0.0; // Tqopp
  if (!data->getValue(GENERATOR_S1, &S10, idx)) S10=0.17; // S10 TBD: check parser
  if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.55; // S12 TBD: check parser
  
  double tmp = sqrt(p_pg*p_pg +p_qg*p_qg);
  if ( tmp > MBase) {
    MBase = tmp*1.2;
  }
}

/**
 * Saturation function
 * @ param x
 */
double gridpack::dynamic_simulation::GensalGenerator::Sat(double x)
{
  double a_ = S12 / S10 - 1.0;
  double b_ = -2 * S12 / S10 + 2.4;
  double c_ = S12 / S10 - 1.44;
  double A = (-b_ - sqrt(b_ * b_ - 4 * a_ * c_)) / (2 * a_);
  double B = S10 / ((1.0 - A) * (1.0 - A));
  
  double tmp = x-A;

  if (tmp<0.0) {
    tmp = 0.0;
  }
  double result = B * tmp * tmp;
  
  return result; // Scaled Quadratic with 1.7.1 equations
}

/**
 * Initialize generator model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 */
void gridpack::dynamic_simulation::GensalGenerator::init(double mag,
    double ang, double ts)
{
  Vterm = mag;
  presentMag = mag;
  Theta = ang;
  presentAng = ang;

  double P = p_pg / MBase;
  double Q = p_qg / MBase;
  
  genP = P;
  genQ = Q;
  
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  Ir = (P * Vrterm + Q * Viterm) / (Vterm * Vterm);
  Ii = (P * Viterm - Q * Vrterm) / (Vterm * Vterm);

  x2w_0 = 0;
  x1d_0 = atan2(Viterm + Ir * Xq + Ii * Ra, Vrterm + Ir * Ra - Ii * Xq);
  Id = Ir * sin(x1d_0) - Ii * cos(x1d_0); // convert values to the dq axis
  Iq = Ir * cos(x1d_0) + Ii * sin(x1d_0); // convert values to the dq axis
  double Vqterm = Vrterm * cos(x1d_0) + Viterm * sin(x1d_0); // convert values to the dq axis
  x5Psiqpp_0 = (Xdpp - Xq) * Iq;

  double Psiq = x5Psiqpp_0 - Iq * Xdpp;
  double Psid = Vqterm + Ra * Iq;
  double Psidpp = Psid + Id * Xdpp;
  x4Psidp_0 = Psidpp - Id * (Xdpp - Xl);
  x3Eqp_0 = x4Psidp_0 + Id * (Xdp - Xl);
  Efd = x3Eqp_0 * (1 + Sat(x3Eqp_0)) + Id * (Xd - Xdp);

  LadIfd = Efd;
  Pmech = Psid * Iq - Psiq * Id;

  Efdinit = Efd;
  Pmechinit = Pmech;

  if (p_hasExciter){
    p_exciter = getExciter();
    p_exciter->setVterminal(Vterm); 
    p_exciter->setVcomp(mag); 
    p_exciter->setFieldVoltage(Efd);
    p_exciter->setFieldCurrent(LadIfd);
    //---yuan add below 20231024---//
    p_exciter->setIri(Ir, Ii); 
    //---yuan add above 20231024---//
    p_exciter->setExtBusNum(p_bus_id);
    p_exciter->init(mag, ang, ts);
  }

  if (p_hasGovernor) {
    p_governor = getGovernor();
    p_governor->setMechanicalPower(Pmech);
    p_governor->setRotorSpeedDeviation(x2w_0); // set Speed Deviation w for wsieg1 
    p_governor->init(mag, ang, ts);
  }
  
  if (p_hasPss) {
    p_pss = getPss();
    p_pss->setOmega(x2w_0);
    
    p_pss->init(mag, ang, ts);
  }
  
  p_Norton_Ya = NortonImpedence();
}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType gridpack::dynamic_simulation::GensalGenerator::INorton()
{
  return p_INorton;
}


/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType gridpack::dynamic_simulation::GensalGenerator::NortonImpedence()
{
  double ra = Ra * p_sbase / MBase;
  double xd = Xdpp * p_sbase / MBase;
  B = -xd / (ra * ra + xd * xd);
  G = ra / (ra * ra + xd * xd);

  gridpack::ComplexType Y_a(G, B);
  return Y_a;
}

/**
 * Predict part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GensalGenerator::predictor_currentInjection(bool flag)
{
  if (!flag) {
    x1d_0 = x1d_1;
    x2w_0 = x2w_1;
    x3Eqp_0 = x3Eqp_1;
    x4Psidp_0 = x4Psidp_1;
    x5Psiqpp_0 = x5Psiqpp_1; 
  }  
  // Calculate INorton_full
  // Admittance
  B = -Xdpp / (Ra * Ra + Xdpp * Xdpp);
  G = Ra / (Ra * Ra + Xdpp * Xdpp);

  // Setup
  double Psiqpp = x5Psiqpp_0; // this will be different for GENROU
  double Psidpp = + x3Eqp_0 * (Xdpp - Xl) / (Xdp - Xl)
    + x4Psidp_0* (Xdp - Xdpp) / (Xdp - Xl);
  double Vd = -Psiqpp;// * (1 + x2w_0);
  double Vq = +Psidpp;// * (1 + x2w_0);
  Vterm = presentMag;
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  
  double Vdterm = Vrterm * sin(x1d_0) - Viterm * cos(x1d_0);
  double Vqterm = Vrterm * cos(x1d_0) + Viterm * sin(x1d_0);
  
  gridpack::ComplexType vt_complex_tmp = gridpack::ComplexType(Vrterm, Viterm); 

  //DQ Axis
  Id = (Vd - Vdterm) * G - (Vq - Vqterm) * B;
  Iq = (Vd - Vdterm) * B + (Vq - Vqterm) * G;
  double Idnorton = Vd * G - Vq * B;
  double Iqnorton = Vd * B + Vq * G;
  //Network
  Ir = + Id * sin(x1d_0) + Iq * cos(x1d_0);
  Ii = - Id * cos(x1d_0) + Iq * sin(x1d_0);
  
  genP = Vrterm*Ir + Viterm*Ii;
  genQ = Viterm*Ir - Vrterm*Ii;
  
  IrNorton = + Idnorton * sin(x1d_0) + Iqnorton * cos(x1d_0);
  IiNorton = - Idnorton * cos(x1d_0) + Iqnorton * sin(x1d_0);
  
  IrNorton = IrNorton * MBase / p_sbase; 
  IiNorton = IiNorton * MBase / p_sbase; 
  
  if (getGenStatus()){
    if (p_tripped){
      p_INorton = p_Norton_Ya*vt_complex_tmp;
    }else{
      p_INorton = gridpack::ComplexType(IrNorton, IiNorton);	
    }		
  } else {
    p_INorton = gridpack::ComplexType(0.0, 0.0);
  }
} 

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GensalGenerator::predictor(
    double t_inc, bool flag)
{
  if (getGenStatus() == false) {
    x1d_0 = 0.0;
    x2w_0 = -1.0;
    x3Eqp_0 = 0.0;
    x4Psidp_0 = 0.0;
    x5Psiqpp_0 = 0.0;
    x1d_1 = 0.0;
    x2w_1 = -1.0;
    x3Eqp_1 = 0.0;
    x4Psidp_1 = 0.0;
    x5Psiqpp_1 = 0.0;
    genP = 0.0;
    genQ = 0.0;
    return;
  }
	    
  if (p_hasExciter){	  
    p_exciter = getExciter();
    Efd = p_exciter->getFieldVoltage();
  } else {
    Efd = Efdinit;
  }

  if (p_hasGovernor){
    p_governor = getGovernor();
    Pmech = p_governor->getMechanicalPower();
  } else {
    Pmech = Pmechinit;
  }
  
  if (!flag) {
    x1d_0 = x1d_1;
    x2w_0 = x2w_1;
    x3Eqp_0 = x3Eqp_1;
    x4Psidp_0 = x4Psidp_1;
    x5Psiqpp_0 = x5Psiqpp_1; 
  }  
  
  B = -Xdpp / (Ra * Ra + Xdpp * Xdpp);
  G = Ra / (Ra * Ra + Xdpp * Xdpp);

  // Setup
  double Psiqpp = x5Psiqpp_0; // this will be different for GENROU
  double Psidpp = + x3Eqp_0 * (Xdpp - Xl) / (Xdp - Xl)
    + x4Psidp_0* (Xdp - Xdpp) / (Xdp - Xl);
  double Vd = -Psiqpp;// * (1 + x2w_0);
  double Vq = +Psidpp;// * (1 + x2w_0);
  Vterm = presentMag;
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  
  double Vdterm = Vrterm * sin(x1d_0) - Viterm * cos(x1d_0);
  double Vqterm = Vrterm * cos(x1d_0) + Viterm * sin(x1d_0);
  
  gridpack::ComplexType vt_complex_tmp = gridpack::ComplexType(Vrterm, Viterm); 
  
  //DQ Axis
  Id = (Vd - Vdterm) * G - (Vq - Vqterm) * B;
  Iq = (Vd - Vdterm) * B + (Vq - Vqterm) * G;

  double pi = 4.0*atan(1.0);
  double Psiq = x5Psiqpp_0 - Iq * Xdpp;
  Psidpp = x3Eqp_0 * (Xdpp - Xl) / (Xdp - Xl) 
                + x4Psidp_0 * (Xdp - Xdpp) / (Xdp - Xl);
  double Psid = Psidpp - Id * Xdpp;
  double Telec = Psid * Iq - Psiq * Id;
  double TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl))
               * ((-x4Psidp_0 - (Xdp - Xl) * Id + x3Eqp_0));
  LadIfd = x3Eqp_0 * (1 + Sat(x3Eqp_0)) + (Xd - Xdp) * (Id + TempD); // update Ifd later
  dx1d_0 = x2w_0 * 2 * pi * 60; // 60 represents the nominal frequency of 60 Hz
  dx2w_0 = 1 / (2 * H) * ((Pmech - D * x2w_0) / (1 + x2w_0) - Telec);
  dx3Eqp_0 = (Efd - LadIfd) / Tdop;
  dx4Psidp_0 = (-x4Psidp_0 - (Xdp - Xl) * Id + x3Eqp_0) / Tdopp;
  dx5Psiqpp_0 = (-x5Psiqpp_0 - (Xq - Xdpp) * Iq) / Tqopp;

  x1d_1 = x1d_0 + dx1d_0 * t_inc;
  x2w_1 = x2w_0 + dx2w_0 * t_inc;
  x3Eqp_1 = x3Eqp_0 + dx3Eqp_0 * t_inc;
  x4Psidp_1 = x4Psidp_0 + dx4Psidp_0 * t_inc;
  x5Psiqpp_1 = x5Psiqpp_0 + dx5Psiqpp_0 * t_inc;

  if (p_hasPss) {
    p_pss = getPss();
    p_pss->setOmega(x2w_1);
    p_pss->predictor(t_inc, flag);
    Vstab = p_pss->getVstab();
  } else {
    Vstab = 0.0;
  }
  
  if (p_hasExciter){
    if (p_hasPss) { 
      p_exciter->setVstab(Vstab);
    }
    p_exciter->setVterminal(presentMag);
    p_exciter->setVcomp(presentMag); //TBD update to Vcomp 
    //  p_exciter->setFieldVoltage(Efd);
    p_exciter->setFieldCurrent(LadIfd);
    
    p_exciter->predictor(t_inc, flag);
  }

  if (p_hasGovernor){
    p_governor->setRotorSpeedDeviation(x2w_0);
    p_governor->predictor(t_inc, flag);
  }

  if (p_tripped){
    x1d_0 = 0.0;
    x2w_0 = -1.0;
    x3Eqp_0 = 0.0;
    x4Psidp_0 = 0.0;
    x5Psiqpp_0 = 0.0;
    x1d_1 = 0.0;
    x2w_1 = -1.0;
    x3Eqp_1 = 0.0;
    x4Psidp_1 = 0.0;
    x5Psiqpp_1 = 0.0;	
    genP = 0.0;
    genQ = 0.0;
  }
}

/**
 * Corrector part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GensalGenerator::corrector_currentInjection(bool flag)
{
  // Calculate INorton_full
  // Admittance
  B = -Xdpp / (Ra * Ra + Xdpp * Xdpp);
  G = Ra / (Ra * Ra + Xdpp * Xdpp);

  // Setup
  double Psiqpp = x5Psiqpp_1; // this will be different for GENROU
  double Psidpp = + x3Eqp_1 * (Xdpp - Xl) / (Xdp - Xl)
    + x4Psidp_1 * (Xdp - Xdpp) / (Xdp - Xl);
  double Vd = -Psiqpp;// * (1 + x2w_1);
  double Vq = +Psidpp;// * (1 + x2w_1);
  Vterm = presentMag;
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  double Vdterm = Vrterm * sin(x1d_1) - Viterm * cos(x1d_1);
  double Vqterm = Vrterm * cos(x1d_1) + Viterm * sin(x1d_1);
  
  gridpack::ComplexType vt_complex_tmp = gridpack::ComplexType(Vrterm, Viterm); 
	
  //DQ Axis
  Id = (Vd - Vdterm) * G - (Vq - Vqterm) * B;
  Iq = (Vd - Vdterm) * B + (Vq - Vqterm) * G;
  double Idnorton = Vd * G - Vq * B;
  double Iqnorton = Vd * B + Vq * G;
  //Network
  Ir = + Id * sin(x1d_1) + Iq * cos(x1d_1);
  Ii = - Id * cos(x1d_1) + Iq * sin(x1d_1);
  
  //genP and genQ should not be updated here to output values, as here the states values are from the predictor, is not accurate enough
  double genP1 = Vrterm*Ir + Viterm*Ii;
  double genQ1 = Viterm*Ir - Vrterm*Ii;
  
  IrNorton = + Idnorton * sin(x1d_1) + Iqnorton * cos(x1d_1);
  IiNorton = - Idnorton * cos(x1d_1) + Iqnorton * sin(x1d_1); 
  
  
  IrNorton = IrNorton * MBase / p_sbase; 
  IiNorton = IiNorton * MBase / p_sbase; 
  
  if (getGenStatus()){	  
    if (p_tripped){
      p_INorton = p_Norton_Ya*vt_complex_tmp;
    }else{
      p_INorton = gridpack::ComplexType(IrNorton, IiNorton);		
    }	 	  
  } else {
    p_INorton = gridpack::ComplexType(0.0, 0.0);
  }
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GensalGenerator::corrector(
    double t_inc, bool flag)
{
  if (getGenStatus() == false) {
    x1d_0 = 0.0;
    x2w_0 = -1.0;
    x3Eqp_0 = 0.0;
    x4Psidp_0 = 0.0;
    x5Psiqpp_0 = 0.0;
    x1d_1 = 0.0;
    x2w_1 = -1.0;
    x3Eqp_1 = 0.0;
    x4Psidp_1 = 0.0;
    x5Psiqpp_1 = 0.0;
    genP = 0.0;
    genQ = 0.0;
    return;
  }

  if (p_hasExciter){	  
    p_exciter = getExciter();
    Efd = p_exciter->getFieldVoltage();
  } else {
    Efd = Efdinit;
  }

  if (p_hasGovernor){
    p_governor = getGovernor();
    Pmech = p_governor->getMechanicalPower();
  } else {
    Pmech = Pmechinit;
  } 

  B = -Xdpp / (Ra * Ra + Xdpp * Xdpp);
  G = Ra / (Ra * Ra + Xdpp * Xdpp);

  // Setup
  double Psiqpp = x5Psiqpp_1; // this will be different for GENROU
  double Psidpp = + x3Eqp_1 * (Xdpp - Xl) / (Xdp - Xl)
    + x4Psidp_1 * (Xdp - Xdpp) / (Xdp - Xl);
  double Vd = -Psiqpp;// * (1 + x2w_1);
  double Vq = +Psidpp;// * (1 + x2w_1);
  Vterm = presentMag;
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  double Vdterm = Vrterm * sin(x1d_1) - Viterm * cos(x1d_1);
  double Vqterm = Vrterm * cos(x1d_1) + Viterm * sin(x1d_1);
  
  gridpack::ComplexType vt_complex_tmp = gridpack::ComplexType(Vrterm, Viterm); 
	
  //DQ Axis
  Id = (Vd - Vdterm) * G - (Vq - Vqterm) * B;
  Iq = (Vd - Vdterm) * B + (Vq - Vqterm) * G;

  double pi = 4.0*atan(1.0);
  double Psiq = x5Psiqpp_1 - Iq * Xdpp;
  Psidpp = x3Eqp_1 * (Xdpp - Xl) / (Xdp - Xl) 
    + x4Psidp_1 * (Xdp - Xdpp) / (Xdp - Xl);
  double Psid = Psidpp - Id * Xdpp;
  double Telec = Psid * Iq - Psiq * Id;
  double TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl))
    * ((-x4Psidp_1 - (Xdp - Xl) * Id + x3Eqp_1));
  LadIfd = x3Eqp_1 * (1 + Sat(x3Eqp_1)) + (Xd - Xdp) * (Id + TempD);
  dx1d_1 = x2w_1 * 2 * pi * 60; // 60 represents the nominal frequency of 60 Hz
  dx2w_1 = 1 / (2 * H) * ((Pmech - D * x2w_1) / (1 + x2w_1) - Telec);
  dx3Eqp_1 = (Efd - LadIfd) / Tdop;
  dx4Psidp_1 = (-x4Psidp_1 - (Xdp - Xl) * Id + x3Eqp_1) / Tdopp;
  dx5Psiqpp_1 = (-x5Psiqpp_1 - (Xq - Xdpp) * Iq) / Tqopp;
  
  x1d_1 = x1d_0 + (dx1d_0 + dx1d_1) / 2.0 * t_inc;
  x2w_1 = x2w_0 + (dx2w_0 + dx2w_1) / 2.0 * t_inc;
  x3Eqp_1 = x3Eqp_0 + (dx3Eqp_0 + dx3Eqp_1) / 2.0 * t_inc;
  x4Psidp_1 = x4Psidp_0 + (dx4Psidp_0 + dx4Psidp_1) / 2.0 * t_inc;
  x5Psiqpp_1 = x5Psiqpp_0 + (dx5Psiqpp_0 + dx5Psiqpp_1) / 2.0 * t_inc;
  
  if (p_hasPss) {
    p_pss = getPss();
    p_pss->setOmega(x2w_1);
    p_pss->corrector(t_inc, flag);
    Vstab = p_pss->getVstab();
  } else {
    Vstab = 0.0;
  }	  
  
  if (p_hasExciter){
    if (p_hasPss) { 
      p_exciter->setVstab(Vstab);
    }
    
    p_exciter->setVterminal(presentMag);
    p_exciter->setVcomp(presentMag); 
    //  p_exciter->setFieldVoltage(Efd);
    p_exciter->setFieldCurrent(LadIfd);
    p_exciter->corrector(t_inc, flag);
  }
  
  if (p_hasGovernor){
    p_governor->setRotorSpeedDeviation(x2w_1); //note previous version here is x2w_0, not correct
    p_governor->corrector(t_inc, flag);
  }
  
  if (p_tripped){
    x1d_0 = 0.0;
    x2w_0 = -1.0;
    x3Eqp_0 = 0.0;
    x4Psidp_0 = 0.0;
    x5Psiqpp_0 = 0.0;
    x1d_1 = 0.0;
    x2w_1 = -1.0;
    x3Eqp_1 = 0.0;
    x4Psidp_1 = 0.0;
    x5Psiqpp_1 = 0.0;
    genP = 0.0;
    genQ = 0.0;		
  } 
}

void gridpack::dynamic_simulation::GensalGenerator::setWideAreaFreqforPSS(double freq)
{
  p_wideareafreq = freq;
  if (p_hasPss){
    p_pss = getPss();
    p_pss->setWideAreaFreqforPSS(freq);	
  }
}

bool gridpack::dynamic_simulation::GensalGenerator::tripGenerator()
{
  p_tripped = true;
  
  return true;
}

/**
 * Set voltage on each generator
 */
void gridpack::dynamic_simulation::GensalGenerator::setVoltage(gridpack::ComplexType voltage)
{
  presentMag = abs(voltage);
  presentAng = atan2(imag(voltage), real(voltage));  
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double gridpack::dynamic_simulation::GensalGenerator::getFieldVoltage()
{
  return Efd;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
bool gridpack::dynamic_simulation::GensalGenerator::serialWrite(
    char* string, const int bufsize, const char *signal)
{
  bool ret = false;
  if (!strcmp(signal,"standard")) {
    sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f	%12.6f\n",
	    p_bus_id, p_ckt.c_str(), x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqpp_1);
    ret = true;
  } else if (!strcmp(signal,"init_debug")) {
    sprintf(string," %8d  %2s Something\n",p_bus_id,p_ckt.c_str());
    ret = true;
  } else if (!strcmp(signal,"debug_initial")) {
    sprintf(string,"Bus: %d PG: %f QG: %f\n",p_bus_id,p_pg,p_qg);
    ret = true;
  } else if(!strcmp(signal,"watch_header")) {
    if(getWatch()) {
      char buf[128];
      std::string tag;
      if(p_ckt[0] != ' ') {
	tag = p_ckt;
      } else {
	tag = p_ckt[1];
      }
      sprintf(buf,", %d_%s_V, %d_%s_Pg, %d_%s_Qg,%d_%s_angle, %d_%s_speed, %d_%s_Efd, %d_%s_Pm",p_bus_id,tag.c_str(),p_bus_id,tag.c_str(),p_bus_id,tag.c_str(),p_bus_id,tag.c_str(),
	      p_bus_id,tag.c_str(),p_bus_id,tag.c_str(),p_bus_id,tag.c_str());
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
      char buf[256];
      sprintf(buf,",%f,%f,%f,%f, %f, %f, %f",
	      Vterm,genP*MBase/p_sbase,genQ*MBase/p_sbase,x1d_1, x2w_1+1.0, Efd,Pmech);
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
    ret = false;
  }
  return ret;
}

/**
 * return a vector containing any generator values that are being
 * watched
 * @param vals vector of watched values
 */
void gridpack::dynamic_simulation::GensalGenerator::getWatchValues(
    std::vector<double> &vals)
{
  vals.clear();
  vals.push_back(x1d_1);
  vals.push_back(x2w_1+1.0);
  
  if (p_generatorObservationPowerSystemBase){
	vals.push_back(genP*MBase/p_sbase);  //output at system mva base
	vals.push_back(genQ*MBase/p_sbase);  //output at system mva base
  }else{
	vals.push_back(genP);  //output at generator mva base
	vals.push_back(genQ);  //output at generator mva base
  }
}

/**
 * Set internal state parameter in generator
 * @param name character string corresponding to state variable
 * @param value new value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::GensalGenerator::setState(std::string name, 
    double value)
{ 
  return false;
}

/**
 * Get internal state parameter in generator
 * @param name character string corresponding to state variable
 * @param value current value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::GensalGenerator::getState(std::string name, 
    double *value)
{ 
  return false;
}
