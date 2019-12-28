/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   gensal.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 24, 2015
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
#include "gensal.hpp"
//#include "exdc1.hpp"

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
  p_sbase = 100.0;

  /*int gstatus;
  double pg, qg, mva, r, dstr, dtr;
  double h, d0;*/

  data->getValue(BUS_NUMBER,&p_bus_id);
  data->getValue(GENERATOR_ID,&p_ckt,idx);
  if (!data->getValue(GENERATOR_PG, &p_pg,idx)) p_pg = 0.0;
  if (!data->getValue(GENERATOR_QG, &p_qg,idx)) p_qg = 0.0;
  if (!data->getValue(GENERATOR_STAT, &p_status,idx)) p_status = 0;
  //printf("load p_pg = %f, p_qg = %f\n", p_pg, p_qg);
  p_pg *= p_sbase;
  p_qg *= p_sbase;

  if (!data->getValue(GENERATOR_MBASE, &MVABase, idx)) MVABase = 0.0; // MVABase
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
  //if (!data->getValue(GENERATOR_S1, &S10, idx)) S10=0.0; // S10
  //if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.0; // S12
  if (!data->getValue(GENERATOR_S1, &S10, idx)) S10=0.17; // S10 TBD: check parser
  if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.55; // S12 TBD: check parser
  //printf("load S10 = %f, S12 = %f\n", S10, S12);
  //if (!data->getValue(GENERATOR_XQP, &Xqp, idx)) Xqp=0.0; // Xqp
  
  // printf("gensal parameters: %12.6f, %12.6f,%12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f  \n", H, D, Ra, Xd, Xq, Xdp, Xdpp, Xl, Tdop, Tdopp, Tqopp, S10, S12);
}

/**
 * Saturation function
 * @ param x
 */
double gridpack::dynamic_simulation::GensalGenerator::Sat(double x)
{
    double a_ = S12 / S10 - 1.0 / 1.2;
    double b_ = -2 * S12 / S10 + 2;
    double c_ = S12 / S10 - 1.2;
    double A = (-b_ - sqrt(b_ * b_ - 4 * a_ * c_)) / (2 * a_);
    double B = S10 / ((1.0 - A) * (1.0 - A));
    double result = B * (x - A) * (x - A) / x;
    //printf("a = %f, b = %f, c = %f, A = %f, B = %f, S12 = %f, S10 = %f\n", a_, b_, c_, A, B, S12, S10);
    //printf("Sat result = %f\n", result); 
	
	// the following is another method for saturation computation, add by renke
	/*
	double a_ = S12 / S10 - 1.0;
    double b_ = -2 * S12 / S10 + 2.4;
    double c_ = S12 / S10 - 1.44;
    double A = (-b_ - sqrt(b_ * b_ - 4 * a_ * c_)) / (2 * a_);
    double B = S10 / ((1.0 - A) * (1.0 - A));
	
	double tmp = x-A;
	//double tmpin = tmp;
	if (tmp<0.0) {
		tmp = 0.0;
	}
    double result = B * tmp * tmp;
	*/
	
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
  //printf("Step0 gen%d mag = %f\n", p_bus_id, mag);
  Vterm = mag;
  presentMag = mag;
  Theta = ang;
  presentAng = ang;
  double P = p_pg / MVABase;
  double Q = p_qg / MVABase;
  //printf("p_pg = %f, p_qg = %f, MVABase = %f\n", p_pg, p_qg, MVABase);
  //printf("Vterm = %f, Theta = %f, P = %f, Q = %f\n", Vterm, Theta, P, Q);
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  Ir = (P * Vrterm + Q * Viterm) / (Vterm * Vterm);
  Ii = (P * Viterm - Q * Vrterm) / (Vterm * Vterm);
  //printf("Ir = %f, Ii = %f\n", Ir, Ii);
  x2w_0 = 0;
  x1d_0 = atan2(Viterm + Ir * Xq + Ii * Ra, Vrterm + Ir * Ra - Ii * Xq);
  Id = Ir * sin(x1d_0) - Ii * cos(x1d_0); // convert values to the dq axis
  Iq = Ir * cos(x1d_0) + Ii * sin(x1d_0); // convert values to the dq axis
  double Vqterm = Vrterm * cos(x1d_0) + Viterm * sin(x1d_0); // convert values to the dq axis
  x5Psiqpp_0 = (Xdpp - Xq) * Iq;
  //printf("Xdpp = %f, Xq = %f, Iq = %f\n", Xdpp, Xq, Iq);
  double Psiq = x5Psiqpp_0 - Iq * Xdpp;
  double Psid = Vqterm + Ra * Iq;
  double Psidpp = Psid + Id * Xdpp;
  x4Psidp_0 = Psidpp - Id * (Xdpp - Xl);
  x3Eqp_0 = x4Psidp_0 + Id * (Xdp - Xl);
  Efd = x3Eqp_0 * (1 + Sat(x3Eqp_0)) + Id * (Xd - Xdp);
  //printf("x3Eqp_0 = %f, Id = %f, Xd = %f, Xdp = %f\n", x3Eqp_0, Id, Xd, Xdp);
  //printf("gensal init Ra=%f \n", Ra);
  LadIfd = Efd;
  Pmech = Psid * Iq - Psiq * Id;

  //printf("gensal init: %f\t%f\t%f\t%f\t%f\n", x1d_0, x2w_0, x3Eqp_0, x4Psidp_0, x5Psiqpp_0);
  //printf("gensal init: Efd = %f, Pmech = %f\n", Efd, Pmech);
  Efdinit = Efd;
  Pmechinit = Pmech;

  if (p_hasExciter){
	p_exciter = getExciter();
	p_exciter->setVterminal(Vterm); 
//  printf ("gensal init Vterm = %f \n", Vterm);
//  printf ("gensal init setVcomp abs(Vterm) = %f \n", abs(Vterm));
//  printf ("gensal init mag = %f \n", mag);
//  p_exciter->setVcomp(abs(Vterm)); //TBD Need to updated later to calculate Vcomp 
	p_exciter->setVcomp(mag); 
//  printf("esst1a Vcomp= %f\n", gridpack::dynamic_simulation::Esst1aModel::p_exciter->Vcomp);
	p_exciter->setFieldVoltage(Efd);
	p_exciter->setFieldCurrent(LadIfd);
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

  // Initialize other variables 
  /*p_mac_ang_s1 = gridpack::ComplexType(0.0,0.0);
  p_mac_spd_s1 = gridpack::ComplexType(0.0,0.0);
  p_dmac_ang_s0 = gridpack::ComplexType(0.0,0.0);
  p_dmac_spd_s0 = gridpack::ComplexType(0.0,0.0);
  p_dmac_ang_s1 = gridpack::ComplexType(0.0,0.0);
  p_dmac_spd_s1 = gridpack::ComplexType(0.0,0.0);
  p_eprime_s1 = gridpack::ComplexType(0.0,0.0);
  p_INorton = gridpack::ComplexType(0.0,0.0);*/
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
  /*double ra = Ra * p_sbase / MVABase;
  double xd = Xdpp * p_sbase / MVABase;
  B = -xd / (ra + xd);
  G = ra / (ra + xd);
  gridpack::ComplexType Y_a(B, G);
  return Y_a;*/
  double ra = Ra * p_sbase / MVABase;
  double xd = Xdpp * p_sbase / MVABase;
  //printf("Ra = %f, Xdpp = %f\n", Ra, Xdpp);
  //printf("ra = %f, xd = %f, p_sbase = %f, MVABase = %f\n", ra, xd, p_sbase, MVABase);
  B = -xd / (ra * ra + xd * xd);
  G = ra / (ra * ra + xd * xd);
  //printf("B = %f, G = %f\n", B, G);
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
  //printf("Xdpp = %f, Ra = %f, B = %f, G = %f\n", Xdpp, Ra, B, G);
  // Setup
  double Psiqpp = x5Psiqpp_0; // this will be different for GENROU
  double Psidpp = + x3Eqp_0 * (Xdpp - Xl) / (Xdp - Xl) 
                + x4Psidp_0* (Xdp - Xdpp) / (Xdp - Xl);
  double Vd = -Psiqpp * (1 + x2w_0);
  double Vq = +Psidpp * (1 + x2w_0);
  Vterm = presentMag;
  //printf("Gensal predictor_currentInjection: %d %f\n", p_bus_id, Vterm);
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  double Vdterm = Vrterm * sin(x1d_0) - Viterm * cos(x1d_0);
  double Vqterm = Vrterm * cos(x1d_0) + Viterm * sin(x1d_0);
  //printf("x5Psiqpp_0 = %f, x3Eqp_0 = %f, Xl = %f, Xdp = %f, Psidpp = %f\n", x5Psiqpp_0, x3Eqp_0, Xl, Xdp, Psidpp);
  //printf("x2w_0 = %f, Vd = %f, Vq = %f, Vrterm = %f, Viterm = %f, Vdterm = %f, Vqterm = %f, Theta = %f\n", x2w_0, Vd, Vq, Vrterm, Viterm, Vdterm, Vqterm, Theta);
  //DQ Axis
  Id = (Vd - Vdterm) * G - (Vq - Vqterm) * B;
  Iq = (Vd - Vdterm) * B + (Vq - Vqterm) * G;
  double Idnorton = Vd * G - Vq * B;
  double Iqnorton = Vd * B + Vq * G;
  //Network
  Ir = + Id * sin(x1d_0) + Iq * cos(x1d_0);
  Ii = - Id * cos(x1d_0) + Iq * sin(x1d_0);
  IrNorton = + Idnorton * sin(x1d_0) + Iqnorton * cos(x1d_0);
  IiNorton = - Idnorton * cos(x1d_0) + Iqnorton * sin(x1d_0);
  IrNorton = IrNorton * MVABase / p_sbase; 
  IiNorton = IiNorton * MVABase / p_sbase; 
  //gridpack::ComplexType INorton(IrNorton, IiNorton);
  if (getGenStatus()){
	  p_INorton = gridpack::ComplexType(IrNorton, IiNorton);	  
  }else {
	  p_INorton = gridpack::ComplexType(0.0, 0.0);
  }
  //printf("gensal::predictor_currentInjection: presentMag = %f, presentAng = %f \n", presentMag, presentAng);
  //printf("gensal::predictor_currentInjuction: p_INorton = %f, %f \n", real(p_INorton), imag(p_INorton));
} 

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GensalGenerator::predictor(
    double t_inc, bool flag)
{
  ///printf("\n***** GEN %d Predicator:\n", p_bus_id);
  if (getGenStatus()){
	    
  if (p_hasExciter){	  
		p_exciter = getExciter();
		Efd = p_exciter->getFieldVoltage();
		// if (p_bus_id == 8022) Efd = Efdinit;  //use this one when relay is added
  }else{
		Efd = Efdinit;
  }

  if (p_hasGovernor){
	p_governor = getGovernor();
	Pmech = p_governor->getMechanicalPower();
	// if (p_bus_id == 8022) Pmech = Pmechinit;
  }else{
	Pmech = Pmechinit;
  }
 
//  printf("predictor: Efd = %f, Pmech = %f\n", Efd, Pmech); 

  if (!flag) {
    x1d_0 = x1d_1;
    x2w_0 = x2w_1;
    x3Eqp_0 = x3Eqp_1;
    x4Psidp_0 = x4Psidp_1;
    x5Psiqpp_0 = x5Psiqpp_1; 
  }  
    
  double pi = 4.0*atan(1.0);
  double Psiq = x5Psiqpp_0 - Iq * Xdpp;
  double Psidpp = x3Eqp_0 * (Xdpp - Xl) / (Xdp - Xl) 
                + x4Psidp_0 * (Xdp - Xdpp) / (Xdp - Xl);
  double Psid = Psidpp - Id * Xdpp;
  double Telec = Psid * Iq - Psiq * Id;
  double TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl))
               * ((-x4Psidp_0 - (Xdp - Xl) * Id + x3Eqp_0));
  LadIfd = x3Eqp_0 * (1 + Sat(x3Eqp_0)) + (Xd - Xdp) * (Id + TempD); // update Ifd later
  dx1d_0 = x2w_0 * 2 * pi * 60; // 60 represents the nominal frequency of 60 Hz
  //printf("H = %f, Pmech = %f, D = %f, x2w_0 = %f, Telec = %f\n", H, Pmech, D, x2w_0, Telec);
  dx2w_0 = 1 / (2 * H) * ((Pmech - D * x2w_0) / (1 + x2w_0) - Telec); //TBD: call Governor for Pmech (Done)
  dx3Eqp_0 = (Efd - LadIfd) / Tdop; //TBD: call Exciter for Efd and LadIfd (Done)
  dx4Psidp_0 = (-x4Psidp_0 - (Xdp - Xl) * Id + x3Eqp_0) / Tdopp;
  dx5Psiqpp_0 = (-x5Psiqpp_0 - (Xq - Xdpp) * Iq) / Tqopp;

  x1d_1 = x1d_0 + dx1d_0 * t_inc;
  x2w_1 = x2w_0 + dx2w_0 * t_inc;
  x3Eqp_1 = x3Eqp_0 + dx3Eqp_0 * t_inc;
  x4Psidp_1 = x4Psidp_0 + dx4Psidp_0 * t_inc;
  x5Psiqpp_1 = x5Psiqpp_0 + dx5Psiqpp_0 * t_inc;
  //printf("gensal pbusid%d, dx: %f\t%f\t%f\t%f\t%f\n", p_bus_id, dx1d_0, dx2w_0, dx3Eqp_0, dx4Psidp_0, dx5Psiqpp_0);
  ///printf("gensal x: %f\t%f\t%f\t%f\t%f\n", x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqpp_1);
  //printf ("gensal predictor presentMag=%f\n", presentMag);
  if (p_hasPss) {
	p_pss = getPss();
	p_pss->setOmega(x2w_1);
	p_pss->predictor(t_inc, flag);
	Vstab = p_pss->getVstab();
	}else{
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
//  printf("predictor gensal: Efd = %f, Pmech = %f\n", Efd, Pmech); 
  }else {
	x1d_0 = 0.0;
    x2w_0 = 0.0;
    x3Eqp_0 = 0.0;
    x4Psidp_0 = 0.0;
    x5Psiqpp_0 = 0.0;
	x1d_1 = 0.0;
    x2w_1 = 0.0;
    x3Eqp_1 = 0.0;
    x4Psidp_1 = 0.0;
    x5Psiqpp_1 = 0.0;
  
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
  //printf("B = %f, G = %f\n", B, G);
  // Setup
  double Psiqpp = x5Psiqpp_1; // this will be different for GENROU
  double Psidpp = + x3Eqp_1 * (Xdpp - Xl) / (Xdp - Xl) 
                + x4Psidp_1 * (Xdp - Xdpp) / (Xdp - Xl);
  double Vd = -Psiqpp * (1 + x2w_1);
  double Vq = +Psidpp * (1 + x2w_1);
  Vterm = presentMag;
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  double Vdterm = Vrterm * sin(x1d_1) - Viterm * cos(x1d_1);
  double Vqterm = Vrterm * cos(x1d_1) + Viterm * sin(x1d_1);
  //DQ Axis
  Id = (Vd - Vdterm) * G - (Vq - Vqterm) * B;
  Iq = (Vd - Vdterm) * B + (Vq - Vqterm) * G;
  double Idnorton = Vd * G - Vq * B;
  double Iqnorton = Vd * B + Vq * G;
  //Network
  Ir = + Id * sin(x1d_1) + Iq * cos(x1d_1);
  Ii = - Id * cos(x1d_1) + Iq * sin(x1d_1);
  IrNorton = + Idnorton * sin(x1d_1) + Iqnorton * cos(x1d_1);
  IiNorton = - Idnorton * cos(x1d_1) + Iqnorton * sin(x1d_1); 
  IrNorton = IrNorton * MVABase / p_sbase; 
  IiNorton = IiNorton * MVABase / p_sbase; 
  //gridpack::ComplexType INorton(IrNorton, IiNorton);
  if (getGenStatus()){
	  p_INorton = gridpack::ComplexType(IrNorton, IiNorton);	 	  
  }else {
	  p_INorton = gridpack::ComplexType(0.0, 0.0);
  }
  //printf("gensal::corrector_currentInjuction: p_INorton = %d, %f, %f \n", p_bus_id, real(p_INorton), imag(p_INorton));
  
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GensalGenerator::corrector(
    double t_inc, bool flag)
{
  ///printf("\n***** GEN %d Corrector:\n", p_bus_id);
  if (getGenStatus()){
  
  if (p_hasExciter){	  
		p_exciter = getExciter();
		Efd = p_exciter->getFieldVoltage();
		// if (p_bus_id == 8022) Efd = Efdinit;  //use this one when relay is added
  }else{
		Efd = Efdinit;
  }

  if (p_hasGovernor){
	p_governor = getGovernor();
	Pmech = p_governor->getMechanicalPower();
	// if (p_bus_id == 8022) Pmech = Pmechinit;
  }else{
	Pmech = Pmechinit;
  }

// printf("Corrector: Efd = %f, Pmech = %f\n", Efd, Pmech); 

  double pi = 4.0*atan(1.0);
  double Psiq = x5Psiqpp_1 - Iq * Xdpp;
  double Psidpp = x3Eqp_1 * (Xdpp - Xl) / (Xdp - Xl) 
                + x4Psidp_1 * (Xdp - Xdpp) / (Xdp - Xl);
  double Psid = Psidpp - Id * Xdpp;
  double Telec = Psid * Iq - Psiq * Id;
  double TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl))
               * ((-x4Psidp_1 - (Xdp - Xl) * Id + x3Eqp_1));
  LadIfd = x3Eqp_1 * (1 + Sat(x3Eqp_1)) + (Xd - Xdp) * (Id + TempD);
  dx1d_1 = x2w_1 * 2 * pi * 60; // 60 represents the nominal frequency of 60 Hz
  //printf("H = %f, Pmech = %f, D = %f, x2w_1 = %f, Telec = %f\n", H, Pmech, D, x2w_1, Telec);
  dx2w_1 = 1 / (2 * H) * ((Pmech - D * x2w_1) / (1 + x2w_1) - Telec); //TBD: call Governor for Pmech
  dx3Eqp_1 = (Efd - LadIfd) / Tdop; //TBD: call Exciter for Efd and LadIfd
  dx4Psidp_1 = (-x4Psidp_1 - (Xdp - Xl) * Id + x3Eqp_1) / Tdopp;
  dx5Psiqpp_1 = (-x5Psiqpp_1 - (Xq - Xdpp) * Iq) / Tqopp;
          
  x1d_1 = x1d_0 + (dx1d_0 + dx1d_1) / 2.0 * t_inc;
  x2w_1 = x2w_0 + (dx2w_0 + dx2w_1) / 2.0 * t_inc;
  x3Eqp_1 = x3Eqp_0 + (dx3Eqp_0 + dx3Eqp_1) / 2.0 * t_inc;
  x4Psidp_1 = x4Psidp_0 + (dx4Psidp_0 + dx4Psidp_1) / 2.0 * t_inc;
  x5Psiqpp_1 = x5Psiqpp_0 + (dx5Psiqpp_0 + dx5Psiqpp_1) / 2.0 * t_inc;
  ///printf("gensal dx: %f\t%f\t%f\t%f\t%f\n", dx1d_1, dx2w_1, dx3Eqp_1, dx4Psidp_1, dx5Psiqpp_1);
  ///printf("gensal x: %f\t%f\t%f\t%f\t%f\n", x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqpp_1);
  
 // p_exciter->setOmega(x2w_1);
 //
   if (p_hasPss) {
	p_pss = getPss();
	p_pss->setOmega(x2w_1);
	p_pss->corrector(t_inc, flag);
	Vstab = p_pss->getVstab();
	}else{
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
  p_governor->setRotorSpeedDeviation(x2w_0);
  p_governor->corrector(t_inc, flag);
 }

  //if (p_bus_id == 1)
    //printf("\t%d          %12.6f   %12.6f   %12.6f   %12.6f   %12.6f\n",    
     //     p_bus_id, x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqpp_1);
  }else {
	x1d_0 = 0.0;
    x2w_0 = 0.0;
    x3Eqp_0 = 0.0;
    x4Psidp_0 = 0.0;
    x5Psiqpp_0 = 0.0;
	x1d_1 = 0.0;
    x2w_1 = 0.0;
    x3Eqp_1 = 0.0;
    x4Psidp_1 = 0.0;
    x5Psiqpp_1 = 0.0;
  
  }
}

void gridpack::dynamic_simulation::GensalGenerator::setWideAreaFreqforPSS(double freq)
{
	p_wideareafreq = freq;
	if (p_hasPss){
		//printf("-----!renke debug: GensalGenerator::setWideAreaFreqforPSS: %12.6f \n", freq);
		p_pss = getPss();
		p_pss->setWideAreaFreqforPSS(freq);	
	}
}

/**
 * Set voltage on each generator
 */
void gridpack::dynamic_simulation::GensalGenerator::setVoltage(
    gridpack::ComplexType voltage)
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
  if (!strcmp(signal,"standard")) {
    //sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f\n",
    //    p_bus_id,p_ckt.c_str(),real(p_mac_ang_s1),real(p_mac_spd_s1),real(p_mech),
    //    real(p_pelect));
    sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f	%12.6f\n",
          p_bus_id, p_ckt.c_str(), x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqpp_1);
    return true;
  } else if (!strcmp(signal,"init_debug")) {
    sprintf(string," %8d  %2s Something\n",p_bus_id,p_ckt.c_str());
    return true;
  } else if (!strcmp(signal,"debug_initial")) {
    sprintf(string,"Bus: %d PG: %f QG: %f\n",p_bus_id,p_pg,p_qg);
    return true;
  } else if (!strcmp(signal,"watch")) {
    if (getWatch()) {
      char buf[128];
//    sprintf(buf,", %f, %f",real(p_mac_ang_s1),real(p_mac_spd_s1));
      sprintf(string,",%8d, %2s, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f,",
          p_bus_id, p_ckt.c_str(), x1d_1, x2w_1+1, x3Eqp_1, x4Psidp_1, x5Psiqpp_1, Vterm, Efd, LadIfd);
      return true;
/*      if (strlen(buf) <= bufsize) {
        sprintf(string,"%s",buf);
        return true; 
      } else {
        printf ("size watch problem at gensal serialWrite() \n");
        return false;
      }
    } else {
        return false;
*/
    }
  } else if (!strcmp(signal,"debug_initial")) {
  return false;
  }
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
  if (getWatch()) {
    vals.push_back(x1d_1);
    vals.push_back(x2w_1+1.0);
  }
}
