
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   genrou.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   Oct 9, 2015
 * 
 * @Modified: November 21, 2022, Shri, Disable saturation if S10 and S12 are 0.
 *
 * @Modified: November 27, 2022, Shri, Fixed the model to validate against PSSE
 *
 * @Modified: Dec 9, 2022, Shri, print voltage and generator power
 * @brief  : Round rotor generator model
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
#include "genrou.hpp"


/**
 * ROTATE - Rotates a 2-d vector Fa + iFb by angle 'ang' to give the rotated
 *          vector Fd + iFq
 */
void rotate(double Fa, double Fb, double ang, double *Fd, double *Fq)
{
  *Fd = Fa*cos(ang) - Fb*sin(ang);
  *Fq = Fa*sin(ang) + Fb*cos(ang);
}

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::GenrouGenerator::GenrouGenerator(void)
{
    dx1d = 0;
    dx2w = 0;
    dx3Eqp = 0;
    dx4Psidp = 0;
    dx5Psiqp = 0;;
    dx6Edp = 0;;
    dx1d_1 = 0;
    dx2w_1 = 0;
    dx3Eqp_1 = 0;
    dx4Psidp_1 = 0;
    dx5Psiqp_1 = 0;;
    dx6Edp_1 = 0;

    p_tripped = false;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::GenrouGenerator::~GenrouGenerator(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::GenrouGenerator::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if(!data->getValue(CASE_SBASE,&p_sbase)) p_sbase = 100.0;
  data->getValue(BUS_NUMBER,&p_bus_id);
  data->getValue(GENERATOR_ID,&p_ckt,idx);
  if (!data->getValue(GENERATOR_PG, &p_pg,idx)) p_pg = 0.0;
  if (!data->getValue(GENERATOR_QG, &p_qg,idx)) p_qg = 0.0;
  if (!data->getValue(GENERATOR_STAT, &p_status,idx)) p_status = 0;
  //printf("load p_pg = %f, p_qg = %f\n", p_pg, p_qg);
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
  Xqpp = Xdpp;
  if (!data->getValue(GENERATOR_XL, &Xl, idx)) Xl=0.0; // Xl
  if (!data->getValue(GENERATOR_TDOP, &Tdop, idx)) Tdop=0.0; // Tdop
  if (!data->getValue(GENERATOR_TDOPP, &Tdopp, idx)) Tdopp=0.0; // Tdopp
  if (!data->getValue(GENERATOR_TQOPP, &Tqopp, idx)) Tqopp=0.0; // Tqopp
  if (!data->getValue(GENERATOR_S1, &S10, idx)) S10=0.05; // S10 TBD: check parser
  if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.3; // S12 TBD: check parser
  if (!data->getValue(GENERATOR_XQP, &Xqp, idx)) Xqp=0.0; // Xqp

  if (!data->getValue(GENERATOR_TQOP, &Tqop, idx)) Tqop=0.75; // Tqop

  if(fabs(S10*S12) < 1e-6) {
    // Zero saturation
    enableSat = false;
  } else enableSat = true;
  printFlag = false;

  double tmp = sqrt(p_pg*p_pg +p_qg*p_qg);
  // Increase Machine base to 1.2*Sgen if Sgen > MBase
  if ( tmp > MBase) {
    MBase = tmp*1.2;
  }
  
}

/**
 * Saturation function
 * @ param x
 */
double gridpack::dynamic_simulation::GenrouGenerator::Sat(double x)
{
  if (enableSat) {
    // the following is another method for saturation computation, add by renke
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
    
    return result; // Scaled Quadratic with 1.7.1 equations
  } else {
    return 0.0;
  }
}

/**
 * Initialize generator model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 */
void gridpack::dynamic_simulation::GenrouGenerator::init(double mag,
    double ang, double ts)
{
  double pi = 4.0*atan(1.0);

  Vterm = mag;
  presentMag = mag;
  Theta = ang;
  presentAng = ang;

  // Generator P and Q in pu on Machine base (MBase)
  double P = p_pg/ MBase;
  double Q = p_qg/ MBase;
  
  genP = P;
  genQ = Q;

  // Terminal voltage in network reference frame
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);

  // Terminal current in network reference frame. Since P and Q
  // are on MBase, Ir and Ii are also on MBase
  Ir = (P * Vrterm + Q * Viterm) / (Vterm * Vterm);
  Ii = (P * Viterm - Q * Vrterm) / (Vterm * Vterm);

  // Generator state variable initialization

  // Speed deviation
  x2w = 0;

  // Machine angle
  x1d = atan2(Viterm + Ir * Xq + Ii * Ra, Vrterm + Ir * Ra - Ii * Xq);

  // axis of rotation is q-axis (lagging behind d-axis by 90 degrees)
  double theta = pi/2.0 - x1d;
  
  // Generator currents in machine dq axis reference frame
  rotate(Ir,Ii,theta, &Id, &Iq);

  // Generator internal voltage in network reference frame
  double Vr = Vrterm + Ra * Ir - Xdpp * Ii; // internal voltage on network reference
  double Vi = Viterm + Ra * Ii + Xdpp * Ir; // internal voltage on network reference

  // Generator internal voltage in machine reference frame
  double Vd, Vq;
  rotate(Vr,Vi,theta, &Vd, &Vq);
  
  Vd = Vr * sin(x1d) - Vi * cos(x1d); // convert to dq reference
  Vq = Vr * cos(x1d) + Vi * sin(x1d); // convert to dq reference
  
  double Psiqpp = -Vd;
  double Psidpp = + Vq;

  // q-axis transient voltage
  x3Eqp = Vq + (Xdp - Xdpp)*Id;

  // d-axis transient voltage
  x6Edp = Vd - (Xqp - Xqpp)*Iq;

  // d-axis flux
  x4Psidp = x3Eqp - (Xdp - Xl)*Id;

  // q-axis flux
  x5Psiqp = x6Edp + (Xqp - Xl)*Iq;
  
  // Field voltage
  Efd = x3Eqp * (1 + Sat(x3Eqp)) + Id * (Xd - Xdp);
  
  // Field current
  LadIfd = Efd;

  // Mechanical power
  Pmech = Psidpp * Iq - Psiqpp * Id;

  Efdinit = Efd;
  Pmechinit = Pmech;

  // Initialize exciter
  if (p_hasExciter) {
    p_exciter = getExciter();
    p_exciter->setVterminal(Vterm);
    p_exciter->setVcomp(mag); 
    p_exciter->setFieldVoltage(Efd);
    p_exciter->setFieldCurrent(LadIfd);
    p_exciter->init(mag, ang, ts);
  }

  // Initialize governor
  if (p_hasGovernor) {
    p_governor = getGovernor();
    p_governor->setMechanicalPower(Pmech);
    p_governor->setRotorSpeedDeviation(x2w); // set Speed Deviation w for wsieg1 
    p_governor->init(mag, ang, ts);
  }

  // Norton impedance
  p_Norton_Ya = NortonImpedence();

  x1d_1     = x1d;
  x2w_1     = x2w;
  x3Eqp_1   = x3Eqp;
  x4Psidp_1 = x4Psidp;
  x5Psiqp_1 = x5Psiqp;
  x6Edp_1   = x6Edp;
}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType gridpack::dynamic_simulation::GenrouGenerator::INorton()
{
  return p_INorton;
}


/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType gridpack::dynamic_simulation::GenrouGenerator::NortonImpedence()
{
  double ra = Ra;
  double xd = Xdpp;
  B = -xd / (ra * ra + xd * xd);
  G = ra / (ra * ra + xd * xd);

  // Conversion from machine base to system base
  B *= MBase/p_sbase;
  G *= MBase/p_sbase;

  gridpack::ComplexType Y_a(G, B);
  return Y_a;
}

/**
 * Predict part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GenrouGenerator::predictor_currentInjection(bool flag)
{
  // Calculate INorton_full
  // Admittance
  B = -Xdpp / (Ra * Ra + Xdpp * Xdpp);
  G = Ra / (Ra * Ra + Xdpp * Xdpp);

  // Setup
  double Psiqpp = - x6Edp * (Xqpp - Xl) / (Xqp - Xl) - x5Psiqp * (Xqp - Xqpp) / (Xqp - Xl); 
  double Psidpp = + x3Eqp * (Xdpp - Xl) / (Xdp - Xl) + x4Psidp * (Xdp - Xdpp) / (Xdp - Xl);

  double Vd = - Psiqpp; //* (1 + x2w);
  double Vq = + Psidpp; //* (1 + x2w);

  Vterm = presentMag;
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);

  
  gridpack::ComplexType vt_complex_tmp = gridpack::ComplexType(Vrterm, Viterm); 

  double Idnorton = Vd * G - Vq * B;
  double Iqnorton = Vd * B + Vq * G;

  //Network
  Ir = + Id * sin(x1d) + Iq * cos(x1d);
  Ii = - Id * cos(x1d) + Iq * sin(x1d);

  genP = Vrterm*Ir + Viterm*Ii;
  genQ = Viterm*Ir - Vrterm*Ii;

  
  IrNorton = + Idnorton * sin(x1d) + Iqnorton * cos(x1d);
  IiNorton = - Idnorton * cos(x1d) + Iqnorton * sin(x1d);

  IrNorton = IrNorton * MBase / p_sbase; 
  IiNorton = IiNorton * MBase / p_sbase; 
  
  if (getGenStatus()){
	  if (p_tripped){
		  p_INorton = p_Norton_Ya*vt_complex_tmp;
	  }else{
		  p_INorton = gridpack::ComplexType(IrNorton, IiNorton);	
	  }		
  }else {
	  p_INorton = gridpack::ComplexType(0.0, 0.0);
  }
  
} 

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GenrouGenerator::predictor(
    double t_inc, bool flag)
{

  if (getGenStatus()){
	  
  if (p_hasExciter) {
    p_exciter = getExciter();
    p_exciter->setOmega(x2w);
    p_exciter->setVterminal(presentMag);
    p_exciter->setVcomp(presentMag);
    p_exciter->setFieldCurrent(LadIfd);

    Efd = p_exciter->getFieldVoltage();
  } else {
    Efd = Efdinit;
  }
  
  if (p_hasGovernor) {
    p_governor = getGovernor();
    p_governor->setRotorSpeedDeviation(x2w);
    Pmech = p_governor->getMechanicalPower();
  } else {
    Pmech = Pmechinit;
  }
  
  double pi = 4.0*atan(1.0);
  double Psiqpp = - x6Edp * (Xqpp - Xl) / (Xqp - Xl) - x5Psiqp * (Xqp - Xqpp) / (Xqp - Xl); 
  double Psidpp = + x3Eqp * (Xdpp - Xl) / (Xdp - Xl) + x4Psidp * (Xdp - Xdpp) / (Xdp - Xl);

  double Vd = - Psiqpp;
  double Vq = + Psidpp;

  Vterm = presentMag;
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  double Vdterm = Vrterm * sin(x1d) - Viterm * cos(x1d);
  double Vqterm = Vrterm * cos(x1d) + Viterm * sin(x1d);

  //DQ Axis currents
  Id = (Vd - Vdterm) * G - (Vq - Vqterm) * B;
  Iq = (Vd - Vdterm) * B + (Vq - Vqterm) * G;

  
  double Telec = Psidpp * Iq - Psiqpp * Id;
  double TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl))
               * (-x4Psidp - (Xdp - Xl) * Id + x3Eqp);
  LadIfd = x3Eqp * (1 + Sat(x3Eqp)) + (Xd - Xdp) * (Id + TempD); // update Ifd later
  //printf("Psiqpp=%f,Psidpp=%f,Telec=%f,TempD=%f,LadIfd=%f\n",Psiqpp,Psidpp,Telec,TempD,LadIfd);
  //printf("Id=%f, Iq=%f\n", Id, Iq);

  dx1d = x2w * 2 * pi * 60; // 60 represents the nominal frequency of 60 Hz
  //printf("H = %f, Pmech = %f, D = %f, x2w = %f, Telec = %f\n", H, Pmech, D, x2w, Telec);
  dx2w = 1 / (2 * H) * ((Pmech - D * x2w) / (1 + x2w) - Telec); //TBD: call Governor for Pmech (Done)
  dx3Eqp = (Efd - LadIfd) / Tdop; //TBD: call Exciter for Efd and LadIfd (Done)
  dx4Psidp = (-x4Psidp - (Xdp - Xl) * Id + x3Eqp) / Tdopp;
  dx5Psiqp = (-x5Psiqp + (Xqp - Xl) * Iq + x6Edp) / Tqopp;
  double TempQ = (Xqp - Xqpp) / ((Xqp - Xl) * (Xqp - Xl))
               * (-x5Psiqp + (Xqp - Xl) * Iq + x6Edp);
  //dx6Edp = (-x6Edp + (Xq - Xqp) * (Iq - TempQ)) / Tqopp;  // SJin: in pdf, its Tqop, I use Tqopp?
  // Yuan modified below 20201002
  dx6Edp = (-x6Edp + (Xq - Xqp) * (Iq - TempQ)) / Tqop;  // SJin: in pdf, its Tqop, I use Tqopp?
  // Yuan modified end

  x1d_1 = x1d + dx1d * t_inc;
  x2w_1 = x2w + dx2w * t_inc;
  x3Eqp_1 = x3Eqp + dx3Eqp * t_inc;
  x4Psidp_1 = x4Psidp + dx4Psidp * t_inc;
  x5Psiqp_1 = x5Psiqp + dx5Psiqp * t_inc;
  x6Edp_1 = x6Edp + dx6Edp * t_inc;
  
  if (printFlag) {
	  printf("genrou dx: %f\t%f\t%f\t%f\t%f\t%f\n", dx1d, dx2w, dx3Eqp, dx4Psidp, dx5Psiqp, x6Edp);
	  printf("genrou x: %f\t%f\t%f\t%f\t%f\t%f\n", x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqp_1, x6Edp_1);
  }
  
  if (p_hasExciter) {
    p_exciter->predictor(t_inc, flag);
  }

  if (p_hasGovernor) {
    p_governor->predictor(t_inc, flag);
  }
  
  	if (p_tripped){
		x1d = 0.0;
		x2w = -1.0;
		x3Eqp = 0.0;
		x4Psidp = 0.0;
		x5Psiqp = 0.0;
		x6Edp = 0.0;
		x1d_1 = 0.0;
		x2w_1 = -1.0;
		x3Eqp_1 = 0.0;
		x4Psidp_1 = 0.0;
		x5Psiqp_1 = 0.0;	
		x6Edp_1 = 0.0;
		genP = 0.0;
		genQ = 0.0;
	}
  }else {
	x1d = 0.0;
    x2w = -1.0;
    x3Eqp = 0.0;
    x4Psidp = 0.0;
    x5Psiqp = 0.0;
	x6Edp = 0.0;
	x1d_1 = 0.0;
    x2w_1 = -1.0;
    x3Eqp_1 = 0.0;
    x4Psidp_1 = 0.0;
    x5Psiqp_1 = 0.0;
	x6Edp_1 = 0.0;
	genP = 0.0;
	genQ = 0.0;
  
  }
  
}

/**
 * Corrector part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GenrouGenerator::corrector_currentInjection(bool flag)
{
  // Calculate INorton_full
  // Admittance
  B = -Xdpp / (Ra * Ra + Xdpp * Xdpp);
  G = Ra / (Ra * Ra + Xdpp * Xdpp);
  //printf("B = %f, G = %f\n", B, G);
  // Setup
  double Psiqpp = - x6Edp_1 * (Xqpp - Xl) / (Xqp - Xl) - x5Psiqp_1 * (Xqp - Xqpp) / (Xqp - Xl); 
  double Psidpp = + x3Eqp_1 * (Xdpp - Xl) / (Xdp - Xl) + x4Psidp_1 * (Xdp - Xdpp) / (Xdp - Xl);
  
  double Vd = -Psiqpp;
  double Vq = +Psidpp;
  
  Vterm = presentMag;
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  double Vdterm = Vrterm * sin(x1d_1) - Viterm * cos(x1d_1);
  double Vqterm = Vrterm * cos(x1d_1) + Viterm * sin(x1d_1);
  
  gridpack::ComplexType vt_complex_tmp = gridpack::ComplexType(Vrterm, Viterm); 
  
  double Idnorton = Vd * G - Vq * B;
  double Iqnorton = Vd * B + Vq * G;

  IrNorton = + Idnorton * sin(x1d_1) + Iqnorton * cos(x1d_1);
  IiNorton = - Idnorton * cos(x1d_1) + Iqnorton * sin(x1d_1);
  
  IrNorton = IrNorton * MBase / p_sbase; 
  IiNorton = IiNorton * MBase / p_sbase; 
  //gridpack::ComplexType INorton(IrNorton, IiNorton);
  //p_INorton = gridpack::ComplexType(IrNorton, IiNorton);
  
  if (getGenStatus()){	  
      if (p_tripped){
		  p_INorton = p_Norton_Ya*vt_complex_tmp;
	  }else{
		  p_INorton = gridpack::ComplexType(IrNorton, IiNorton);		
	  }	 	  
  }else {
	  p_INorton = gridpack::ComplexType(0.0, 0.0);
  }
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GenrouGenerator::corrector(
    double t_inc, bool flag)
{
 
  if (getGenStatus()){
	  
  if (p_hasExciter) {
    p_exciter = getExciter();
    p_exciter->setOmega(x2w_1);
    p_exciter->setVterminal(presentMag);
    p_exciter->setVcomp(presentMag);
    p_exciter->setFieldCurrent(LadIfd);
    
    Efd = p_exciter->getFieldVoltage();
  } else {
    Efd = Efdinit;
  }
  
  if (p_hasGovernor) {
    p_governor = getGovernor();
    p_governor->setRotorSpeedDeviation(x2w_1);
    Pmech = p_governor->getMechanicalPower();
  } else {
    Pmech = Pmechinit;
  }

  double pi = 4.0*atan(1.0);
  double Psiqpp = - x6Edp_1 * (Xqpp - Xl) / (Xqp - Xl) - x5Psiqp_1 * (Xqp - Xqpp) / (Xqp - Xl); 
  double Psidpp = + x3Eqp_1 * (Xdpp - Xl) / (Xdp - Xl) + x4Psidp_1 * (Xdp - Xdpp) / (Xdp - Xl);

  double Vd = - Psiqpp;
  double Vq = + Psidpp;

  
  Vterm = presentMag;
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  double Vdterm = Vrterm * sin(x1d_1) - Viterm * cos(x1d_1);
  double Vqterm = Vrterm * cos(x1d_1) + Viterm * sin(x1d_1);

  //DQ Axis
  Id = (Vd - Vdterm) * G - (Vq - Vqterm) * B;
  Iq = (Vd - Vdterm) * B + (Vq - Vqterm) * G;
  
  double Telec = Psidpp * Iq - Psiqpp * Id;
  double TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl))
               * (-x4Psidp_1 - (Xdp - Xl) * Id + x3Eqp_1);
  LadIfd = x3Eqp_1 * (1 + Sat(x3Eqp_1)) + (Xd - Xdp) * (Id + TempD); // update Ifd later
  dx1d_1 = x2w_1 * 2 * pi * 60; // 60 represents the nominal frequency of 60 Hz
  //printf("H = %f, Pmech = %f, D = %f, x2w_1 = %f, Telec = %f\n", H, Pmech, D, x2w_1, Telec);
  dx2w_1 = 1 / (2 * H) * ((Pmech - D * x2w_1) / (1 + x2w_1) - Telec); //TBD: call Governor for Pmech (Done)
  dx3Eqp_1 = (Efd - LadIfd) / Tdop; //TBD: call Exciter for Efd and LadIfd (Done)
  dx4Psidp_1 = (-x4Psidp_1 - (Xdp - Xl) * Id + x3Eqp_1) / Tdopp;
  dx5Psiqp_1 = (-x5Psiqp_1 + (Xqp - Xl) * Iq + x6Edp_1) / Tqopp;
  double TempQ = (Xqp - Xqpp) / ((Xqp - Xl) * (Xqp - Xl))
               * (-x5Psiqp_1 + (Xqp - Xl) * Iq + x6Edp_1);

  dx6Edp_1 = (-x6Edp_1 + (Xq - Xqp) * (Iq - TempQ)) / Tqop;

  x1d = x1d + (dx1d + dx1d_1) / 2.0 * t_inc;
  x2w = x2w + (dx2w + dx2w_1) / 2.0 * t_inc;
  x3Eqp = x3Eqp + (dx3Eqp + dx3Eqp_1) / 2.0 * t_inc;
  x4Psidp = x4Psidp + (dx4Psidp + dx4Psidp_1) / 2.0 * t_inc;
  x5Psiqp = x5Psiqp + (dx5Psiqp + dx5Psiqp_1) / 2.0 * t_inc;
  x6Edp = x6Edp + (dx6Edp + dx6Edp_1) / 2.0 * t_inc;

  if (p_hasExciter) {
    p_exciter->corrector(t_inc, flag);
  }

  if (p_hasGovernor) {
    p_governor->corrector(t_inc, flag);
  }

 	if (p_tripped){
		x1d = 0.0;
		x2w = -1.0;
		x3Eqp = 0.0;
		x4Psidp = 0.0;
		x5Psiqp = 0.0;
		x6Edp = 0.0;
		x1d_1 = 0.0;
		x2w_1 = -1.0;
		x3Eqp_1 = 0.0;
		x4Psidp_1 = 0.0;
		x5Psiqp_1 = 0.0;
		x6Edp_1 = 0.0;
		genP = 0.0;
		genQ = 0.0;	
	}

  }else {
	x1d = 0.0;
    x2w = -1.0;
    x3Eqp = 0.0;
    x4Psidp = 0.0;
    x5Psiqp = 0.0;
	x6Edp = 0.0;
	x1d_1 = 0.0;
    x2w_1 = -1.0;
    x3Eqp_1 = 0.0;
    x4Psidp_1 = 0.0;
    x5Psiqp_1 = 0.0;
	x6Edp_1 = 0.0;
	genP = 0.0;
	genQ = 0.0; 
  }
  
}

bool gridpack::dynamic_simulation::GenrouGenerator::tripGenerator()
{
	p_tripped = true;
	
	return true;
}

/**
 * Set voltage on each generator
 */
void gridpack::dynamic_simulation::GenrouGenerator::setVoltage(
    gridpack::ComplexType voltage)
{
  presentMag = abs(voltage);
  presentAng = atan2(imag(voltage), real(voltage));  
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double gridpack::dynamic_simulation::GenrouGenerator::getFieldVoltage()
{
  return Efd;
}


/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
// Yuan commented out the below 20201011
/*
void gridpack::dynamic_simulation::GenrouGenerator::write(
    const char* signal, char *string)
{
  if (!strcmp(signal,"standard")) {
    //sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f\n",
    //    p_bus_id,p_ckt.c_str(),real(p_mac_ang_s1),real(p_mac_spd_s1),real(p_mech),
    //    real(p_pelect));
    sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f	%12.6f	%12.6f\n",
          p_bus_id, p_ckt.c_str(), x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqp_1, x6Edp_1);
  }
}
*/
// Yuan commented end

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
bool gridpack::dynamic_simulation::GenrouGenerator::serialWrite(
    char* string, const int bufsize, const char *signal)
{
  bool ret = false;
  if (!strcmp(signal,"standard")) {
    sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f	%12.6f  %12.6f\n",
          p_bus_id, p_ckt.c_str(), x1d_1, x2w_1+1.0, x3Eqp_1, x4Psidp_1, x5Psiqp_1, x6Edp_1);
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
void gridpack::dynamic_simulation::GenrouGenerator::getWatchValues(
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
bool gridpack::dynamic_simulation::GenrouGenerator::setState(std::string name,
    double value)
{
  if(name == "ANGLE") {
    x1d = x1d_1 = value;
    return true;
  } else if(name == "SPEED_DEV") {
    x2w = x2w_1 = value;
    return true;
  } else {
    return false;
  }
}

/**
 * Get internal state parameter in generator
 * @param name character string corresponding to state variable
 * @return value current value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::GenrouGenerator::getState(std::string name,
    double *value)
{
  if(name == "ANGLE") {
    *value = x1d;
    return true;
  } else if(name == "SPEED_DEV") {
    *value = x2w;
    return true;
  } else {
    return false;
  }
}
