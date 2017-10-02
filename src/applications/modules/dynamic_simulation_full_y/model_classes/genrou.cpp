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
#include "genrou.hpp"
//#include "exdc1.hpp"

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
    dx6Edp_1 = 0;;
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
  p_sbase = 100.0;

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
  if (!data->getValue(GENERATOR_TQOPP, &Tqopp, idx)) Tdopp=0.0; // Tqopp
  //if (!data->getValue(GENERATOR_S1, &S10, idx)) S10=0.0; // S10
  //if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.0; // S12
  if (!data->getValue(GENERATOR_S1, &S10, idx)) S10=0.17; // S10 TBD: check parser
  if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.55; // S12 TBD: check parser
  //printf("load S10 = %f, S12 = %f\n", S10, S12);
  if (!data->getValue(GENERATOR_XQP, &Xqp, idx)) Xqp=0.0; // Xqp
  //if (!data->getValue(GENERATOR_XQPP, &Xqp, idx)) Xqpp=0.0; // Xqpp // SJin: no GENERATOR_XQPP yet
  if (!data->getValue(GENERATOR_XDPP, &Xqp, idx)) Xqpp=0.0; // Xqpp // SJin: use Xdpp for compile
}

/**
 * Saturation function
 * @ param x
 */
double gridpack::dynamic_simulation::GenrouGenerator::Sat(double x)
{
    double a_ = S12 / S10 - 1.0 / 1.2;
    double b_ = -2 * S12 / S10 + 2;
    double c_ = S12 / S10 - 1.2;
    double A = (-b_ - sqrt(b_ * b_ - 4 * a_ * c_)) / (2 * a_);
    double B = S10 / ((1.0 - A) * (1.0 - A));
    double result = B * (x - A) * (x - A) / x;
    //printf("a = %f, b = %f, c = %f, A = %f, B = %f, S12 = %f, S10 = %f\n", a_, b_, c_, A, B, S12, S10);
    //printf("Sat result = %f\n", result); 
    return result; // Scaled Quadratic with 1.7.1 equations
}

/**
 * Initialize generator model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 */
void gridpack::dynamic_simulation::GenrouGenerator::init(double mag,
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
  x2w = 0;
  x1d = atan2(Viterm + Ir * Xq + Ii * Ra, Vrterm + Ir * Ra - Ii * Xq);
  Id = Ir * sin(x1d) - Ii * cos(x1d); // convert values to the dq axis
  Iq = Ir * cos(x1d) + Ii * sin(x1d); // convert values to the dq axis
  double Vr = Vrterm + Ra * Ir - Xdpp * Ii; // internal voltage on network reference
  double Vi = Viterm + Ra * Ii + Xdpp * Ir; // internal voltage on network reference
  // SJin: in pdf, its Vrint and Viint for Vd and Vq calculation?
  double Vd = Vr * sin(x1d) - Vi * sin(x1d); // convert to dq reference
  double Vq = Vr * cos(x1d) + Vi * cos(x1d); // convert to dq reference
  double Psiqpp = -Vd / (1 + x2w);
  double Psidpp = + Vq / (1 + x2w);
  x4Psidp = Psidpp - Id * (Xdpp - Xl);
  x3Eqp = x4Psidp + Id * (Xdp - Xl);
  x6Edp = Iq * (Xq - Xqp); 
  x5Psiqp = x6Edp + Iq * (Xqp - Xl);
  //printf("Xdpp = %f, Xq = %f, Iq = %f\n", Xdpp, Xq, Iq);
  Efd = x3Eqp * (1 + Sat(x3Eqp)) + Id * (Xd - Xdp);
  //printf("x3Eqp = %f, Id = %f, Xd = %f, Xdp = %f\n", x3Eqp, Id, Xd, Xdp);
  LadIfd = Efd;
  Pmech = Psidpp * Iq - Psiqpp * Id;

  printf("genrou init: %f\t%f\t%f\t%f\t%f\t%f\n", x1d, x2w, x3Eqp, x4Psidp, x5Psiqp, x6Edp);
  printf("gensal init: Efd = %f, Pmech = %f\n", Efd, Pmech);
  //Efdinit = Efd;
  //Pmechinit = Pmech;

  p_exciter = getExciter();
  p_exciter->setFieldVoltage(Efd);
  p_exciter->init(mag, ang, ts);

  p_governor = getGovernor();
  p_governor->setMechanicalPower(Pmech);
  p_governor->setRotorSpeedDeviation(x2w); // set Speed Deviation w for wsieg1 
  p_governor->init(mag, ang, ts);

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
void gridpack::dynamic_simulation::GenrouGenerator::predictor_currentInjection(bool flag)
{
  if (!flag) {
    x1d = x1d_1;
    x2w = x2w_1;
    x3Eqp = x3Eqp_1;
    x4Psidp = x4Psidp_1;
    x5Psiqp = x5Psiqp_1;
    x6Edp = x6Edp_1; 
  }  
  // Calculate INorton_full
  // Admittance
  B = -Xdpp / (Ra * Ra + Xdpp * Xdpp);
  G = Ra / (Ra * Ra + Xdpp * Xdpp);
  //printf("Xdpp = %f, Ra = %f, B = %f, G = %f\n", Xdpp, Ra, B, G);
  // Setup
  double Psiqpp = - x6Edp * (Xqpp - Xl) / (Xqp - Xl) - x5Psiqp * (Xqp - Xqpp) / (Xqp - Xl); 
  double Psidpp = + x3Eqp * (Xdpp - Xl) / (Xdp - Xl) + x4Psidp * (Xdp - Xdpp) / (Xdp - Xl); 
  double Vd = - Psiqpp * (1 + x2w);
  double Vq = + Psidpp * (1 + x2w);
  Vterm = presentMag;
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  double Vdterm = Vrterm * sin(x1d) - Viterm * cos(x1d);
  double Vqterm = Vrterm * cos(x1d) + Viterm * sin(x1d);
  //printf("x2w = %f, Vd = %f, Vq = %f, Vrterm = %f, Viterm = %f, Vdterm = %f, Vqterm = %f, Theta = %f\n", x2w, Vd, Vq, Vrterm, Viterm, Vdterm, Vqterm, Theta);
  //DQ Axis
  Id = (Vd - Vdterm) * G - (Vq - Vqterm) * B;
  Iq = (Vd - Vdterm) * B + (Vq - Vqterm) * G;
  double Idnorton = Vd * G - Vq * B;
  double Iqnorton = Vd * B + Vq * G;
  //Network
  Ir = + Id * sin(x1d) + Iq * cos(x1d);
  Ii = - Id * cos(x1d) + Iq * sin(x1d);
  IrNorton = + Idnorton * sin(x1d) + Iqnorton * cos(x1d);
  IiNorton = - Idnorton * cos(x1d) + Iqnorton * sin(x1d);
  IrNorton = IrNorton * MVABase / p_sbase; 
  IiNorton = IiNorton * MVABase / p_sbase; 
  //gridpack::ComplexType INorton(IrNorton, IiNorton);
  p_INorton = gridpack::ComplexType(IrNorton, IiNorton);
} 

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GenrouGenerator::predictor(
    double t_inc, bool flag)
{
  printf("\n***** GEN %d Predicator:\n", p_bus_id);

  p_exciter = getExciter();
  Efd = p_exciter->getFieldVoltage();
  //Efd = Efdinit;

  p_governor = getGovernor();
  Pmech = p_governor->getMechanicalPower();
  //Pmech = Pmechinit;
 
  printf("Efd = %f, Pmech = %f\n", Efd, Pmech); 

  if (!flag) {
    x1d = x1d_1;
    x2w = x2w_1;
    x3Eqp = x3Eqp_1;
    x4Psidp = x4Psidp_1;
    x5Psiqp = x5Psiqp_1;
    x6Edp = x6Edp_1; 
  }  
    
  double pi = 4.0*atan(1.0);
  double Psiqpp = - x6Edp * (Xqpp - Xl) / (Xqp - Xl) - x5Psiqp * (Xqp - Xqpp) / (Xqp - Xl); 
  double Psidpp = + x3Eqp * (Xdpp - Xl) / (Xdp - Xl) + x4Psidp * (Xdp - Xdpp) / (Xdp - Xl); 
  double Telec = Psidpp * Iq - Psiqpp * Id;
  double TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl))
               * (-x4Psidp - (Xdp - Xl) * Id + x3Eqp);
  LadIfd = x3Eqp * (1 + Sat(x3Eqp)) + (Xd - Xdp) * (Id + TempD); // update Ifd later
  dx1d = x2w * 2 * pi * 60; // 60 represents the nominal frequency of 60 Hz
  //printf("H = %f, Pmech = %f, D = %f, x2w = %f, Telec = %f\n", H, Pmech, D, x2w, Telec);
  dx2w = 1 / (2 * H) * ((Pmech - D * x2w) / (1 + x2w) - Telec); //TBD: call Governor for Pmech (Done)
  dx3Eqp = (Efd - LadIfd) / Tdop; //TBD: call Exciter for Efd and LadIfd (Done)
  dx4Psidp = (-x4Psidp - (Xdp - Xl) * Id + x3Eqp) / Tdopp;
  dx5Psiqp = (-x5Psiqp + (Xqp - Xl) * Iq + x6Edp) / Tqopp;
  double TempQ = (Xqp - Xqpp) / ((Xqp - Xl) * (Xqp - Xl))
               * (-x5Psiqp + (Xqp - Xl) * Iq + x6Edp);
  dx6Edp = (-x6Edp + (Xq - Xqp) * (Iq - TempQ)) / Tqopp;  // SJin: in pdf, its Tqop, I use Tqopp?

  x1d_1 = x1d + dx1d * t_inc;
  x2w_1 = x2w + dx2w * t_inc;
  x3Eqp_1 = x3Eqp + dx3Eqp * t_inc;
  x4Psidp_1 = x4Psidp + dx4Psidp * t_inc;
  x5Psiqp_1 = x5Psiqp + dx5Psiqp * t_inc;
  x6Edp_1 = x6Edp + dx6Edp * t_inc;
  printf("genrou dx: %f\t%f\t%f\t%f\t%f\t%f\n", dx1d, dx2w, dx3Eqp, dx4Psidp, dx5Psiqp, x6Edp);
  printf("genrou x: %f\t%f\t%f\t%f\t%f\t%f\n", x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqp_1, x6Edp_1);
  
  p_exciter->setOmega(x2w_1);
  p_exciter->setVterminal(presentMag);
  p_exciter->predictor(t_inc, flag);

  p_governor->setRotorSpeedDeviation(x2w);
  p_governor->predictor(t_inc, flag);
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
  double Psiqpp = - x6Edp * (Xqpp - Xl) / (Xqp - Xl) - x5Psiqp * (Xqp - Xqpp) / (Xqp - Xl); 
  double Psidpp = + x3Eqp * (Xdpp - Xl) / (Xdp - Xl) + x4Psidp * (Xdp - Xdpp) / (Xdp - Xl); 
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
  p_INorton = gridpack::ComplexType(IrNorton, IiNorton);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GenrouGenerator::corrector(
    double t_inc, bool flag)
{
  printf("\n***** GEN %d Corrector:\n", p_bus_id);

  p_exciter = getExciter();
  Efd = p_exciter->getFieldVoltage(); 
  //Efd = Efdinit;

  p_governor = getGovernor();
  Pmech = p_governor->getMechanicalPower();
  //Pmech = Pmechinit; 

  printf("Efd = %f, Pmech = %f\n", Efd, Pmech); 

  double pi = 4.0*atan(1.0);
  double Psiqpp = - x6Edp_1 * (Xqpp - Xl) / (Xqp - Xl) - x5Psiqp_1 * (Xqp - Xqpp) / (Xqp - Xl); 
  double Psidpp = + x3Eqp_1 * (Xdpp - Xl) / (Xdp - Xl) + x4Psidp_1 * (Xdp - Xdpp) / (Xdp - Xl); 
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
  dx6Edp_1 = (-x6Edp_1 + (Xq - Xqp) * (Iq - TempQ)) / Tqopp; // SJin: Tqop -> Tqopp?

  x1d_1 = x1d + (dx1d + dx1d_1) / 2.0 * t_inc;
  x2w_1 = x2w + (dx2w + dx2w_1) / 2.0 * t_inc;
  x3Eqp_1 = x3Eqp + (dx3Eqp + dx3Eqp_1) / 2.0 * t_inc;
  x4Psidp_1 = x4Psidp + (dx4Psidp + dx4Psidp_1) / 2.0 * t_inc;
  x5Psiqp_1 = x5Psiqp + (dx5Psiqp + dx5Psiqp_1) / 2.0 * t_inc;
  x6Edp_1 = x6Edp + (dx6Edp + dx6Edp_1) / 2.0 * t_inc;
  printf("genrou dx: %f\t%f\t%f\t%f\t%f\t%f\n", dx1d_1, dx2w_1, dx3Eqp_1, dx4Psidp_1, dx5Psiqp_1, dx6Edp);
  printf("genrou x: %f\t%f\t%f\t%f\t%f\t%f\n", x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqp_1, x6Edp);
  
  p_exciter->setOmega(x2w_1);
  p_exciter->setVterminal(presentMag);
  p_exciter->corrector(t_inc, flag);
 
  p_governor->setRotorSpeedDeviation(x2w);
  p_governor->corrector(t_inc, flag);

  //if (p_bus_id == 1)
    //printf("\t%d          %12.6f   %12.6f   %12.6f   %12.6f   %12.6f	%12.6f\n",    
     //     p_bus_id, x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqp_1, x6Edp_1);
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

/**
 * return a vector containing any generator values that are being
 * watched
 * @param vals vector of watched values
 */
void gridpack::dynamic_simulation::GenrouGenerator::getWatchValues(
    std::vector<double> &vals)
{
  vals.clear();
  if (getWatch()) {
    vals.push_back(x1d_1+1.0);
    vals.push_back(x2w_1);
  }
}
