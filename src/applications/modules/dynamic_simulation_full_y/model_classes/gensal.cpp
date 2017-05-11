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

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/include/gridpack.hpp"
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
  dx5Psiqpp_1 = 0;;
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

  data->getValue(BUS_NUMBER,&p_bus_id);
  data->getValue(GENERATOR_ID,&p_ckt,idx);
  if (!data->getValue(GENERATOR_PG, &p_pg,idx)) p_pg = 0.0;
  if (!data->getValue(GENERATOR_QG, &p_qg,idx)) p_qg = 0.0;
  if (!data->getValue(GENERATOR_STAT, &p_status,idx)) p_status = 0;
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
  if (!data->getValue(GENERATOR_S1, &S10, idx)) S10=0.17; // S10 TBD: check parser
  if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.55; // S12 TBD: check parser
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
  double P = p_pg / MVABase;
  double Q = p_qg / MVABase;
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

  p_exciter = getExciter();
  p_exciter->setVterminal(Vterm); 
  //  p_exciter->setVcomp(abs(Vterm)); //TBD Need to updated later to calculate Vcomp 
  p_exciter->setVcomp(mag); 
  p_exciter->setFieldVoltage(Efd);
  p_exciter->setFieldCurrent(LadIfd);
  p_exciter->init(mag, ang, ts);

  p_governor = getGovernor();
  p_governor->setMechanicalPower(Pmech);
  p_governor->setRotorSpeedDeviation(x2w_0); // set Speed Deviation w for wsieg1 
  p_governor->init(mag, ang, ts);
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
  double ra = Ra * p_sbase / MVABase;
  double xd = Xdpp * p_sbase / MVABase;
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
  double Vd = -Psiqpp * (1 + x2w_0);
  double Vq = +Psidpp * (1 + x2w_0);
  Vterm = presentMag;
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  double Vdterm = Vrterm * sin(x1d_0) - Viterm * cos(x1d_0);
  double Vqterm = Vrterm * cos(x1d_0) + Viterm * sin(x1d_0);
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
  if (getGenStatus()){
    p_INorton = gridpack::ComplexType(IrNorton, IiNorton);	  
  }else {
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
  if (getGenStatus()){
    p_exciter = getExciter();
    Efd = p_exciter->getFieldVoltage();

    p_governor = getGovernor();
    Pmech = p_governor->getMechanicalPower();

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
    dx2w_0 = 1 / (2 * H) * ((Pmech - D * x2w_0) / (1 + x2w_0) - Telec); //TBD: call Governor for Pmech (Done)
    dx3Eqp_0 = (Efd - LadIfd) / Tdop; //TBD: call Exciter for Efd and LadIfd (Done)
    dx4Psidp_0 = (-x4Psidp_0 - (Xdp - Xl) * Id + x3Eqp_0) / Tdopp;
    dx5Psiqpp_0 = (-x5Psiqpp_0 - (Xq - Xdpp) * Iq) / Tqopp;

    x1d_1 = x1d_0 + dx1d_0 * t_inc;
    x2w_1 = x2w_0 + dx2w_0 * t_inc;
    x3Eqp_1 = x3Eqp_0 + dx3Eqp_0 * t_inc;
    x4Psidp_1 = x4Psidp_0 + dx4Psidp_0 * t_inc;
    x5Psiqpp_1 = x5Psiqpp_0 + dx5Psiqpp_0 * t_inc;
    p_exciter->setVterminal(presentMag);
    p_exciter->setVcomp(presentMag); //TBD update to Vcomp 
    p_exciter->setFieldCurrent(LadIfd);


    p_exciter->predictor(t_inc, flag);

    p_governor->setRotorSpeedDeviation(x2w_0);
    p_governor->predictor(t_inc, flag);
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

}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GensalGenerator::corrector(
    double t_inc, bool flag)
{
  if (getGenStatus()){
    p_exciter = getExciter();
    Efd = p_exciter->getFieldVoltage(); 

    p_governor = getGovernor();
    Pmech = p_governor->getMechanicalPower();

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
    dx2w_1 = 1 / (2 * H) * ((Pmech - D * x2w_1) / (1 + x2w_1) - Telec); //TBD: call Governor for Pmech
    dx3Eqp_1 = (Efd - LadIfd) / Tdop; //TBD: call Exciter for Efd and LadIfd
    dx4Psidp_1 = (-x4Psidp_1 - (Xdp - Xl) * Id + x3Eqp_1) / Tdopp;
    dx5Psiqpp_1 = (-x5Psiqpp_1 - (Xq - Xdpp) * Iq) / Tqopp;

    x1d_1 = x1d_0 + (dx1d_0 + dx1d_1) / 2.0 * t_inc;
    x2w_1 = x2w_0 + (dx2w_0 + dx2w_1) / 2.0 * t_inc;
    x3Eqp_1 = x3Eqp_0 + (dx3Eqp_0 + dx3Eqp_1) / 2.0 * t_inc;
    x4Psidp_1 = x4Psidp_0 + (dx4Psidp_0 + dx4Psidp_1) / 2.0 * t_inc;
    x5Psiqpp_1 = x5Psiqpp_0 + (dx5Psiqpp_0 + dx5Psiqpp_1) / 2.0 * t_inc;

    // p_exciter->setOmega(x2w_1);
    //
    p_exciter->setVterminal(presentMag);
    p_exciter->setVcomp(presentMag); 
    //  p_exciter->setFieldVoltage(Efd);
    p_exciter->setFieldCurrent(LadIfd);
    p_exciter->corrector(t_inc, flag);

    p_governor->setRotorSpeedDeviation(x2w_0);
    p_governor->corrector(t_inc, flag);

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
      sprintf(string,",%8d, %2s, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f,",
          p_bus_id, p_ckt.c_str(), x1d_1, x2w_1+1, x3Eqp_1, x4Psidp_1, x5Psiqpp_1, Vterm);
      return true;
    }
  } else if (!strcmp(signal,"debug_initial")) {
    return false;
  }
}
