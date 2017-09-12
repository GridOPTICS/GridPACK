/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   motorw.cpp
 * @author Shuangshuang Jin 
 * @Last modified:  Dec 12, 2016 
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_load_model.hpp"
#include "motorw.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::MotorwLoad::MotorwLoad(void)
{
  loadFactor = 0.8; // load factor for calculating MVA base using motor power
  Pul = 0.0; // percentage of load represented by this motor model
  
  rs = 0.01;
  lls = 0.1;
  lm = 3.0;
  lmr = 0.0;
  rr1 = 0.005;
  llr1 = 0.08;
  rr2 = 0.15;
  llr2 = 0.175;
  
  H = 0.4 ; // s
  MVABase = 30 ; // Motor mva base , MVA
  A =0.0; // Mechanical torque equation coefficient 
  B =0.0; // Mechanical torque equation coefficient 
  C0=0.0; // Mechanical torque equation coefficient 
  D =0.0; // Mechanical torque equation coefficient
  E =0.0; // Mechanical torque equation coefficient
  
  Ls = 0.0 ;
  lmp = 0.0;
  Lp = 0.0;
  lmpp = 0.0 ;
  Lpp =0.0 ;
  tpo = 0.0 ;
  tppo = 0.0 ;
  
  // boundary variables
  volt = 0.0 ; // pu
  freq = 0.0 ;  // pu
  Id = 0.0 ; // d-axis line current
  Iq = 0.0 ; // q-axis line current
  
  // state variables for predictor
  epq0 = 0.0 ;
  epd0 = 0.0 ;
  eppq0 = 0.0 ;
  eppd0 = 0.0 ;
  slip0 = 0.0 ;
  
  // state variables for corrector
  epq = 0.0 ;
  epd = 0.0 ;
  eppq = 0.0 ;
  eppd = 0.0 ;
  slip = 0.0 ;
  
  // derivatives of state variables for predictor
  depq_dt0 = 0.0;
  depd_dt0 = 0.0 ;
  deppq_dt0 = 0.0 ;
  deppd_dt0 = 0.0 ;
  dslip_dt0 = 0.0 ;
  
  // derivatives of state variables for corrector
  depq_dt = 0.0;
  depd_dt = 0.0 ;
  deppq_dt = 0.0 ;
  deppd_dt = 0.0 ;
  dslip_dt = 0.0 ;
  
  // other variables
  double pi = 4.0*atan(1.0);
  w0 = 2.0 * pi * 60.0 ; // rad/s

  TL = 0.0 ; // mechanical load torque, pu
  Tm0 =0.0;  // mechanical load torque with slip = 0;
  p = -1.0 ; // Assume internal voltage source to be Es = p*eppd + j q*eppq, p = 1 or -1,         
  q = -1.0 ; // q = 1 or -1
  Pmotor = 0.0 ;  // real power consumed by motor, MW
  Qmotor = 0.0 ; // reactive power consumed by motor, MVAr
  Qmotor_init = 0.0;
  sysMVABase = 100.0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::MotorwLoad::~MotorwLoad(void)
{
}

/**
 * Load parameters from DataCollection object into load model
 * @param data collection of load parameters from input files
 * @param index of load on bus
 */
void gridpack::dynamic_simulation::MotorwLoad::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx, double dloadP, double dloadQ, int ibCMPL)
{
  p_sbase = 100.0;

  data->getValue(BUS_NUMBER,&p_bus_id);
  //data->getValue(LOAD_ID,&load_id,idx);
  //if (!data->getValue(LOAD_PL, &p_pl, idx)) p_pl = 0.0;
  //if (!data->getValue(LOAD_QL, &p_ql, idx)) p_ql = 0.0;
  p_pl = dloadP;
  p_ql = dloadQ;
  
  if (ibCMPL == 1){ // if the parameters are defined from composite load model
  
    //if (!data->getValue(LOAD_MVA,  &MVABase)) 
	MVABase = 0.0;  
	if (!data->getValue(LOAD_LFM,  &loadFactor, idx)) loadFactor = 0.8;
	if (!data->getValue(LOAD_RS,   &rs, idx)) rs = 0.08;
	if (!data->getValue(LOAD_LS,   &Ls, idx)) Ls = 0.15;
	if (!data->getValue(LOAD_LP,   &Lp, idx)) Lp = 0.175;
	if (!data->getValue(LOAD_LPP,  &Lpp, idx))  Lpp = 0.08;
	if (!data->getValue(LOAD_TPO,  &tpo, idx))  tpo = 0.15;
	if (!data->getValue(LOAD_TPPO, &tppo, idx)) tppo = 0.175;
	if (!data->getValue(LOAD_H,   &H, idx)) H = 0.08;
	if (!data->getValue(LOAD_ETRQ, &E, idx)) E = 0.175;
	D = 1.0;
	A = 0.0;
	B = 0.0;
	C0 = 0.0;
	
	setDynLoadP(p_pl);
	data->getValue(LOAD_ID,&p_loadid,idx);
	setDynLoadID(p_loadid);
	
	printf ("MotorwLoad::load(), motorw with composite load model, ID: %s, loadFactor: %f, rs: %f, Ls: %f, Lp: %f, Lpp: %f, tpo: %f, tppo: %f, H: %f, A: %f, B: %f, C0: %f, D: %f, E: %f \n", 
	    p_loadid.c_str(), loadFactor, rs, Ls, Lp, Lpp, tpo, tppo, H, A, B, C0, D, E);

  } else { // if the model is CIM6BL
  
    if (!data->getValue(LOAD_PMULT, &Pul, idx)) Pul = 0.0;
    if (!data->getValue(LOAD_RA, &rs, idx)) rs = 0.08;
	if (!data->getValue(LOAD_XA, &lls, idx)) lls = 0.0;
	if (!data->getValue(LOAD_XM, &lm, idx)) lm = 0.0;
	
	if (!data->getValue(LOAD_R1, &rr1, idx)) rr1 = 0.0;
	if (!data->getValue(LOAD_X1, &llr1, idx)) llr1 = 0.0;
	if (!data->getValue(LOAD_R2, &rr2, idx)) rr2 = 0.0;
	if (!data->getValue(LOAD_X2, &llr2, idx)) llr2 = 0.0;
	if (!data->getValue(LOAD_MBASE, &MVABase, idx)) MVABase = 30.0;
	if (!data->getValue(LOAD_H, &H, idx)) H = 0.0;
	if (!data->getValue(LOAD_A, &A, idx)) A = 0.0;
	if (!data->getValue(LOAD_B, &B, idx)) B = 0.0;
	if (!data->getValue(LOAD_C0, &C0, idx)) C0 = 0.0;
	if (!data->getValue(LOAD_D, &D, idx)) D = 0.0;
	if (!data->getValue(LOAD_E, &E, idx)) E = 0.0;
	
	//tmp code please remove this
	llr1 = 0.08;
	
	data->getValue(LOAD_ID,&p_loadid,idx);
	setDynLoadP(p_pl);
    setDynLoadID(p_loadid);
	
	printf (" MotorwLoad::load(), motorw with CIM6BL, ID: %s, Pul: %f, rs: %f, lls: %f, lm: %f, rr1: %f, llr1: %f, rr2: %f, llr2: %f, MVABase: %f, \n", p_loadid.c_str(), Pul, rs, lls, lm, rr1, llr1, rr2, llr2, MVABase);
	printf (" MotorwLoad::load(), motorw with CIM6BL, ID: %s, H: %f, A: %f, B: %f, C0: %f, D: %f, E: %f, \n", p_loadid.c_str(), H, A, B, C0, D, E);
  }
}

/**
 * Initialize load model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 */
void gridpack::dynamic_simulation::MotorwLoad::init(double mag,
    double ang, double ts)
{
  //SJin: What are parameters systemMVABase, Pini, vt, and wt? From where to get them?
  //Fake declaration; 
  double systemMVABase, Pini, wt;
  gridpack::ComplexType vt(mag*cos(ang), mag*sin(ang));
  
  double pi = 4.0*atan(1.0);

  Pini = p_pl;
  wt = 2.0*60.0*pi;
  w0 = wt;
  systemMVABase = 100.0;
  sysMVABase = systemMVABase;
  
  printf("MotorwLoad::init(), vt: %12.6f +j*%12.6f, Pini: %12.6f, w0: %12.6f\n", real(vt), imag(vt), Pini, w0);
  
  // initialize the paramters
  // if the data is input in the form of motor equivalent circuit
  // paramters
  if (Ls <= 0.0 && lls > 0.0 && lm > 0.0) { 
    Ls = lls + lm ;
    lmp = lm * llr1 / (lm + llr1) ;
    Lp = lmp + lls ;
    lmpp = 1.0 / (1.0/lm + 1.0/llr1 + 1.0/llr2) ;
    Lpp = lmpp + lls ;
    tpo = llr1 * lm / (w0 * rr1 * lmp) ;
    tppo = llr2 * lmp / (w0 * rr2 * lmpp) ;
	
	printf(" MotorwLoad::init(), Ls: %12.6f, lmp: %12.6f, Lp: %12.6f, lmpp: %12.6f, Lpp: %12.6f, tpo: %12.6f, tppo: %12.6f, \n", Ls, lmp, Lp, lmpp, Lpp, tpo, tppo);
  // else if the data is input in the form of
  // subtransient/transient reactance and time constant paramters
  // need to obtain the corresponding equivalent circuit paramters
  } else if (Ls > 0.0 && Lp > 0.0) {
    if (lls ==0.0) {  
      lls = 0.1; // default paramter
    }
    // make sure the motorW is 
    if(Lpp>0.0 && lls > Lpp) {
      lls =0.8*Lpp; 
    }
    lm = Ls-lls;

    // the first cage - r1
    double xm1 = Lp - lls;
    llr1 = 1/(1/xm1-1/lm);
    rr1 = (llr1+lm)/w0/tpo;

    // the second cage - r2
    if (tppo > 0.0 && Lpp < Lp && Lpp > 0.0) {
      double xm2 = Lpp - lls;
      llr2 = 1/(1/xm2 - 1/lm - 1/llr1);
      rr2 = (llr2+xm1)/w0/tppo;
    }
  }

  // check if MVABase was properly specified, if not, set it based
  // on load factor = 0.8 
  if (MVABase <= 0.0) { 
    if (Pul>0.1) loadFactor = 1.0/Pul;
    MVABase= Pini/loadFactor;
  }

  double Vd0 = real(vt); // initial d-axis terminal voltage
  double Vq0 = imag(vt);
  //Vs0 = Vd0 + 1j * Vq0;  // pu
  gridpack::ComplexType Vs0(Vd0, Vq0);
  
  printf("    MotorwLoad::init(), Vs0: %12.6f + j*%12.6f, MVABase: %12.6f, loadFactor: %12.6f, \n", real(Vs0), imag(Vs0), MVABase, loadFactor);

  double Pe[1001];  // electrical power, MW
  double sl[1001]; // slip, pu
  gridpack::ComplexType cur[1001];  // terminal current for different slip, pu

  // sweep over the range of slip (0.0 - 1.0), calculate Pe-slip
  // relation
  sl[0] = 0.01 + 0.0011 * (1-100);
  for (int k = 1; k < 1000; k++) {
    sl[k] = 0.01 + 0.0011 * (k-100);
    if ( k < 100 ) {
      sl[k] = 0.0001*k;
    }
    //z1 = rr1/sl[k] + 1j*llr1;
    gridpack::ComplexType z1(rr1/sl[k], llr1);
    //z2 = rr2/sl[k] + 1j*llr2;
    gridpack::ComplexType z2(rr2/sl[k], llr2);
    //zr = 1j*lmr + ( z1 * z2 / (z1 + z2) );
    gridpack::ComplexType lmr_img(0.0, lmr);
    gridpack::ComplexType zr = lmr_img + (z1 * z2 / (z1+z2));
    //zs = (rs + 1j*lls) + ( 1j*lm * zr / ( zr + 1j*lm ) );
    gridpack::ComplexType a(rs, lls);
    gridpack::ComplexType b(0.0, lm);
    gridpack::ComplexType c = zr + b;
    gridpack::ComplexType zs = a + (b * zr / c);
    cur[k] = Vs0/zs;
    Pe[k] = real( Vs0 * conj(cur[k]) ) * MVABase;  // MW
	//printf("    MotorwLoad::init(), Pe[%d]: %12.6f, Sl[%d]: %12.6f, \n", k, Pe[k], k, sl[k]);
  }

  // find specific slip for initial power (Pint)
  double sl0;
  for (int m = 1; m < 1000; m++) {
    if (Pe[m] <= Pini && Pe[m+1] > Pini) {
      sl0 = sl[m] ; // found slip0
      break;
    }
  }

  slip = sl0;  // initial slip

  // calculate coefficients of 4 linear equations, state variable
  // slip is an input.
  double A1 = -1.0 ;
  double B1 = - wt * slip * tpo ;
  double C1 = (Ls - Lp) * q * Lpp / (rs*rs + Lpp*Lpp) ;
  double D1 = (Ls - Lp) * p * rs / (rs*rs + Lpp*Lpp) ;
  double E1 = -(Ls - Lp) * (Vd0 * rs + Vq0 * Lpp) / (rs*rs + Lpp*Lpp) ;

  double A2 = wt * slip * tpo ;
  double B2 = -1 ;
  double C2 = -(Ls - Lp) * q * rs / (rs*rs + Lpp*Lpp) ;
  double D2 = (Ls - Lp) * p * Lpp / (rs*rs + Lpp*Lpp) ;
  double E2 = (Ls - Lp) * (Vq0 * rs - Vd0 * Lpp) / (rs*rs + Lpp*Lpp) ;

  double A3 = 1.0 / tppo ;
  double B3 = wt * slip ;
  double C3 = (Lp - Lpp) / tppo * q * Lpp / (rs*rs + Lpp*Lpp) - 1.0 / tppo ;
  double D3 = (Lp - Lpp) / tppo * p * rs / (rs*rs + Lpp*Lpp) - wt * slip ;
  double E3 = -(Lp - Lpp) / tppo * (Vd0 * rs + Vq0 * Lpp) / (rs*rs + Lpp*Lpp) ;

  double A4 = -wt * slip ;
  double B4 = 1.0 / tppo ;
  double C4 = wt * slip - (Lp - Lpp) / tppo * q * rs / (rs*rs + Lpp*Lpp) ;
  double D4 = (Lp - Lpp) / tppo * p * Lpp / (rs*rs + Lpp*Lpp) - 1.0 / tppo ;
  double E4 = (Lp - Lpp) / tppo * (Vq0 * rs - Vd0 * Lpp) / (rs*rs + Lpp*Lpp) ;
  
  printf(" MotorwLoad::init(), A1: %12.6f, B1: %12.6f, C1: %12.6f, D1: %12.6f, E1: %12.6f, \n", A1, B1, C1, D1, E1);
  printf(" MotorwLoad::init(), A2: %12.6f, B2: %12.6f, C2: %12.6f, D2: %12.6f, E2: %12.6f, \n", A2, B2, C2, D2, E2);
  printf(" MotorwLoad::init(), A3: %12.6f, B3: %12.6f, C3: %12.6f, D3: %12.6f, E3: %12.6f, \n", A3, B3, C3, D3, E3);
  printf(" MotorwLoad::init(), A4: %12.6f, B4: %12.6f, C4: %12.6f, D4: %12.6f, E4: %12.6f, \n", A4, B4, C4, D4, E4);

  // solve the 4 linear equations to obtain 4 state variables
  epq = (B1*C2*D3*E4 - B1*C2*D4*E3 - B1*C3*D2*E4 + B1*C3*D4*E2 + B1*C4*D2*E3
  - B1*C4*D3*E2 - B2*C1*D3*E4 + B2*C1*D4*E3 + B2*C3*D1*E4 - B2*C3*D4*E1
  - B2*C4*D1*E3 + B2*C4*D3*E1 + B3*C1*D2*E4 - B3*C1*D4*E2 - B3*C2*D1*E4
  + B3*C2*D4*E1 + B3*C4*D1*E2 - B3*C4*D2*E1 - B4*C1*D2*E3 + B4*C1*D3*E2
  + B4*C2*D1*E3 - B4*C2*D3*E1 - B4*C3*D1*E2 + B4*C3*D2*E1)/
  (A1*B2*C3*D4 - A1*B2*C4*D3 - A1*B3*C2*D4 + A1*B3*C4*D2 + A1*B4*C2*D3 - 
  A1*B4*C3*D2 - A2*B1*C3*D4 + A2*B1*C4*D3 + A2*B3*C1*D4 - A2*B3*C4*D1 - 
  A2*B4*C1*D3 + A2*B4*C3*D1 + A3*B1*C2*D4 - A3*B1*C4*D2 - A3*B2*C1*D4 + 
  A3*B2*C4*D1 + A3*B4*C1*D2 - A3*B4*C2*D1 - A4*B1*C2*D3 + A4*B1*C3*D2 + 
  A4*B2*C1*D3 - A4*B2*C3*D1 - A4*B3*C1*D2 + A4*B3*C2*D1) ;

  epd =  -(A1*C2*D3*E4 - A1*C2*D4*E3 - A1*C3*D2*E4 + A1*C3*D4*E2 + A1*C4*D2*E3
  - A1*C4*D3*E2 - A2*C1*D3*E4 + A2*C1*D4*E3 + A2*C3*D1*E4 - A2*C3*D4*E1
  - A2*C4*D1*E3 + A2*C4*D3*E1 + A3*C1*D2*E4 - A3*C1*D4*E2 - A3*C2*D1*E4
  + A3*C2*D4*E1 + A3*C4*D1*E2 - A3*C4*D2*E1 - A4*C1*D2*E3 + A4*C1*D3*E2
  + A4*C2*D1*E3 - A4*C2*D3*E1 - A4*C3*D1*E2 + A4*C3*D2*E1)/
  (A1*B2*C3*D4 - A1*B2*C4*D3 - A1*B3*C2*D4 + A1*B3*C4*D2 + A1*B4*C2*D3 -
  A1*B4*C3*D2 - A2*B1*C3*D4 + A2*B1*C4*D3 + A2*B3*C1*D4 - A2*B3*C4*D1 -
  A2*B4*C1*D3 + A2*B4*C3*D1 + A3*B1*C2*D4 - A3*B1*C4*D2 - A3*B2*C1*D4 +
  A3*B2*C4*D1 + A3*B4*C1*D2 - A3*B4*C2*D1 - A4*B1*C2*D3 + A4*B1*C3*D2 +
  A4*B2*C1*D3 - A4*B2*C3*D1 - A4*B3*C1*D2 + A4*B3*C2*D1) ;

  eppq = (A1*B2*D3*E4 - A1*B2*D4*E3 - A1*B3*D2*E4 + A1*B3*D4*E2 + A1*B4*D2*E3
  - A1*B4*D3*E2 - A2*B1*D3*E4 + A2*B1*D4*E3 + A2*B3*D1*E4 - A2*B3*D4*E1
  - A2*B4*D1*E3 + A2*B4*D3*E1 + A3*B1*D2*E4 - A3*B1*D4*E2 - A3*B2*D1*E4
  + A3*B2*D4*E1 + A3*B4*D1*E2 - A3*B4*D2*E1 - A4*B1*D2*E3 + A4*B1*D3*E2
  + A4*B2*D1*E3 - A4*B2*D3*E1 - A4*B3*D1*E2 + A4*B3*D2*E1)/
  (A1*B2*C3*D4 - A1*B2*C4*D3 - A1*B3*C2*D4 + A1*B3*C4*D2 + A1*B4*C2*D3 -
  A1*B4*C3*D2 - A2*B1*C3*D4 + A2*B1*C4*D3 + A2*B3*C1*D4 - A2*B3*C4*D1 -
  A2*B4*C1*D3 + A2*B4*C3*D1 + A3*B1*C2*D4 - A3*B1*C4*D2 - A3*B2*C1*D4 +
  A3*B2*C4*D1 + A3*B4*C1*D2 - A3*B4*C2*D1 - A4*B1*C2*D3 + A4*B1*C3*D2 +
  A4*B2*C1*D3 - A4*B2*C3*D1 - A4*B3*C1*D2 + A4*B3*C2*D1) ;

  eppd = -(A1*B2*C3*E4 - A1*B2*C4*E3 - A1*B3*C2*E4 + A1*B3*C4*E2 + A1*B4*C2*E3
  - A1*B4*C3*E2 - A2*B1*C3*E4 + A2*B1*C4*E3 + A2*B3*C1*E4 - A2*B3*C4*E1
  - A2*B4*C1*E3 + A2*B4*C3*E1 + A3*B1*C2*E4 - A3*B1*C4*E2 - A3*B2*C1*E4
  + A3*B2*C4*E1 + A3*B4*C1*E2 - A3*B4*C2*E1 - A4*B1*C2*E3 + A4*B1*C3*E2
  + A4*B2*C1*E3 - A4*B2*C3*E1 - A4*B3*C1*E2 + A4*B3*C2*E1)/
  (A1*B2*C3*D4 - A1*B2*C4*D3 - A1*B3*C2*D4 + A1*B3*C4*D2 + A1*B4*C2*D3 -
  A1*B4*C3*D2 - A2*B1*C3*D4 + A2*B1*C4*D3 + A2*B3*C1*D4 - A2*B3*C4*D1 -
  A2*B4*C1*D3 + A2*B4*C3*D1 + A3*B1*C2*D4 - A3*B1*C4*D2 - A3*B2*C1*D4 +
  A3*B2*C4*D1 + A3*B4*C1*D2 - A3*B4*C2*D1 - A4*B1*C2*D3 + A4*B1*C3*D2 +
  A4*B2*C1*D3 - A4*B2*C3*D1 - A4*B3*C1*D2 + A4*B3*C2*D1) ;

  //             s = [epq, epd, eppq, eppd, slip] ;

  Id = ( (Vd0 - p * eppd) * rs + (Vq0 - q * eppq) * Lpp ) / (rs*rs + Lpp*Lpp) ;
  Iq = ( (Vq0 - q * eppq) * rs - (Vd0 - p * eppd) * Lpp ) / (rs*rs + Lpp*Lpp) ;
  TL = p * eppd * Id + q * eppq * Iq ;
  
  printf("MotorwLoad::init(), states: epq: %f, epd: %f, eppq: %f, eppd: %f, slip: %f, Id: %f, Iq: %f, TL: %f, \n", epq, epd, eppq, eppd, slip, Id, Iq, TL);

  double w = 1.0 - slip ; // rotor speed, pu
  C0 = 1.0 - A*w*w - B*w - D*(pow(w, E));
  Tm0 = TL;

  gridpack::ComplexType tmp(Id, Iq);
  //Pmotor = real( vt * conj(Id + 1j * Iq) ) * MVABase;
  //Qmotor = imag( vt * conj(Id + 1j * Iq) ) * MVABase;
  Pmotor = real( vt * conj(tmp) ) * MVABase;
  Qmotor = imag( vt * conj(tmp) ) * MVABase;
  
  printf("MotorwLoad::init(), C0: %f, Tm0: %f, Pmotor: %f, Qmotor: %f, \n", C0, Tm0, Pmotor, Qmotor);

  // slightly adjust slip to accurately match resulting real power from the
  // power flow solution
  double error = 0.0001 ;  // MW
  int flag = 0 ; // do not adjust
  if (Pini - Pmotor > 0.0) {
    flag = 1 ;  // increment slip
  } else if ((Pini - Pmotor) < 0.0) {
    flag = 2 ; // decrement slip
  }

  if (flag == 1) {
    while ((Pini - Pmotor) > error) {
      slip = slip + sl0 * error;  // increment initial slip

      // calculate coefficients of 4 linear equations, state variable
      // slip is an input.
      A1 = -1.0 ;
      B1 = - wt * slip * tpo ;
      C1 = (Ls - Lp) * q * Lpp / (rs*rs + Lpp*Lpp) ;
      D1 = (Ls - Lp) * p * rs / (rs*rs + Lpp*Lpp) ;
      E1 = -(Ls - Lp) * (Vd0 * rs + Vq0 * Lpp) / (rs*rs + Lpp*Lpp) ;

      A2 = wt * slip * tpo ;
      B2 = -1.0 ;
      C2 = -(Ls - Lp) * q * rs / (rs*rs + Lpp*Lpp) ;
      D2 = (Ls - Lp) * p * Lpp / (rs*rs + Lpp*Lpp) ;
      E2 = (Ls - Lp) * (Vq0 * rs - Vd0 * Lpp) / (rs*rs + Lpp*Lpp) ;

      A3 = 1.0 / tppo ;
      B3 = wt * slip ;
      C3 = (Lp - Lpp) / tppo * q * Lpp / (rs*rs + Lpp*Lpp) - 1.0 / tppo ;
      D3 = (Lp - Lpp) / tppo * p * rs / (rs*rs + Lpp*Lpp) - wt * slip ;
      E3 = -(Lp - Lpp) / tppo * (Vd0 * rs + Vq0 * Lpp) / (rs*rs + Lpp*Lpp) ;

      A4 = -wt * slip ;
      B4 = 1.0 / tppo ;
      C4 = wt * slip - (Lp - Lpp) / tppo * q * rs / (rs*rs + Lpp*Lpp) ;
      D4 = (Lp - Lpp) / tppo * p * Lpp / (rs*rs + Lpp*Lpp) - 1.0 / tppo ;
      E4 = (Lp - Lpp) / tppo * (Vq0 * rs - Vd0 * Lpp) / (rs*rs + Lpp*Lpp) ;

      // solve the 4 linear equations to obtain 4 state variables
      epq = (B1*C2*D3*E4 - B1*C2*D4*E3 - B1*C3*D2*E4 + B1*C3*D4*E2 + B1*C4*D2*E3
	  - B1*C4*D3*E2 - B2*C1*D3*E4 + B2*C1*D4*E3 + B2*C3*D1*E4 - B2*C3*D4*E1
	  - B2*C4*D1*E3 + B2*C4*D3*E1 + B3*C1*D2*E4 - B3*C1*D4*E2 - B3*C2*D1*E4
	  + B3*C2*D4*E1 + B3*C4*D1*E2 - B3*C4*D2*E1 - B4*C1*D2*E3 + B4*C1*D3*E2
	  + B4*C2*D1*E3 - B4*C2*D3*E1 - B4*C3*D1*E2 + B4*C3*D2*E1)/
	  (A1*B2*C3*D4 - A1*B2*C4*D3 - A1*B3*C2*D4 + A1*B3*C4*D2 + A1*B4*C2*D3 - 
	  A1*B4*C3*D2 - A2*B1*C3*D4 + A2*B1*C4*D3 + A2*B3*C1*D4 - A2*B3*C4*D1 - 
	  A2*B4*C1*D3 + A2*B4*C3*D1 + A3*B1*C2*D4 - A3*B1*C4*D2 - A3*B2*C1*D4 + 
	  A3*B2*C4*D1 + A3*B4*C1*D2 - A3*B4*C2*D1 - A4*B1*C2*D3 + A4*B1*C3*D2 + 
	  A4*B2*C1*D3 - A4*B2*C3*D1 - A4*B3*C1*D2 + A4*B3*C2*D1) ;

      epd =  -(A1*C2*D3*E4 - A1*C2*D4*E3 - A1*C3*D2*E4 + A1*C3*D4*E2 + A1*C4*D2*E3
	  - A1*C4*D3*E2 - A2*C1*D3*E4 + A2*C1*D4*E3 + A2*C3*D1*E4 - A2*C3*D4*E1
	  - A2*C4*D1*E3 + A2*C4*D3*E1 + A3*C1*D2*E4 - A3*C1*D4*E2 - A3*C2*D1*E4
	  + A3*C2*D4*E1 + A3*C4*D1*E2 - A3*C4*D2*E1 - A4*C1*D2*E3 + A4*C1*D3*E2
	  + A4*C2*D1*E3 - A4*C2*D3*E1 - A4*C3*D1*E2 + A4*C3*D2*E1)/
	  (A1*B2*C3*D4 - A1*B2*C4*D3 - A1*B3*C2*D4 + A1*B3*C4*D2 + A1*B4*C2*D3 -
	  A1*B4*C3*D2 - A2*B1*C3*D4 + A2*B1*C4*D3 + A2*B3*C1*D4 - A2*B3*C4*D1 -
	  A2*B4*C1*D3 + A2*B4*C3*D1 + A3*B1*C2*D4 - A3*B1*C4*D2 - A3*B2*C1*D4 +
	  A3*B2*C4*D1 + A3*B4*C1*D2 - A3*B4*C2*D1 - A4*B1*C2*D3 + A4*B1*C3*D2 +
	  A4*B2*C1*D3 - A4*B2*C3*D1 - A4*B3*C1*D2 + A4*B3*C2*D1) ;

      eppq = (A1*B2*D3*E4 - A1*B2*D4*E3 - A1*B3*D2*E4 + A1*B3*D4*E2 + A1*B4*D2*E3
	  - A1*B4*D3*E2 - A2*B1*D3*E4 + A2*B1*D4*E3 + A2*B3*D1*E4 - A2*B3*D4*E1
	  - A2*B4*D1*E3 + A2*B4*D3*E1 + A3*B1*D2*E4 - A3*B1*D4*E2 - A3*B2*D1*E4
	  + A3*B2*D4*E1 + A3*B4*D1*E2 - A3*B4*D2*E1 - A4*B1*D2*E3 + A4*B1*D3*E2
	  + A4*B2*D1*E3 - A4*B2*D3*E1 - A4*B3*D1*E2 + A4*B3*D2*E1)/
	  (A1*B2*C3*D4 - A1*B2*C4*D3 - A1*B3*C2*D4 + A1*B3*C4*D2 + A1*B4*C2*D3 -
	  A1*B4*C3*D2 - A2*B1*C3*D4 + A2*B1*C4*D3 + A2*B3*C1*D4 - A2*B3*C4*D1 -
	  A2*B4*C1*D3 + A2*B4*C3*D1 + A3*B1*C2*D4 - A3*B1*C4*D2 - A3*B2*C1*D4 +
	  A3*B2*C4*D1 + A3*B4*C1*D2 - A3*B4*C2*D1 - A4*B1*C2*D3 + A4*B1*C3*D2 +
	  A4*B2*C1*D3 - A4*B2*C3*D1 - A4*B3*C1*D2 + A4*B3*C2*D1) ;

      eppd = -(A1*B2*C3*E4 - A1*B2*C4*E3 - A1*B3*C2*E4 + A1*B3*C4*E2 + A1*B4*C2*E3
	  - A1*B4*C3*E2 - A2*B1*C3*E4 + A2*B1*C4*E3 + A2*B3*C1*E4 - A2*B3*C4*E1
	  - A2*B4*C1*E3 + A2*B4*C3*E1 + A3*B1*C2*E4 - A3*B1*C4*E2 - A3*B2*C1*E4
	  + A3*B2*C4*E1 + A3*B4*C1*E2 - A3*B4*C2*E1 - A4*B1*C2*E3 + A4*B1*C3*E2
	  + A4*B2*C1*E3 - A4*B2*C3*E1 - A4*B3*C1*E2 + A4*B3*C2*E1)/
	  (A1*B2*C3*D4 - A1*B2*C4*D3 - A1*B3*C2*D4 + A1*B3*C4*D2 + A1*B4*C2*D3 -
	  A1*B4*C3*D2 - A2*B1*C3*D4 + A2*B1*C4*D3 + A2*B3*C1*D4 - A2*B3*C4*D1 -
	  A2*B4*C1*D3 + A2*B4*C3*D1 + A3*B1*C2*D4 - A3*B1*C4*D2 - A3*B2*C1*D4 +
	  A3*B2*C4*D1 + A3*B4*C1*D2 - A3*B4*C2*D1 - A4*B1*C2*D3 + A4*B1*C3*D2 +
	  A4*B2*C1*D3 - A4*B2*C3*D1 - A4*B3*C1*D2 + A4*B3*C2*D1) ;

      Id = ( (Vd0 - p * eppd) * rs + (Vq0 - q * eppq) * Lpp ) / (rs*rs + Lpp*Lpp) ;
      Iq = ( (Vq0 - q * eppq) * rs - (Vd0 - p * eppd) * Lpp ) / (rs*rs + Lpp*Lpp) ;
      TL = p * eppd * Id + q * eppq * Iq ;

      w = 1.0 - slip ; // rotor speed, pu
      C0 = 1.0 - A*w*w - B*w - D*(pow(w,E));
      Tm0 = TL;

      gridpack::ComplexType tmp2(Id, Iq);
      //Pmotor = real( vt * conj(Id + 1j * Iq) ) * MVABase;
      //Qmotor = imag( vt * conj(Id + 1j * Iq) ) * MVABase;
      Pmotor = real( vt * conj(tmp2) ) * MVABase;
      Qmotor = imag( vt * conj(tmp2) ) * MVABase;
    }

  } else if (flag == 2) {
    while ((Pini - Pmotor) < -error) {
      slip = slip - sl0 * error;  // increment initial slip

      // calculate coefficients of 4 linear equations, state variable
      // slip is an input.
      A1 = -1.0 ;
      B1 = - wt * slip * tpo ;
      C1 = (Ls - Lp) * q * Lpp / (rs*rs + Lpp*Lpp) ;
      D1 = (Ls - Lp) * p * rs / (rs*rs + Lpp*Lpp) ;
      E1 = -(Ls - Lp) * (Vd0 * rs + Vq0 * Lpp) / (rs*rs + Lpp*Lpp) ;

      A2 = wt * slip * tpo ;
      B2 = -1.0 ;
      C2 = -(Ls - Lp) * q * rs / (rs*rs + Lpp*Lpp) ;
      D2 = (Ls - Lp) * p * Lpp / (rs*rs + Lpp*Lpp) ;
      E2 = (Ls - Lp) * (Vq0 * rs - Vd0 * Lpp) / (rs*rs + Lpp*Lpp) ;

      A3 = 1.0 / tppo ;
      B3 = wt * slip ;
      C3 = (Lp - Lpp) / tppo * q * Lpp / (rs*rs + Lpp*Lpp) - 1.0 / tppo ;
      D3 = (Lp - Lpp) / tppo * p * rs / (rs*rs + Lpp*Lpp) - wt * slip ;
      E3 = -(Lp - Lpp) / tppo * (Vd0 * rs + Vq0 * Lpp) / (rs*rs + Lpp*Lpp) ;

      A4 = -wt * slip ;
      B4 = 1.0 / tppo ;
      C4 = wt * slip - (Lp - Lpp) / tppo * q * rs / (rs*rs + Lpp*Lpp) ;
      D4 = (Lp - Lpp) / tppo * p * Lpp / (rs*rs + Lpp*Lpp) - 1.0 / tppo ;
      E4 = (Lp - Lpp) / tppo * (Vq0 * rs - Vd0 * Lpp) / (rs*rs + Lpp*Lpp) ;

      // solve the 4 linear equations to obtain 4 state variables
      epq = (B1*C2*D3*E4 - B1*C2*D4*E3 - B1*C3*D2*E4 + B1*C3*D4*E2 + B1*C4*D2*E3
	  - B1*C4*D3*E2 - B2*C1*D3*E4 + B2*C1*D4*E3 + B2*C3*D1*E4 - B2*C3*D4*E1
	  - B2*C4*D1*E3 + B2*C4*D3*E1 + B3*C1*D2*E4 - B3*C1*D4*E2 - B3*C2*D1*E4
	  + B3*C2*D4*E1 + B3*C4*D1*E2 - B3*C4*D2*E1 - B4*C1*D2*E3 + B4*C1*D3*E2
	  + B4*C2*D1*E3 - B4*C2*D3*E1 - B4*C3*D1*E2 + B4*C3*D2*E1)/
	  (A1*B2*C3*D4 - A1*B2*C4*D3 - A1*B3*C2*D4 + A1*B3*C4*D2 + A1*B4*C2*D3 - 
	  A1*B4*C3*D2 - A2*B1*C3*D4 + A2*B1*C4*D3 + A2*B3*C1*D4 - A2*B3*C4*D1 - 
	  A2*B4*C1*D3 + A2*B4*C3*D1 + A3*B1*C2*D4 - A3*B1*C4*D2 - A3*B2*C1*D4 + 
	  A3*B2*C4*D1 + A3*B4*C1*D2 - A3*B4*C2*D1 - A4*B1*C2*D3 + A4*B1*C3*D2 + 
	  A4*B2*C1*D3 - A4*B2*C3*D1 - A4*B3*C1*D2 + A4*B3*C2*D1) ;

      epd =  -(A1*C2*D3*E4 - A1*C2*D4*E3 - A1*C3*D2*E4 + A1*C3*D4*E2 + A1*C4*D2*E3
	  - A1*C4*D3*E2 - A2*C1*D3*E4 + A2*C1*D4*E3 + A2*C3*D1*E4 - A2*C3*D4*E1
	  - A2*C4*D1*E3 + A2*C4*D3*E1 + A3*C1*D2*E4 - A3*C1*D4*E2 - A3*C2*D1*E4
	  + A3*C2*D4*E1 + A3*C4*D1*E2 - A3*C4*D2*E1 - A4*C1*D2*E3 + A4*C1*D3*E2
	  + A4*C2*D1*E3 - A4*C2*D3*E1 - A4*C3*D1*E2 + A4*C3*D2*E1)/
	  (A1*B2*C3*D4 - A1*B2*C4*D3 - A1*B3*C2*D4 + A1*B3*C4*D2 + A1*B4*C2*D3 -
	  A1*B4*C3*D2 - A2*B1*C3*D4 + A2*B1*C4*D3 + A2*B3*C1*D4 - A2*B3*C4*D1 -
	  A2*B4*C1*D3 + A2*B4*C3*D1 + A3*B1*C2*D4 - A3*B1*C4*D2 - A3*B2*C1*D4 +
	  A3*B2*C4*D1 + A3*B4*C1*D2 - A3*B4*C2*D1 - A4*B1*C2*D3 + A4*B1*C3*D2 +
	  A4*B2*C1*D3 - A4*B2*C3*D1 - A4*B3*C1*D2 + A4*B3*C2*D1) ;

      eppq = (A1*B2*D3*E4 - A1*B2*D4*E3 - A1*B3*D2*E4 + A1*B3*D4*E2 + A1*B4*D2*E3
	  - A1*B4*D3*E2 - A2*B1*D3*E4 + A2*B1*D4*E3 + A2*B3*D1*E4 - A2*B3*D4*E1
	  - A2*B4*D1*E3 + A2*B4*D3*E1 + A3*B1*D2*E4 - A3*B1*D4*E2 - A3*B2*D1*E4
	  + A3*B2*D4*E1 + A3*B4*D1*E2 - A3*B4*D2*E1 - A4*B1*D2*E3 + A4*B1*D3*E2
	  + A4*B2*D1*E3 - A4*B2*D3*E1 - A4*B3*D1*E2 + A4*B3*D2*E1)/
	  (A1*B2*C3*D4 - A1*B2*C4*D3 - A1*B3*C2*D4 + A1*B3*C4*D2 + A1*B4*C2*D3 -
	  A1*B4*C3*D2 - A2*B1*C3*D4 + A2*B1*C4*D3 + A2*B3*C1*D4 - A2*B3*C4*D1 -
	  A2*B4*C1*D3 + A2*B4*C3*D1 + A3*B1*C2*D4 - A3*B1*C4*D2 - A3*B2*C1*D4 +
	  A3*B2*C4*D1 + A3*B4*C1*D2 - A3*B4*C2*D1 - A4*B1*C2*D3 + A4*B1*C3*D2 +
	  A4*B2*C1*D3 - A4*B2*C3*D1 - A4*B3*C1*D2 + A4*B3*C2*D1) ;

      eppd = -(A1*B2*C3*E4 - A1*B2*C4*E3 - A1*B3*C2*E4 + A1*B3*C4*E2 + A1*B4*C2*E3
	  - A1*B4*C3*E2 - A2*B1*C3*E4 + A2*B1*C4*E3 + A2*B3*C1*E4 - A2*B3*C4*E1
	  - A2*B4*C1*E3 + A2*B4*C3*E1 + A3*B1*C2*E4 - A3*B1*C4*E2 - A3*B2*C1*E4
	  + A3*B2*C4*E1 + A3*B4*C1*E2 - A3*B4*C2*E1 - A4*B1*C2*E3 + A4*B1*C3*E2
	  + A4*B2*C1*E3 - A4*B2*C3*E1 - A4*B3*C1*E2 + A4*B3*C2*E1)/
	  (A1*B2*C3*D4 - A1*B2*C4*D3 - A1*B3*C2*D4 + A1*B3*C4*D2 + A1*B4*C2*D3 -
	  A1*B4*C3*D2 - A2*B1*C3*D4 + A2*B1*C4*D3 + A2*B3*C1*D4 - A2*B3*C4*D1 -
	  A2*B4*C1*D3 + A2*B4*C3*D1 + A3*B1*C2*D4 - A3*B1*C4*D2 - A3*B2*C1*D4 +
	  A3*B2*C4*D1 + A3*B4*C1*D2 - A3*B4*C2*D1 - A4*B1*C2*D3 + A4*B1*C3*D2 +
	  A4*B2*C1*D3 - A4*B2*C3*D1 - A4*B3*C1*D2 + A4*B3*C2*D1) ;

      //             s = [epq, epd, eppq, eppd, slip] ;

      Id = ( (Vd0 - p * eppd) * rs + (Vq0 - q * eppq) * Lpp ) / (rs*rs + Lpp*Lpp) ;
      Iq = ( (Vq0 - q * eppq) * rs - (Vd0 - p * eppd) * Lpp ) / (rs*rs + Lpp*Lpp) ;
      TL = p * eppd * Id + q * eppq * Iq ;

      w = 1.0 - slip ; // rotor speed, pu
      C0 = 1,0 - A*w*w - B*w - D*(pow(w,E));
      Tm0 = TL;

      gridpack::ComplexType tmp2(Id, Iq);
      //Pmotor = real( vt * conj(Id + 1j * Iq) ) * MVABase;
      //Qmotor = imag( vt * conj(Id + 1j * Iq) ) * MVABase;
      Pmotor = real( vt * conj(tmp2) ) * MVABase;
      Qmotor = imag( vt * conj(tmp2) ) * MVABase;
    }
  }

  printf(" MotorwLoad::init(), states after slightly adjust, epq: %f, epd: %f, eppq: %f, eppd: %f, slip: %f, Id: %f, Iq: %f, TL: %f, \n", epq, epd, eppq, eppd, slip, Id, Iq, TL);
  printf(" MotorwLoad::init(),  after slightly adjust, C0: %f, Tm0: %f, Pmotor: %f, Qmotor: %f, \n", C0, Tm0, Pmotor, Qmotor);
    
  epq0  = epq;
  epd0  = epd; 
  eppq0 = eppq;
  eppd0 = eppd;
  slip0 = slip;
  
  Qmotor_init = Qmotor; 
  setDynLoadQ(Qmotor_init);
  
}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType  gridpack::dynamic_simulation::MotorwLoad::INorton()
{
  //SJIN: Matlab getINorton?
  return p_INorton;
}


/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType  gridpack::dynamic_simulation::MotorwLoad::NortonImpedence()
{
  //SJIN: Matlab getNortonImpedance?
  //return ??; // refer to gensal.cpp
  gridpack::ComplexType  Yn(0.0, 0.0);
  gridpack::ComplexType  temp(rs, Lpp);
  Yn = 1.0 / temp;
  Yn = Yn * MVABase / sysMVABase;
  printf("MotorwLoad::NortonImpedence(), Yn: %12.6f + j*%12.6f \n", real(Yn), imag(Yn));
  return Yn;
}

/**
 * Predict part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::MotorwLoad::predictor_currentInjection(bool flag)
{
  if (!flag) {
    epq0 = epq ;
    epd0 = epd ;
    eppq0 = eppq ;
    eppd0 = eppd ;
    slip0 = slip ;
  }

  // calculate current injection for post-predictor network
  // solution, notice the direction. In this case, Inject !
  gridpack::ComplexType  In(0.0, 0.0);
  gridpack::ComplexType  a(p*eppd, q*eppq);
  gridpack::ComplexType  b(rs, Lpp);
  //In = ( p * eppd + 1j * q * eppq ) / ( rs + 1j * Lpp ) ; // pu
  In = a / b;
  In = In * MVABase / sysMVABase ;  // convert Norton injection current from motor base to system base
  p_INorton = In;
  printf("MotorwLoad::predictor_currentInjection(), p_INorton: %12.6f + j*%12.6f \n", real(In), imag(In));
} 

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::MotorwLoad::predictor(
    double t_inc, bool flag)
{
  //SJin: What are parameters vt, wt, and dt? From where to get them?
  //Fake declaration; 
  gridpack::ComplexType vt = vt_complex;
  double pi = 4.0*atan(1.0);
  double wt = presentFreq*2.0*60.0*pi;
  double dt = t_inc;
  printf("MotorwLoad::predictor(), vt: %12.6f + j*%12.6f, wt: %12.6f \n", real(vt), imag(vt), wt);

  // Step-1: update predictor state variables using corrector
  // state variables;
  // s0 = s ;
  if (!flag) {
	epq0 = epq ;
	epd0 = epd ;
	eppq0 = eppq ;
	eppd0 = eppd ;
	slip0 = slip ;
  }

  // Step-2: calculate predictor dx/dt
  //Es = p * eppd0 + 1j * q * eppq0 ;  // p * eppd + j q * eppq
  gridpack::ComplexType Es(p*eppd0, q*eppq0);
  gridpack::ComplexType tmp(rs, Lpp);
  //Id = real( ( vt - Es ) / (rs + 1j * Lpp) ) ;
  //Iq = imag( ( vt - Es ) / (rs + 1j * Lpp) ) ;
  Id = real( ( vt - Es ) / tmp ) ;
  Iq = imag( ( vt - Es ) / tmp ) ;

  double w = 1.0 - slip0 ; // rotor speed, pu

  TL=  Tm0*(A*w*w + B*w + C0+D*(pow(w,E)));

  depq_dt0 =  ( -( epq0 + (Ls - Lp)* Id ) - wt * slip0 * tpo * epd0 ) / tpo ;  // d_epq
  depd_dt0 = ( -( epd0 - (Ls - Lp)* Iq ) + wt * slip0 * tpo * epq0 ) / tpo ;  // d_epd
  deppq_dt0 = ( epq0 - eppq0 - (Lp - Lpp) * Id )/ tppo + wt * slip0 * ( epd0 - eppd0 ) + depq_dt0 ;  // d_eppq
  deppd_dt0 = ( epd0 - eppd0 + (Lp - Lpp) * Iq )/ tppo - wt * slip0 * ( epq0 - eppq0 ) + depd_dt0 ;  // d_eppd
  dslip_dt0 = -( p * eppd0 * Id + q * eppq0 * Iq - TL ) / (2.0 * H * w) ;

  // Step-3: integrate
  epq = epq0 +  depq_dt0 * dt;
  epd = epd0 + depd_dt0 * dt;
  eppq = eppq0 + deppq_dt0 * dt ;
  eppd = eppd0 + deppd_dt0 * dt ;
  slip = slip0 + dslip_dt0 * dt ;

  // Step-4: update outputs
  //Es = p * eppd + 1j * q * eppq ;  // p * eppd + j q * eppq
  //real(Es) = p*eppd;
  //imag(Es) = q*eppq;
  Es = gridpack::ComplexType(p*eppd,q*eppq);
  gridpack::ComplexType tmp2(rs, Lpp);
  //Id = real( ( vt - Es ) / (rs + 1j * Lpp) ) ;
  //Iq = imag( ( vt - Es ) / (rs + 1j * Lpp) ) ;
  Id = real( ( vt - Es ) / tmp2 ) ;
  Iq = imag( ( vt - Es ) / tmp2 ) ;
  gridpack::ComplexType tmp3(Id, Iq);
  //Pmotor = real( vt * conj(Id + 1j * Iq) ) * MVABase;
  //Qmotor = imag( vt * conj(Id + 1j * Iq) ) * MVABase;
  Pmotor = real( vt * conj(tmp3) ) * MVABase;
  Qmotor = imag( vt * conj(tmp3) ) * MVABase;
  
  printf(" MotorwLoad::predictor(), %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f\n", presentMag, presentAng, presentFreq, epq, epd, eppq, eppd, 1.0-slip, Pmotor, Qmotor, Id, Iq );
    
}

/**
 * Corrector part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::MotorwLoad::corrector_currentInjection(bool flag)
{
  // calculate current injection for post-predictor network
  // solution, notice the direction. In this case, Inject !
  gridpack::ComplexType  In(0.0, 0.0);
  gridpack::ComplexType  a(p*eppd, q*eppq);
  gridpack::ComplexType  b(rs, Lpp);
  //In = ( p * eppd + 1j * q * eppq ) / ( rs + 1j * Lpp ) ; // pu
  In = a / b;
  In = In * MVABase / sysMVABase ;  // convert Norton injection current from motor base to system base
  p_INorton = In;
  printf("MotorwLoad::corrector_currentInjection(), p_INorton: %12.6f + j*%12.6f \n", real(In), imag(In));
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::MotorwLoad::corrector(
    double t_inc, bool flag)
{
  //SJin: What are parameters vt, wt, and dt? From where to get them?
  //Fake declaration; 
  gridpack::ComplexType vt = vt_complex;
  double pi = 4.0*atan(1.0);
  double wt = presentFreq*2.0*60.0*pi;
  double dt = t_inc;

  printf("MotorwLoad::corrector(), vt: %12.6f + j*%12.6f, wt: %12.6f \n", real(vt), imag(vt), wt);
  
  //g Step-1: calculate corrector dx'/dt
  //Es = p * eppd + 1j * q * eppq ;  //g p * eppd + j q * eppq
  gridpack::ComplexType Es(p*eppd, q*eppq);
  gridpack::ComplexType tmp(rs, Lpp);
  //Id = real( ( vt - Es ) / (rs + 1j * Lpp) ) ;
  //Iq = imag( ( vt - Es ) / (rs + 1j * Lpp) ) ;
  Id = real( ( vt - Es ) / tmp ) ;
  Iq = imag( ( vt - Es ) / tmp ) ;

  double w = 1.0 - slip;  //g rotor speed, pu
  TL=  Tm0*(A*w*w + B*w + C0 + D*(pow(w,E)));

  depq_dt =  ( -( epq + (Ls - Lp)* Id ) - wt * slip * tpo * epd ) / tpo ;  //g d_epq
  depd_dt = ( -( epd - (Ls - Lp)* Iq ) + wt * slip * tpo * epq ) / tpo ;  //g d_epd
  deppq_dt = ( epq - eppq - (Lp - Lpp) * Id )/ tppo + wt * slip * ( epd - eppd ) + depq_dt ;  //g d_eppq
  deppd_dt = ( epd - eppd + (Lp - Lpp) * Iq )/ tppo - wt * slip * ( epq - eppq ) + depd_dt ;  //g d_eppd
  dslip_dt = -( p * eppd * Id + q * eppq * Iq - TL ) / (2.0 * H * w) ;

  //g Step-2: integrate
  epq = epq0 +  0.5 * ( depq_dt0 + depq_dt  )* dt;
  epd = epd0 + 0.5 * ( depd_dt0 + depd_dt ) * dt;
  eppq = eppq0 + 0.5 * ( deppq_dt0 + deppq_dt ) * dt ;
  eppd = eppd0 + 0.5 * ( deppd_dt0 + deppd_dt ) * dt ;
  slip = slip0 + 0.5 * ( dslip_dt0 + dslip_dt ) * dt ;

  //g Step-3: update outputs
  //Es = p * eppd + 1j * q * eppq ;  //g p * eppd + j q * eppq
  //gridpack::ComplexType Es(p*eppd, q*eppq);
  //real(Es) = p*eppd;
  //imag(Es) = q*eppq;
  Es = gridpack::ComplexType(p*eppd,q*eppq);
  gridpack::ComplexType tmp2(rs, Lpp);
  //Id = real( ( vt - Es ) / (rs + 1j * Lpp) ) ;
  //Iq = imag( ( vt - Es ) / (rs + 1j * Lpp) ) ;
  Id = real( ( vt - Es ) / tmp2 ) ;
  Iq = imag( ( vt - Es ) / tmp2 ) ;
  gridpack::ComplexType tmp3(Id, Iq);
  //Pmotor = real( vt * conj(Id + 1j * Iq) ) * MVABase;
  //Qmotor = imag( vt * conj(Id + 1j * Iq) ) * MVABase;
  Pmotor = real( vt * conj(tmp3) ) * MVABase;
  Qmotor = imag( vt * conj(tmp3) ) * MVABase;
  
  printf(" Output MotorwLoad::corrector(), bus: %d,  %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f \n", p_bus_id, presentMag, presentAng, presentFreq, epq, epd, eppq, eppd, 1.0-slip, Pmotor, Qmotor, TL, Id, Iq );
    
}

/**
 * Set voltage on each load
 */
void gridpack::dynamic_simulation::MotorwLoad::setVoltage(
    gridpack::ComplexType voltage)
{
  printf("MotorwLoad::setVoltage, %12.6f + j %12.6f \n", real(voltage), imag(voltage));	
  presentMag = abs(voltage);
  presentAng = atan2(imag(voltage), real(voltage));  
  vt_complex = voltage;
  
}

/**
 * Set terminal voltage frequency on each load
 */
void gridpack::dynamic_simulation::MotorwLoad::setFreq(double dFreq)
{
  presentFreq = dFreq;
}

/**
 * Write out load state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
bool gridpack::dynamic_simulation::MotorwLoad::serialWrite(
    char* string, const int bufsize, const char *signal)
{
  return false;
}

double gridpack::dynamic_simulation::MotorwLoad::getMotorQ() 
{
  return Qmotor;
}
/**
     * get intialized reactive power of the dynamic load model
     */
double gridpack::dynamic_simulation::MotorwLoad::getInitReactivePower() 
{
  return Qmotor_init;
}


