
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   genrou.cpp
 * 
 * @Modified: November 21, 2022, Disable saturation if S10 and S12 are 0.
 *
 * @Modified: November 27, 2022, Fixed the model to validate against PSSE
 *
 * @Modified: Dec 9, 2022, print voltage and generator power
 *
 * @Modified: Mar 2026, Yousu Chen
 * - Fixed saturation: unscaled quadratic Se=B*(x-A)^2, Sat at Psi_ag,
 *   iterative q-axis saturation in init/predictor/corrector.
 *
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

    Vstab = 0.0;
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

  // Note: MBASE from RAW file is the per-unit base for DYR parameters.
  // Do not auto-adjust MBASE when |Sgen| > MBASE; this is valid operation.

  // GENROU machine parameter auto-correction (matches PowerWorld rules)
  // 1. Transient reactances
  if (Xdp > Xd) {
    printf("GENROU bus %d: Auto-correct Xd'=%.6f > Xd=%.6f -> Xd'=%.6f\n",
           p_bus_id, Xdp, Xd, 0.8*Xd);
    Xdp = 0.8 * Xd;
  }
  if (Xqp > Xq) {
    printf("GENROU bus %d: Auto-correct Xq'=%.6f > Xq=%.6f -> Xq'=%.6f\n",
           p_bus_id, Xqp, Xq, Xq);
    Xqp = Xq;
  }
  // 2. Subtransient reactances (Xdpp = Xqpp for GENROU)
  double Xmin = (Xdp < Xqp) ? Xdp : Xqp;
  if (Xdpp > Xmin) {
    printf("GENROU bus %d: Auto-correct Xd''=%.6f > min(Xd',Xq')=%.6f -> Xd''=%.6f\n",
           p_bus_id, Xdpp, Xmin, 0.8*Xmin);
    Xdpp = 0.8 * Xmin;
    Xqpp = Xdpp;
  }
  if (Xdpp < 0.05) {
    printf("GENROU bus %d: Auto-correct Xd''=%.6f < 0.05 -> Xd''=0.05\n",
           p_bus_id, Xdpp);
    Xdpp = 0.05;
    Xqpp = Xdpp;
  }
  // 3. Leakage reactance
  if (Xl > Xdpp) {
    printf("GENROU bus %d: Auto-correct Xl=%.6f > Xd''=%.6f -> Xl=%.6f\n",
           p_bus_id, Xl, Xdpp, 0.8*Xdpp);
    Xl = 0.8 * Xdpp;
  }
  // 4. Inertia
  if (H < 0.1) {
    printf("GENROU bus %d: Auto-correct H=%.6f < 0.1 -> H=0.1\n",
           p_bus_id, H);
    H = 0.1;
  }

}

/**
 * Update parameters in DataCollection object with current values from
 * generator
 * @param data collection object for bus that hosts generator
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::GenrouGenerator::updateData(
    boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  if (!data->setValue(GENERATOR_PG_CURRENT, genP*MBase/p_sbase, idx)) {
    data->addValue(GENERATOR_PG_CURRENT, genP*MBase/p_sbase, idx);
  }
  if (!data->setValue(GENERATOR_QG_CURRENT, genQ*MBase/p_sbase, idx)) {
    data->addValue(GENERATOR_QG_CURRENT, genQ*MBase/p_sbase, idx);
  }
  if (p_exciter.get() != NULL) {
    p_exciter->updateData(data, idx);
  }
  if (p_governor.get() != NULL) {
    p_governor->updateData(data, idx);
  }
}

/**
 * Saturation function
 * @ param x
 */
double gridpack::dynamic_simulation::GenrouGenerator::Sat(double x)
{
  if (enableSat && x > 1e-6) {
    // PowerWorld standard scaled saturation: Se(x) = B*(x-A)^2 / x
    // Fitted from: Se(1.0)=S10, Se(1.2)=S12:
    //   B*(1-A)^2/1.0 = S10,  B*(1.2-A)^2/1.2 = S12
    // => R = 1.2*S12/S10 = (1.2-A)^2/(1-A)^2
    double R = 1.2 * S12 / S10;
    double sqrtR = sqrt(R);
    double A = (1.2 - sqrtR) / (1.0 - sqrtR);
    double B = S10 / ((1.0 - A) * (1.0 - A));

    double tmp = x - A;
    if (tmp < 0.0) {
      tmp = 0.0;
    }
    double result = B * tmp * tmp / x;

    return result;
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

  // Iterative angle initialization with q-axis saturation.
  // At steady state, iron saturation reduces the effective q-axis
  // synchronous reactance.  Saturation reduces the mutual inductance
  // but not the leakage, so:  Xq_eff = Xl + (Xq - Xl) / (1 + Sq)
  // where Sq = Se(Psi_ag) * (Xq - Xl) / (Xd - Xl).
  // We iterate until the angle converges.
  double Xq_eff = Xq;
  double Vd, Vq;
  double Psiqpp, Psidpp, Psi_ag;

  for (int sat_iter = 0; sat_iter < 10; sat_iter++) {
    // Machine angle from voltage behind (Ra + j*Xq_eff)
    x1d = atan2(Viterm + Ir * Xq_eff + Ii * Ra, Vrterm + Ir * Ra - Ii * Xq_eff);

    double theta = pi/2.0 - x1d;
    rotate(Ir, Ii, theta, &Id, &Iq);

    // Norton internal voltage (behind Ra + jXdpp) in network frame
    double Vr = Vrterm + Ra * Ir - Xdpp * Ii;
    double Vi = Viterm + Ra * Ii + Xdpp * Ir;

    // Rotate to machine dq reference frame
    Vd = Vr * sin(x1d) - Vi * cos(x1d);
    Vq = Vr * cos(x1d) + Vi * sin(x1d);

    Psiqpp = -Vd;
    Psidpp = +Vq;
    Psi_ag = sqrt(Psidpp * Psidpp + Psiqpp * Psiqpp);

    double Se = Sat(Psi_ag);
    if (Se < 1e-10) break;  // no saturation, angle is exact

    double Sq = Se * (Xq - Xl) / (Xd - Xl);
    double Xq_eff_new = Xl + (Xq - Xl) / (1.0 + Sq);

    if (fabs(Xq_eff_new - Xq_eff) < 1e-8) break;  // converged
    Xq_eff = Xq_eff_new;
  }

  // q-axis transient voltage
  x3Eqp = Vq + (Xdp - Xdpp)*Id;

  // d-axis transient voltage (consistent with saturated angle)
  x6Edp = Vd - (Xqp - Xqpp)*Iq;

  // d-axis flux
  x4Psidp = x3Eqp - (Xdp - Xl)*Id;

  // q-axis flux
  x5Psiqp = x6Edp + (Xqp - Xl)*Iq;

  // Field voltage — evaluate saturation at air-gap flux magnitude
  // Standard: LadIfd = E'q + Se(Psi_ag)*Psidpp + (Xd-Xdp)*Id  (TempD=0 at SS)
  Efd = x3Eqp + Sat(Psi_ag) * Psidpp + Id * (Xd - Xdp);

  // Field current
  LadIfd = Efd;

  // Mechanical power
  Pmech = Psidpp * Iq - Psiqpp * Id;

  Efdinit = Efd;
  Pmechinit = Pmech;

  // printf("print: yuan debug here, inside genrou model, Ir=%f, Ii=%f\n", Ir, Ii);
  // Initialize exciter
  if (p_hasExciter) {
    p_exciter = getExciter();
    p_exciter->setVterminal(Vterm);
    p_exciter->setVcomp(mag); 
    p_exciter->setFieldVoltage(Efd);
    p_exciter->setFieldCurrent(LadIfd);
    //---yuan add below 20231024---//
    p_exciter->setIri(Ir, Ii); 
    //---yuan add above 20231024---//
    p_exciter->init(mag, ang, ts);
  }

  // Initialize governor
  if (p_hasGovernor) {
    p_governor = getGovernor();
    p_governor->setMechanicalPower(Pmech);
    p_governor->setRotorSpeedDeviation(x2w); // set Speed Deviation w for wsieg1
    p_governor->init(mag, ang, ts);
  }

  if (p_hasPss) {
    boost::shared_ptr<BasePssModel> pss = getPss();
    pss->setOmega(x2w);
    pss->init(mag, ang, ts);
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
 * Rebalance equilibrium after network initialization solve.
 *
 * At t=0 the Norton-augmented Y-bus produces a bus voltage that differs
 * from the power-flow voltage.  The generator states (x4,x5,x6) were
 * initialized at the PF voltage and are NOT at equilibrium for the
 * Norton voltage.  This function solves for the exact equilibrium
 * Id,Iq at the Norton voltage (a 2x2 linear system), then sets ALL
 * damper winding states (x4,x5,x6), Pmech, and Efd to their
 * equilibrium values so that every state derivative is zero.
 *
 * Derivation:  At steady state the GENROU ODEs give
 *   x4 = x3 - (Xdp-Xl)*Id,   x5 = x6 + (Xqp-Xl)*Iq,
 *   x6 = (Xq-Xqp)*Iq,        TempD = TempQ = 0.
 * Substituting into the subtransient flux expressions yields
 *   Psidpp = x3 - (Xdp-Xdpp)*Id,   Psiqpp = -(Xq-Xdpp)*Iq,
 * whence the internal voltage is Vd = (Xq-Xdpp)*Iq, Vq = x3-(Xdp-Xdpp)*Id.
 * Combined with the stator algebraic equations this gives a 2x2 linear
 * system in (Id, Iq).
 */
void gridpack::dynamic_simulation::GenrouGenerator::rebalanceEquilibrium()
{
  if (!getGenStatus()) return;

  // Norton admittance components (machine base)
  B = -Xdpp / (Ra * Ra + Xdpp * Xdpp);
  G = Ra / (Ra * Ra + Xdpp * Xdpp);

  // Terminal voltage from network solve (already set by setVolt)
  Vterm = presentMag;
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  double Vdterm = Vrterm * sin(x1d) - Viterm * cos(x1d);
  double Vqterm = Vrterm * cos(x1d) + Viterm * sin(x1d);

  // --- Solve for equilibrium Id, Iq at the Norton voltage ---
  // At equilibrium, the internal voltages become:
  //   Vd = gamma*Iq,  Vq = x3 - alpha*Id
  // where alpha = Xdp - Xdpp, gamma = Xq_eff - Xdpp.
  // With q-axis saturation: Xq_eff = Xl + (Xq-Xl)/(1+Sq),
  // so gamma depends on saturation which depends on Id,Iq.
  // We iterate until gamma converges.
  double alpha = Xdp - Xdpp;
  double gamma = Xq  - Xdpp;  // initial unsaturated value
  double Psidpp, Psiqpp, Psi_ag;

  for (int sat_iter = 0; sat_iter < 10; sat_iter++) {
    double a11 = 1.0 - alpha * B;
    double a12 = -gamma * G;
    double a21 = alpha * G;
    double a22 = 1.0 - gamma * B;

    double rhs1 = -Vdterm * G - (x3Eqp - Vqterm) * B;
    double rhs2 = -Vdterm * B + (x3Eqp - Vqterm) * G;

    double det = a11 * a22 - a12 * a21;
    Id = (rhs1 * a22 - rhs2 * a12) / det;
    Iq = (a11 * rhs2 - a21 * rhs1) / det;

    Psidpp = x3Eqp - (Xdp - Xdpp) * Id;
    Psiqpp = -gamma * Iq;
    Psi_ag = sqrt(Psidpp * Psidpp + Psiqpp * Psiqpp);

    double Se = Sat(Psi_ag);
    if (Se < 1e-10) break;

    double Sq = Se * (Xq - Xl) / (Xd - Xl);
    double gamma_new = Xl + (Xq - Xl) / (1.0 + Sq) - Xdpp;

    if (fabs(gamma_new - gamma) < 1e-8) break;
    gamma = gamma_new;
  }

  // --- Update flux states to equilibrium ---
  x4Psidp = x3Eqp - (Xdp - Xl) * Id;
  // E'd from circuit (Vd = gamma*Iq, x6Edp = Vd - (Xqp-Xqpp)*Iq)
  x6Edp   = gamma * Iq - (Xqp - Xqpp) * Iq;
  x5Psiqp = x6Edp + (Xqp - Xl) * Iq;

  // Electrical torque and field current (TempD = 0 at equilibrium)
  double Telec = Psidpp * Iq - Psiqpp * Id;
  LadIfd = x3Eqp + Sat(Psi_ag) * Psidpp + (Xd - Xdp) * Id;

  // Set Pmech = Telec so dx2w = 0
  Pmech = Telec;
  Pmechinit = Pmech;

  // Set Efd = LadIfd so dx3Eqp = 0
  Efd = LadIfd;
  Efdinit = Efd;

  // Update _1 copies to match
  x4Psidp_1 = x4Psidp;
  x5Psiqp_1 = x5Psiqp;
  x6Edp_1   = x6Edp;

  // Update network currents for watch output
  Ir = +Id * sin(x1d) + Iq * cos(x1d);
  Ii = -Id * cos(x1d) + Iq * sin(x1d);
  genP = Vrterm * Ir + Viterm * Ii;
  genQ = Viterm * Ir - Vrterm * Ii;

  // Re-initialize exciter with updated Efd
  if (p_hasExciter) {
    p_exciter = getExciter();
    p_exciter->setVterminal(Vterm);
    p_exciter->setVcomp(Vterm);
    p_exciter->setFieldVoltage(Efd);
    p_exciter->setFieldCurrent(LadIfd);
    p_exciter->init(Vterm, Theta, 0.0);
  }

  // Re-initialize governor with updated Pmech
  if (p_hasGovernor) {
    p_governor = getGovernor();
    p_governor->setMechanicalPower(Pmech);
    p_governor->setRotorSpeedDeviation(x2w);
    p_governor->init(Vterm, Theta, 0.0);
  }
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

  double Vd = - Psiqpp;
  double Vq = + Psidpp;

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

  if (getGenStatus()) {
    
    if (p_hasExciter) {
      p_exciter = getExciter();
      p_exciter->setOmega(x2w);
      p_exciter->setVterminal(presentMag);
      p_exciter->setVcomp(presentMag);

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

    double Psi_ag = sqrt(Psidpp * Psidpp + Psiqpp * Psiqpp);
    double Se = Sat(Psi_ag);
    // Standard: LadIfd = E'q + Se*Psidpp + (Xd-Xdp)*(Id+TempD)
    LadIfd = x3Eqp + Se * Psidpp + (Xd - Xdp) * (Id + TempD);

    dx1d = x2w * 2 * pi * 60; // 60 represents the nominal frequency of 60 Hz
    dx2w = 1 / (2 * H) * ((Pmech - D * x2w) / (1 + x2w) - Telec);
    dx3Eqp = (Efd - LadIfd) / Tdop;
    dx4Psidp = (-x4Psidp - (Xdp - Xl) * Id + x3Eqp) / Tdopp;
    dx5Psiqp = (-x5Psiqp + (Xqp - Xl) * Iq + x6Edp) / Tqopp;
    double TempQ = (Xqp - Xqpp) / ((Xqp - Xl) * (Xqp - Xl))
      * (-x5Psiqp + (Xqp - Xl) * Iq + x6Edp);
    // Standard: q-axis saturation with gamma_qd = (Xq-Xl)/(Xd-Xl)
    // Note: our Psiqpp = -psi_q'' (NREL), so sign is +
    dx6Edp = (-x6Edp + (Xq - Xqp) * (Iq - TempQ)
      + Se * Psiqpp * (Xq - Xl) / (Xd - Xl)) / Tqop;

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

    // PSS: run before exciter predictor so Vstab is available
    if (p_hasPss) {
      boost::shared_ptr<BasePssModel> pss = getPss();
      pss->setOmega(x2w_1);
      pss->predictor(t_inc, flag);
      Vstab = pss->getVstab();
    } else {
      Vstab = 0.0;
    }

    if (p_hasExciter) {
      if (p_hasPss) {
        p_exciter->setVstab(Vstab);
      }
      p_exciter->setFieldCurrent(LadIfd);
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
  } else {
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
  
  if (getGenStatus()) {	  
    if (p_tripped){
      p_INorton = p_Norton_Ya*vt_complex_tmp;
    } else {
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
  if (getGenStatus()) {
    if (p_hasExciter) {
      p_exciter = getExciter();
      p_exciter->setOmega(x2w_1);
      p_exciter->setVterminal(presentMag);
      p_exciter->setVcomp(presentMag);

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

    double Psi_ag = sqrt(Psidpp * Psidpp + Psiqpp * Psiqpp);
    double Se = Sat(Psi_ag);
    // Standard: LadIfd = E'q + Se*Psidpp + (Xd-Xdp)*(Id+TempD)
    LadIfd = x3Eqp_1 + Se * Psidpp + (Xd - Xdp) * (Id + TempD);
    dx1d_1 = x2w_1 * 2 * pi * 60; // 60 represents the nominal frequency of 60 Hz
    dx2w_1 = 1 / (2 * H) * ((Pmech - D * x2w_1) / (1 + x2w_1) - Telec);
    dx3Eqp_1 = (Efd - LadIfd) / Tdop;
    dx4Psidp_1 = (-x4Psidp_1 - (Xdp - Xl) * Id + x3Eqp_1) / Tdopp;
    dx5Psiqp_1 = (-x5Psiqp_1 + (Xqp - Xl) * Iq + x6Edp_1) / Tqopp;
    double TempQ = (Xqp - Xqpp) / ((Xqp - Xl) * (Xqp - Xl))
      * (-x5Psiqp_1 + (Xqp - Xl) * Iq + x6Edp_1);
    // Standard: q-axis saturation with gamma_qd = (Xq-Xl)/(Xd-Xl)
    dx6Edp_1 = (-x6Edp_1 + (Xq - Xqp) * (Iq - TempQ)
      + Se * Psiqpp * (Xq - Xl) / (Xd - Xl)) / Tqop;

    x1d = x1d + (dx1d + dx1d_1) / 2.0 * t_inc;
    x2w = x2w + (dx2w + dx2w_1) / 2.0 * t_inc;
    x3Eqp = x3Eqp + (dx3Eqp + dx3Eqp_1) / 2.0 * t_inc;
    x4Psidp = x4Psidp + (dx4Psidp + dx4Psidp_1) / 2.0 * t_inc;
    x5Psiqp = x5Psiqp + (dx5Psiqp + dx5Psiqp_1) / 2.0 * t_inc;
    x6Edp = x6Edp + (dx6Edp + dx6Edp_1) / 2.0 * t_inc;

    // PSS: run before exciter corrector so Vstab is available
    if (p_hasPss) {
      boost::shared_ptr<BasePssModel> pss = getPss();
      pss->setOmega(x2w_1);
      pss->corrector(t_inc, flag);
      Vstab = pss->getVstab();
    } else {
      Vstab = 0.0;
    }

    if (p_hasExciter) {
      if (p_hasPss) {
        p_exciter->setVstab(Vstab);
      }
      p_exciter->setFieldCurrent(LadIfd);
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
    
  } else {
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

void gridpack::dynamic_simulation::GenrouGenerator::setWideAreaFreqforPSS(double freq)
{
  p_wideareafreq = freq;
  if (p_hasPss) {
    boost::shared_ptr<BasePssModel> pss = getPss();
    pss->setWideAreaFreqforPSS(freq);
  }
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
bool gridpack::dynamic_simulation::GenrouGenerator::serialWrite(
    char* string, const int bufsize, const char *signal)
{
  bool ret = false;
  string[0] = '\0';
  if (!strcmp(signal,"standard")) {
    sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f"
        "    %12.6f	%12.6f  %12.6f\n", p_bus_id, p_ckt.c_str(), x1d_1,
        x2w_1+1.0, x3Eqp_1, x4Psidp_1, x5Psiqp_1, x6Edp_1);
    ret = true;
  } else if (!strcmp(signal,"init_debug")) {
    sprintf(string," %8d  %2s Something\n",p_bus_id,p_ckt.c_str());
    ret = true;
  } else if (!strcmp(signal,"debug_initial")) {
    sprintf(string,"Bus: %d PG: %f QG: %f\n",p_bus_id,p_pg,p_qg);
    ret = true;
  } else if(!strcmp(signal,"watch_header")) {
    if(getWatch()) {
      char buf[256];
      std::string tag;
      if(p_ckt[1] != ' ') {
        tag = p_ckt;
      } else {
        tag = p_ckt[0];
      }
      sprintf(buf,", %d_%s_V, %d_%s_Pg, %d_%s_Qg,%d_%s_angle, %d_%s_speed,"
          " %d_%s_Efd, %d_%s_Pm, %d_%s_PowerAngle",p_bus_id,tag.c_str(),p_bus_id,tag.c_str(),
          p_bus_id,tag.c_str(),p_bus_id,tag.c_str(),p_bus_id,tag.c_str(),
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
      double powerAngle = (x1d_1 - presentAng) * 180.0 / 3.14159265358979323846;
      sprintf(buf,",%f,%f,%f,%f, %f, %f, %f, %f",Vterm,genP*MBase/p_sbase,
          genQ*MBase/p_sbase,x1d_1, x2w_1+1.0, Efd,Pmech,powerAngle);
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
  if (p_generatorObservationPowerSystemBase) {
    vals.push_back(genP*MBase/p_sbase);  //output at system mva base
    vals.push_back(genQ*MBase/p_sbase);  //output at system mva base
  } else {
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
