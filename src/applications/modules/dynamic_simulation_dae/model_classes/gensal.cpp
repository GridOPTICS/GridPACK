/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   genrou.cpp
 * @author Zakaria El Mrabet
 * @Created:   10/15/19
 * @Modified:  12/04/19 - Shri
 *  
 * @brief  
 *
 *
 */

#include <gensal.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

GensalGen::GensalGen(void)
{
  delta = 0.0; // Rotor angle
  dw = 0.0; // Rotor speed
  Eqp = 0.0; // Transient Q axis Eq 
  Psidp = 0.0; // Transient D axis flux
  Psiqpp = 0.0; // Transient Q axis flux
  ddelta = 0.0;
  ddw = 0.0;
  dEqp = 0.0;
  dPsidp = 0.0;
  dPsiqpp = 0.0;
  Ra = 0.0; // Machine stator resistance
  H = 0.0; // Machine inertia constant
  D = 0.0; // Machine damping coefficient
  Pmech = 0.0; // Mechanical power
  Xd = 0.0;
  Xq = 0.0;
  Xdp = 0.0; // Machine transient reactance
  Xdpp = 0.0;
  Xl = 0.0;
  Tdop = 0.0;
  Tdopp = 0.0;
  Tqopp = 0.0;
  S10 = 0.0;
  S12 = 0.0;
  Tqop = 0.0;

  B = 0.0;
  G = 0.0;

  nxgen = 5;
}

GensalGen::~GensalGen(void)
{
  if(p_hasExciter) {
    free(xexc_loc);
    free(dEfd_dxexc);
    free(dEfd_dxgen);
  }

  if(p_hasGovernor) {
    free(xgov_loc);
    free(dPmech_dxgov);
  }
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseGeneratorModel
 */
void GensalGen::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseGenModel::load(data,idx); // load parameters in base generator model

  // load parameters for the model type

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
  if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.55; // S12 
  if (!data->getValue(GENERATOR_TQOP, &Tqop, idx)) Tqop=0.0; // Tqop

  // Values given on mbase. Convert to sbase
  double mult = mbase/sbase;
  H *= mult;
  D *= mult;
  Ra /= mult;
  Xd /= mult;
  Xq /= mult;
  Xdp /= mult;
  Xdpp /= mult;
  Xl   /= mult;

  // Saturation constants
  double temp = sqrt(S10/(1.2*S12));
  sat_A = (1.0 - temp*1.2)/(1-temp);
  sat_B = S10/((1.0 - sat_A)*(1.0 - sat_A));

  // Set up arrays for generator exciter Jacobian coupling
  if(p_hasExciter) {
    int nxexc;
    p_exciter = getExciter();
    p_exciter->vectorSize(&nxexc);
    xexc_loc = (int*)malloc(nxexc*sizeof(int));
    dEfd_dxexc = (double*)malloc(nxexc*sizeof(double));
    dEfd_dxgen  = (double*)malloc(nxgen*sizeof(double));
  }

  if(p_hasGovernor) {
    int nxgov;
    p_governor = getGovernor();
    p_governor->vectorSize(&nxgov);
    xgov_loc = (int*)malloc(nxgov*sizeof(int));
    dPmech_dxgov = (double*)malloc(nxgov*sizeof(double));
  }
}

/**
 * Saturation function
 * @ param x
 */
double GensalGen::Sat(double Eqp)
{
    double result = sat_B * (Eqp - sat_A) * (Eqp - sat_A) /Eqp;
    if(Eqp < sat_A) result = 0.0;
    //    return result; // Scaled Quadratic with 1.7.1 equations\
    return 0.0; // disabled saturation effect
}

/**
 * Initialize generator model before calculation
 * @param [output] values - array where initialized generator variables should be set
 */
void GensalGen::init(gridpack::ComplexType* values) 
{
  double Vterm = sqrt(VD*VD + VQ*VQ); 
  double P, Q; // Generator real and reactive power
  double Vrterm = VD;
  double Viterm = VQ;
  double Ir,Ii,Id,Iq;
  double Eppr,Eppi,Vd,Vq,Vdterm,Vqterm;
  double Psid,Psiq,Psidpp,Telec;


  // Admittance
  B = -Xdpp / (Ra * Ra + Xdpp * Xdpp);
  G = Ra / (Ra * Ra + Xdpp * Xdpp);

  P = pg / sbase;
  Q = qg / sbase;

  Ir = (P * Vrterm + Q * Viterm) / (Vterm * Vterm);
  Ii = (P * Viterm - Q * Vrterm) / (Vterm * Vterm);
  dw = 0.0;
  delta = atan2(Viterm + Ir * Xq + Ii * Ra, Vrterm + Ir * Ra - Ii * Xq);

  Eppr = Vrterm + Ra * Ir - Xdpp * Ii; // internal voltage on network reference
  Eppi = Viterm + Ra * Ii + Xdpp * Ir; // internal voltage on network reference
  Vd = Eppr * sin(delta) - Eppi * cos(delta); // convert to dq reference
  Vq = Eppr * cos(delta) + Eppi * sin(delta); // convert to dq reference
  Vdterm = VD*sin(delta) - VQ*cos(delta);
  Vqterm = VD*cos(delta) + VQ*sin(delta);
  Id = (Vd-Vdterm)*G - (Vq-Vqterm)*B;
  Iq = (Vd-Vdterm)*B + (Vq-Vqterm)*G;

  Psiqpp = (Xdpp - Xq) * Iq;

  Psiq = Psiqpp - Iq*Xdpp;
  Psidpp = Vq;
  Psid = Psidpp - Id*Xdpp;
  Telec = Psid* Iq - Psiq * Id;

  Psidp = Psidpp - Id * (Xdpp - Xl);
  Eqp   = Psidp + Id * (Xdp - Xl);

  Efd = Eqp * (1 + Sat(Eqp)) + Id * (Xd - Xdp); 
  LadIfd = Efd;
  Pmech = Psid * Iq - Psiq * Id;

  values[0] = delta;
  values[1] = dw;
  values[2] = Eqp;
  values[3] = Psidp;
  values[4] = Psiqpp;
  
  // Initialize exciters
  if (p_hasExciter) {
    p_exciter = getExciter();
    p_exciter->setInitialFieldVoltage(Efd);
  }

  if (p_hasGovernor) {
    p_governor = getGovernor();
    p_governor->setInitialMechanicalPower(Pmech);
  }

}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */


double GensalGen::getAngle(void)
{
  return delta;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void GensalGen::write(const char* signal, char* string)
{
}

/**
 *  Set the number of variables for this generator model
 *  @param [output] number of variables for this model
 */
bool GensalGen::vectorSize(int *nvar) const
{
  *nvar = nxgen;
  return true;
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto generators
 * @param values array containing generator state variables
*/
void GensalGen::setValues(gridpack::ComplexType *values)
{
  if(p_mode == XVECTOBUS) {
    delta = real(values[0]);
    dw = real(values[1]);
    Eqp = real(values[2]);
    Psidp = real(values[3]);
    Psiqpp = real(values[4]);
  } else if(p_mode == XDOTVECTOBUS) {
    ddelta = real(values[0]);
    ddw = real(values[1]);
    dEqp = real(values[2]);
    dPsidp = real(values[3]);
    dPsiqpp = real(values[4]);
  } else if(p_mode == XVECPRETOBUS) {
    deltaprev = real(values[0]);
    dwprev = real(values[1]);
    Eqpprev = real(values[2]);
    Psidpprev = real(values[3]);
    Psiqppprev = real(values[4]);
  }    
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
bool GensalGen::vectorValues(gridpack::ComplexType *values)
{
  int delta_idx = 0;
  int dw_idx = 1;
  int Eqp_idx = 2;
  int Psidp_idx = 3;
  int Psiqpp_idx = 4;
  double Vd,Vq,Vdterm,Vqterm;
  double Id,Iq;
  double Psid,Psiq,Psidpp;
  double Telec,TempD;
  // On fault (p_mode == FAULT_EVAL flag), the generator variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    values[delta_idx] = delta - deltaprev;
    values[dw_idx] = dw - dwprev;
    values[Eqp_idx] = Eqp - Eqpprev;
    values[Psidp_idx] = Psidp - Psidpprev;
    values[Psiqpp_idx] = Psiqpp - Psiqppprev;

  } else if(p_mode == RESIDUAL_EVAL) {
    if (p_hasExciter) {
      p_exciter = getExciter();
      Efd = p_exciter->getFieldVoltage(); // Efd obtained from exciter
    }
    
    if (p_hasGovernor) {
      p_governor = getGovernor();
      Pmech = p_governor->getMechanicalPower();
    }

    // Generator equations
    Psidpp = + Eqp * (Xdpp - Xl) / (Xdp - Xl) + Psidp * (Xdp - Xdpp) / (Xdp - Xl);
    Vd = -Psiqpp * (1 + dw);
    Vq = +Psidpp * (1 + dw);
    Vdterm = VD*sin(delta) - VQ*cos(delta);
    Vqterm = VD*cos(delta) + VQ*sin(delta);
    Id = (Vd-Vdterm)*G - (Vq-Vqterm)*B;
    Iq = (Vd-Vdterm)*B + (Vq-Vqterm)*G;

    Psiq = Psiqpp - Iq * Xdpp;
    Psid  = Psidpp - Id * Xdpp;
    Telec = Psid * Iq - Psiq * Id;
    TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl)) * (-Psidp - (Xdp - Xl) * Id + Eqp);
    LadIfd = Eqp * (1 + Sat(Eqp)) + (Xd - Xdp) * (Id + TempD);

    // RESIDUAL_EVAL for state 1 to 5
    values[delta_idx] = dw * OMEGA_S - ddelta;
    values[dw_idx] = 1 / (2 * H) * ((Pmech - D * dw) / (1 + dw) - Telec) - ddw; 
    values[Eqp_idx] = (Efd - LadIfd) / Tdop - dEqp; 
    values[Psidp_idx] = (-Psidp - (Xdp - Xl) * Id + Eqp) / Tdopp - dPsidp;
    values[Psiqpp_idx] = (-Psiqpp - (Xq - Xdpp) * Iq ) / Tqopp - dPsiqpp;

  }

  return true;
}

/**
 * Set Jacobian values
 * @param values a 2-d array of Jacobian block for the bus
 */
bool GensalGen::setJacobian(gridpack::ComplexType **values) 
{

  int VD_idx = 0; /* Row/col number for bus voltage VD variable */
  int VQ_idx = 1; /* Row/col number for bus voltage VQ variable */
  int IGQ_idx = 0; /* Row/col location for IGQ equations */
  int IGD_idx = 1; /* Row/col location for IGD equations */
  int delta_idx = offsetb;
  int dw_idx    = offsetb+1;
  int Eqp_idx   = offsetb+2;
  int Psidp_idx = offsetb+3;
  int Psiqpp_idx = offsetb+4;
  double Vd,Vq,Vdterm,Vqterm;
  double Id,Iq;
  double Psid,Psiq,Psidpp;
  double Telec,TempD;

  if (p_hasExciter) {
    p_exciter = getExciter();
    Efd = p_exciter->getFieldVoltage(); // Efd obtained from exciter
  }
  
  if (p_hasGovernor) {
    p_governor = getGovernor();
    Pmech = p_governor->getMechanicalPower();
  }
  
  Psidpp = + Eqp * (Xdpp - Xl) / (Xdp - Xl) + Psidp * (Xdp - Xdpp) / (Xdp - Xl);
  Vd = -Psiqpp * (1 + dw);
  Vq = +Psidpp * (1 + dw);
  Vdterm = VD*sin(delta) - VQ*cos(delta);
  Vqterm = VD*cos(delta) + VQ*sin(delta);
  Id = (Vd-Vdterm)*G - (Vq-Vqterm)*B;
  Iq = (Vd-Vdterm)*B + (Vq-Vqterm)*G;
  
  double dPsidpp_ddelta = 0.0;
  double dPsidpp_ddw    = 0.0;
  double dPsidpp_dEqp   = (Xdpp - Xl)/(Xdp - Xl);
  double dPsidpp_dPsidp = (Xdp - Xdpp)/(Xdp - Xl);
  double dPsidpp_dPsiqpp = 0.0;
  
  double dVd_dPsiqpp = -(1 + dw);
  double dVd_ddw    = -Psiqpp;
  
  double dVq_dEqp   = dPsidpp_dEqp*(1 + dw);
  double dVq_dPsidp = dPsidpp_dPsidp*(1 + dw);
  double dVq_ddw    = Psidpp;
  
  double dVdterm_dVD = sin(delta), dVdterm_dVQ = -cos(delta);
  double dVqterm_dVD = cos(delta), dVqterm_dVQ =  sin(delta);
  
  double dVdterm_ddelta =  VD*cos(delta) + VQ*sin(delta);
  double dVqterm_ddelta = -VD*sin(delta) + VQ*cos(delta);
  
  double dId_dVd     =  G, dId_dVq     = -B;
  double dId_dVdterm = -G, dId_dVqterm =  B;
  
  double dIq_dVd     =  B, dIq_dVq     =  G;
  double dIq_dVdterm = -B, dIq_dVqterm = -G;
  
  double dId_ddelta = dId_dVdterm*dVdterm_ddelta + dId_dVqterm*dVqterm_ddelta;
  double dIq_ddelta = dIq_dVdterm*dVdterm_ddelta + dIq_dVqterm*dVqterm_ddelta;
  
  double dId_ddw = dId_dVd*dVd_ddw + dId_dVq*dVq_ddw;
  double dIq_ddw = dIq_dVd*dVd_ddw + dIq_dVq*dVq_ddw;
  
  double dId_dEqp = dId_dVq*dVq_dEqp;
  double dIq_dEqp = dIq_dVq*dVq_dEqp;
  
  double dId_dPsidp = dId_dVq*dVq_dPsidp;
  double dIq_dPsidp = dIq_dVq*dVq_dPsidp;
  
  double dId_dPsiqpp = dId_dVd*dVd_dPsiqpp;
  double dIq_dPsiqpp = dIq_dVd*dVd_dPsiqpp;
  
  double dId_dVD = dId_dVdterm*dVdterm_dVD + dId_dVqterm*dVqterm_dVD;
  double dId_dVQ = dId_dVdterm*dVdterm_dVQ + dId_dVqterm*dVqterm_dVQ;
  
  double dIq_dVD = dIq_dVdterm*dVdterm_dVD + dIq_dVqterm*dVqterm_dVD;
  double dIq_dVQ = dIq_dVdterm*dVdterm_dVQ + dIq_dVqterm*dVqterm_dVQ;
  
  Psiq = Psiqpp - Iq * Xdpp;
  Psid  = Psidpp - Id * Xdpp;
  
  double dPsiq_ddelta  = -Xdpp*dIq_ddelta;
  double dPsiq_ddw     = -Xdpp*dIq_ddw;
  double dPsiq_dEqp    = -Xdpp*dIq_dEqp;
  double dPsiq_dPsidp  = -Xdpp*dIq_dPsidp;
  double dPsiq_dPsiqpp = 1.0 - Xdpp*dIq_dPsiqpp;
  
  double dPsiq_dVD     = -Xdpp*dIq_dVD;
  double dPsiq_dVQ     = -Xdpp*dIq_dVQ;
  
  double dPsid_ddelta  = dPsidpp_ddelta - Xdpp*dId_ddelta;
  double dPsid_ddw     = dPsidpp_ddw - Xdpp*dId_ddw;
  double dPsid_dEqp    = dPsidpp_dEqp - Xdpp*dId_dEqp;
  double dPsid_dPsidp  = dPsidpp_dPsidp - Xdpp*dId_dPsidp;
  double dPsid_dPsiqpp = dPsidpp_dPsiqpp - Xdpp*dId_dPsiqpp;
  
  double dPsid_dVD     =  -Xdpp*dId_dVD;
  double dPsid_dVQ     =  -Xdpp*dId_dVQ;
  
  double dTelec_ddelta = dPsid_ddelta*Iq + Psid*dIq_ddelta
    - dPsiq_ddelta*Id - Psiq*dId_ddelta;
  double dTelec_ddw    = dPsid_ddw*Iq + Psid*dIq_ddw
    - dPsiq_ddw*Id - Psiq*dId_ddw;
  double dTelec_dEqp   = dPsid_dEqp*Iq + Psid*dIq_dEqp
    - dPsiq_dEqp*Id - Psiq*dId_dEqp;
  double dTelec_dPsidp = dPsid_dPsidp*Iq + Psid*dIq_dPsidp
    - dPsiq_dPsidp*Id - Psiq*dId_dPsidp;
  double dTelec_dPsiqpp= dPsid_dPsiqpp*Iq + Psid*dIq_dPsiqpp
    - dPsiq_dPsiqpp*Id - Psiq*dId_dPsiqpp;
  double dTelec_dVD    = dPsid_dVD*Iq + Psid*dIq_dVD
    - dPsiq_dVD*Id - Psiq*dId_dVD;
  double dTelec_dVQ    = dPsid_dVQ*Iq + Psid*dIq_dVQ
    - dPsiq_dVQ*Id - Psiq*dId_dVQ;
  
  if(p_mode == FAULT_EVAL) {
    // Generator variables held constant
    // dF_dX
    // Set diagonal values to 1.0
    values[delta_idx][delta_idx] = 1.0;
    values[dw_idx][dw_idx] = 1.0;
    values[Eqp_idx][Eqp_idx] = 1.0;
    values[Psidp_idx][dw_idx] = 1.0;
    values[Psiqpp_idx][Psiqpp_idx] = 1.0;
  } else {
    
    values[delta_idx][delta_idx] = -shift;
    values[dw_idx][delta_idx]    = OMEGA_S;

    // F2  values[dw_idx] = 1 / (2 * H) * ((Pmech - D * dw) / (1 + dw) - Telec) - ddw; // Pmech can be called from Governor
    // dF2_dx
    double const1 = 1/(2*H);
    
    values[delta_idx][dw_idx] = const1*-dTelec_ddelta;
    values[dw_idx][dw_idx] = const1*(-(Pmech-D*dw)/((1+dw)*(1+dw)) - D/(1+dw) - dTelec_ddw) - shift;
    values[Eqp_idx][dw_idx] = const1*-dTelec_dEqp;
    values[Psidp_idx][dw_idx] = const1*-dTelec_dPsidp;
    values[Psiqpp_idx][dw_idx] = const1*-dTelec_dPsiqpp;
    
    values[VD_idx][dw_idx] = const1*-dTelec_dVD;
    values[VQ_idx][dw_idx] = const1*-dTelec_dVQ;
    
    // Add Pmech contributions
    if(p_hasGovernor) {
      int nxgov,i;
      p_governor->vectorSize(&nxgov);
      p_governor->getMechanicalPowerPartialDerivatives(xgov_loc,dPmech_dxgov);
      
      /* Partials w.r.t. governor mechanical power Pmech */
      for(i=0; i < nxgov; i++) {
	values[xgov_loc[i]][dw_idx] = const1*dPmech_dxgov[i]/(1+dw);
      }
    }
    
    double const2 =  (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl));
    double dLadIfd_ddelta =  (Xd - Xdp)*(dId_ddelta + const2*(-(Xdp - Xl)*dId_ddelta));
    double dLadIfd_ddw    =  (Xd - Xdp)*(dId_ddw + const2*(-(Xdp - Xl)*dId_ddw));
    double dLadIfd_dEqp   = 1.0 + (Xd - Xdp)*(dId_dEqp + const2*(1.0 - (Xdp - Xl)*dId_dEqp)); // No saturation considered yet
    double dLadIfd_dPsidp =  (Xd - Xdp)*(dId_dPsidp + const2*(-1.0 -(Xdp - Xl)*dId_dPsidp));
    double dLadIfd_dPsiqpp=  (Xd - Xdp)*(dId_dPsiqpp + const2*(-(Xdp - Xl)*dId_dPsiqpp));
    
    double dLadIfd_dVD   =  (Xd - Xdp)*(dId_dVD + const2*(-(Xdp - Xl)*dId_dVD));
    double dLadIfd_dVQ   =  (Xd - Xdp)*(dId_dVQ + const2*(-(Xdp - Xl)*dId_dVQ));
    
    values[delta_idx][Eqp_idx] = -dLadIfd_ddelta/Tdop;
    values[dw_idx][Eqp_idx]    = -dLadIfd_ddw/Tdop;
    values[Eqp_idx][Eqp_idx]   = -dLadIfd_dEqp/Tdop - shift;
    values[Psidp_idx][Eqp_idx] = -dLadIfd_dPsidp/Tdop;
    values[Psiqpp_idx][Eqp_idx]= -dLadIfd_dPsiqpp/Tdop;
    
    values[VD_idx][Eqp_idx]    = -dLadIfd_dVD/Tdop;
    values[VQ_idx][Eqp_idx]    = -dLadIfd_dVQ/Tdop;
    
    // Add Efd contributions
    if(p_hasExciter) {
      int nexc,i;
      p_exciter->vectorSize(&nexc);
      p_exciter->getFieldVoltagePartialDerivatives(xexc_loc,dEfd_dxexc,dEfd_dxgen);
      
      /* Partials w.r.t. exciter variables */
      for(i=0; i < nexc; i++) {
	values[xexc_loc[i]][Eqp_idx] = dEfd_dxexc[i]/Tdop;
      }

      /* Partials contributions for Efd w.r.t. generator variables. Note that this may be because
	 the exciter Efd calculation uses the field current LadIfd (which is a function of generator variables)
      */
      for(i=0; i < nxgen; i++) {
	values[offsetb+i][Eqp_idx] += dEfd_dxgen[i]/Tdop;
      }
    }

    // State 4
    values[delta_idx][Psidp_idx] = -(Xdp - Xl)*dId_ddelta/Tdopp;
    values[dw_idx][Psidp_idx]    = -(Xdp - Xl)*dId_ddw/Tdopp;
    values[Eqp_idx][Psidp_idx]   = (1.0 -(Xdp - Xl)*dId_dEqp)/Tdopp;
    values[Psidp_idx][Psidp_idx] = (-1.0 -(Xdp - Xl)*dId_dPsidp)/Tdopp - shift;
    values[Psiqpp_idx][Psidp_idx]= (-(Xdp - Xl)*dId_dPsiqpp)/Tdopp;

    values[VD_idx][Psidp_idx] = -(Xdp - Xl)*dId_dVD/Tdopp;
    values[VQ_idx][Psidp_idx] = -(Xdp - Xl)*dId_dVQ/Tdopp;

    // State 5
    values[delta_idx][Psiqpp_idx] = -(Xq - Xdpp)*dIq_ddelta/Tqopp;
    values[dw_idx][Psiqpp_idx]    = -(Xq - Xdpp)*dIq_ddw/Tqopp;
    values[Eqp_idx][Psiqpp_idx]   = -(Xq - Xdpp)*dIq_dEqp/Tqopp;
    values[Psidp_idx][Psiqpp_idx] = -(Xq - Xdpp)*dIq_dPsidp/Tqopp;
    values[Psiqpp_idx][Psiqpp_idx]= (-1.0 -(Xq - Xdpp)*dIq_dPsiqpp)/Tqopp - shift;
    
    values[VD_idx][Psiqpp_idx]    = -(Xq - Xdpp)*dIq_dVD/Tqopp;
    values[VQ_idx][Psiqpp_idx]    = -(Xq - Xdpp)*dIq_dVQ/Tqopp;
  }

  // Partial of generator currents w.r.t. x and V
  values[delta_idx][IGD_idx]    = dId_ddelta*sin(delta) + Id*cos(delta)
    + dIq_ddelta*cos(delta) - Iq*sin(delta);
  values[dw_idx][IGD_idx]       = dId_ddw*sin(delta) + dIq_ddw*cos(delta);
  values[Eqp_idx][IGD_idx]      = dId_dEqp*sin(delta) + dIq_dEqp*cos(delta);
  values[Psidp_idx][IGD_idx]    = dId_dPsidp*sin(delta) + dIq_dPsidp*cos(delta);
  values[Psiqpp_idx][IGD_idx]   = dId_dPsiqpp*sin(delta) + dIq_dPsiqpp*cos(delta);
  values[VD_idx][IGD_idx]      += dId_dVD*sin(delta) + dIq_dVD*cos(delta);
  values[VQ_idx][IGD_idx]      += dId_dVQ*sin(delta) + dIq_dVQ*cos(delta);
  
  values[delta_idx][IGQ_idx]    = -dId_ddelta*cos(delta) + Id*sin(delta)
    + dIq_ddelta*sin(delta) + Iq*cos(delta);
  values[dw_idx][IGQ_idx]       = -dId_ddw*cos(delta) + dIq_ddw*sin(delta);
  values[Eqp_idx][IGQ_idx]      = -dId_dEqp*cos(delta) + dIq_dEqp*sin(delta);
  values[Psidp_idx][IGQ_idx]     = -dId_dPsidp*cos(delta) + dIq_dPsidp*sin(delta);
  values[Psiqpp_idx][IGQ_idx]   = -dId_dPsiqpp*cos(delta) + dIq_dPsiqpp*sin(delta);
  values[VD_idx][IGQ_idx]      += -dId_dVD*cos(delta) + dIq_dVD*sin(delta);
  values[VQ_idx][IGQ_idx]      += -dId_dVQ*cos(delta) + dIq_dVQ*sin(delta);
      
  return true;
}


/**
 * Return the generator current injection (in rectangular form) 
 * @param [output] IGD - real part of the generator current 
 * @param [output] IGQ - imaginary part of the generator current 
*/
void GensalGen::getCurrent(double *IGD, double *IGQ)
{
  double Psidpp = + Eqp * (Xdpp - Xl) / (Xdp - Xl) + Psidp * (Xdp - Xdpp) / (Xdp - Xl); 
  double Vd = -Psiqpp * (1 + dw);
  double Vq = +Psidpp * (1 + dw);
  double Vdterm = VD*sin(delta) - VQ*cos(delta);
  double Vqterm = VD*cos(delta) + VQ*sin(delta);
  double Id = (Vd-Vdterm)*G - (Vq-Vqterm)*B;
  double Iq = (Vd-Vdterm)*B + (Vq-Vqterm)*G;

  // Generator current injections in the network
  *IGD =   Id * sin(delta) + Iq * cos(delta);
  *IGQ =  -Id * cos(delta) + Iq * sin(delta);
}

double GensalGen::getRotorSpeedDeviation()
{
  return dw;
}

int GensalGen::getRotorSpeedDeviationLocation()
{
  return offsetb+1;
}


double GensalGen::getFieldCurrent()
{
  double Psidpp = + Eqp * (Xdpp - Xl) / (Xdp - Xl) + Psidp * (Xdp - Xdpp) / (Xdp - Xl); 
  double Vd = -Psiqpp * (1 + dw);
  double Vq = +Psidpp * (1 + dw);
  double Vdterm = VD*sin(delta) - VQ*cos(delta);
  double Vqterm = VD*cos(delta) + VQ*sin(delta);
  double Id = (Vd-Vdterm)*G - (Vq-Vqterm)*B;
  double Iq = (Vd-Vdterm)*B + (Vq-Vqterm)*G;
  double Psiq = Psiqpp - Iq * Xdpp;
  double Psid  = Psidpp - Id * Xdpp;
  double TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl)) * (-Psidp - (Xdp - Xl) * Id + Eqp);
  double LadIfd = Eqp * (1 + Sat(Eqp)) + (Xd - Xdp) * (Id + TempD);

  return LadIfd;

}

/**
 * Return the matrix entries
 * @param [output] nval - number of values set
 * @param [output] row - row indics for matrix entries
 * @param [output] col - col indices for matrix entries
 * @param [output] values - matrix entries
 * return true when matrix entries set
 */

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
bool GensalGen::serialWrite(
    char* string, const int bufsize, const char *signal)
{
  return false;
}
