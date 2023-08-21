/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   genrou.cpp
 *  
 * @brief  Implements the GENROU generator model
 *
 *
 */

#include <genrou.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

GenrouGen::GenrouGen(void)
{
  delta = 0.0; // Rotor angle
  dw = 0.0; // Rotor speed
  Eqp = 0.0; // Transient Q axis Eq 
  Psidp = 0.0; // Transient D axis flux
  Psiqp = 0.0; // Transient Q axis flux
  Edp = 0.0; // Transient D axis Ed
  ddelta = 0.0;
  ddw = 0.0;
  dEqp = 0.0;
  dPsidp = 0.0;
  dPsiqp = 0.0;
  dEdp = 0.0;
  Ra = 0.0; // Machine stator resistance
  H = 0.0; // Machine inertia constant
  D = 0.0; // Machine damping coefficient
  Pmech = 0.0; // Mechanical power
  Id = 0.0; // Generator current on d axis
  Iq = 0.0; // Generator current on q axis
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
  Xqp = 0.0;
  Xqpp = 0.0;
  Tqop = 0.0;

  B = 0.0;
  G = 0.0;

  nxgen = 6; // Number of variables for this model 
}

GenrouGen::~GenrouGen(void)
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
void GenrouGen::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseGenModel::load(data,idx); // load parameters in base generator model

  // load parameters for the model type
  data->getValue(BUS_NUMBER, &busnum);
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
  if (!data->getValue(GENERATOR_S1, &S10, idx)) S10=0.067; // S10 TBD: check parser
  if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.579; // S12 TBD: check parser
  if (!data->getValue(GENERATOR_XQP, &Xqp, idx)) Xqp=0.0; // Xqp
  if (!data->getValue(GENERATOR_XDPP, &Xqpp, idx)) Xqpp=Xdpp; // Xqpp 
  if (!data->getValue(GENERATOR_TQOP, &Tqop, idx)) Tqop=0.0; // Tqop

  // Convert generator parameters from machine base to MVA base
  double mult = mbase/sbase;
  H *= mult;
  D *= mult;

  Xdp /= mult;
  Xd /= mult;
  Xq /= mult;
  Xdpp /= mult;
  Xl /= mult;
  Xqp /= mult;
  Xqpp /= mult;

  B = -Xdpp/(Ra*Ra + Xdpp*Xdpp);
  G = Ra/(Ra*Ra + Xdpp*Xdpp);

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

  // Set up arrays for generator governor coupling
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
/** Note: The saturation function is currently disabled as it needs a slightly
    complicated process for initialization. This process involves solving a set
    of nonlinear equations simultaneously to initialize the generator variables.
    This initialization is not implemented, hence the saturation function is disabled.
**/
double GenrouGen::Sat(double Psidpp,double Psiqpp)
{
  double Psipp = sqrt(Psidpp*Psidpp + Psiqpp*Psiqpp);
  double result = sat_B * (Psipp - sat_A) * (Psipp - sat_A) /Psipp;
  if(Psipp < sat_A) result = 0.0;
  //  return result;
  return 0.0; // disabled the saturation effect
}

/**
 * Initialize generator model before calculation
 * @param [output] values - array where initialized generator variables should be set
 */
void GenrouGen::init(gridpack::ComplexType* values) 
{
  double Vm = sqrt(VD*VD + VQ*VQ); 
  double theta;
  double Ir,Ii;
  double P, Q; // Generator real and reactive power
  double Vd,Vq;
  double Psid,Psiq,Psidpp,Psiqpp;
  double Eppr,Eppi,Vdterm,Vqterm;

  P = pg / sbase;
  Q = qg / sbase;

  // Generator currents in network reference frame
  Ir = (P * VD + Q * VQ) / (Vm * Vm);
  Ii = (P * VQ - Q * VD) / (Vm * Vm);

  // Machine angle and speed deviation
  delta = atan2(VQ + Ir * Xq + Ii * Ra, VD + Ir * Ra - Ii * Xq);
  dw = 0.0;

  theta = delta - PI/2.0; // Rotor angle lags quadrature axis by 90 deg.

  Eppr = VD + Ra * Ir - Xdpp * Ii; // internal voltage on network reference
  Eppi = VQ + Ra * Ii + Xdpp * Ir; // internal voltage on network reference
  Vd = Eppr * cos(theta) + Eppi * sin(theta); // convert to dq reference
  Vq = Eppr *-sin(theta) + Eppi * cos(theta); // convert to dq reference
  Vdterm = VD*cos(theta) + VQ*sin(theta);
  Vqterm = VD*-sin(theta) + VQ*cos(theta);
  Id = (Vd-Vdterm)*G - (Vq-Vqterm)*B;
  Iq = (Vd-Vdterm)*B + (Vq-Vqterm)*G;

  Psidpp = Vq;
  Psiqpp = -Vd;
  double Psipp = sqrt(Psidpp*Psidpp + Psiqpp*Psiqpp);

  Psid = Psidpp - Id*Xdpp;
  Psiq = Psiqpp - Iq*Xdpp;

  Psidp = Psidpp - (Xdpp - Xl)*Id;
  Eqp   = Psidp  + (Xdp - Xl)*Id;

  Edp = (Xq - Xqp)*Iq - (Xq - Xl)/(Xd - Xl)*Psiqpp*Sat(Psidpp,Psiqpp)/Psipp;
  Psiqp = Edp + (Xqp - Xl)*Iq;

  Efd = Eqp + (Xd - Xdp)*Id  + Psidpp*Sat(Psidpp,Psiqpp)/Psipp;;

  LadIfd = Efd;
  Pmech = P;

  double Telec = Psidpp*Iq - Psiqpp*Id;
  values[0] = delta;
  values[1] = dw;
  values[2] = Eqp;
  values[3] = Psidp;
  values[4] = Psiqp;
  values[5] = Edp;

  // Initialize exciter field voltage and current 
  if (p_hasExciter) {
    p_exciter->setInitialFieldVoltage(Efd);
  }
    
  // Initialize governor mechanical power and speed deviation
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
bool GenrouGen::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

double GenrouGen::getAngle(void)
{
  return delta;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void GenrouGen::write(const char* signal, char* string)
{

}

/**
 *  Set the number of variables for this generator model
 *  @param [output] number of variables for this model
 */
bool GenrouGen::vectorSize(int *nvar) const
{
  *nvar = nxgen;
  return true;
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto generators
 * @param values array containing generator state variables
*/
void GenrouGen::setValues(gridpack::ComplexType *values)
{
  if(p_mode == XVECTOBUS) {
    delta = real(values[0]);
    dw = real(values[1]);
    Eqp = real(values[2]);
    Psidp = real(values[3]);
    Psiqp = real(values[4]);
    Edp = real(values[5]);
  } else if(p_mode == XDOTVECTOBUS) {
    ddelta = real(values[0]);
    ddw = real(values[1]);
    dEqp = real(values[2]);
    dPsidp = real(values[3]);
    dPsiqp = real(values[4]);
    dEdp = real(values[5]);
  } else if(p_mode == XVECPRETOBUS) {
    /* This state is only called once during fault on/off condition */
    deltaprev = real(values[0]);
    dwprev = real(values[1]);
    Eqpprev = real(values[2]);
    Psidpprev = real(values[3]);
    Psiqpprev = real(values[4]);
    Edpprev = real(values[5]);
  }
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
bool GenrouGen::vectorValues(gridpack::ComplexType *values)
{
  int delta_idx = 0;
  int dw_idx = 1;
  int Eqp_idx = 2;
  int Psidp_idx = 3;
  int Psiqp_idx = 4;
  int Edp_idx = 5;
  // On fault (p_mode == FAULT_EVAL flag), the generator variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    values[delta_idx] = delta - deltaprev;
    values[dw_idx] = dw - dwprev;
    values[Eqp_idx] = Eqp - Eqpprev;
    values[Psidp_idx] = Psidp - Psidpprev;
    values[Psiqp_idx] = Psiqp - Psiqpprev;
    values[Edp_idx] = Edp - Edpprev;

  } else if(p_mode == RESIDUAL_EVAL) {

    if (p_hasExciter) {
      p_exciter = getExciter();
      Efd = p_exciter->getFieldVoltage(); // Efd are called from Exciter
    }
    
    if (p_hasGovernor) {
      p_governor = getGovernor();
      Pmech = p_governor->getMechanicalPower();
    }

    // Generator equations
    double theta = delta - PI/2.0;
    double Vdterm = VD * cos(theta) + VQ * sin(theta);
    double Vqterm = VD *-sin(theta) + VQ * cos(theta);

    double tempd1,tempd2,tempq1,tempq2;
    tempd1 = (Xdpp - Xl)/(Xdp - Xl);
    tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
    tempq1 = (Xdpp - Xl)/(Xqp - Xl);
    tempq2 = (Xqp - Xdpp)/(Xqp - Xl);

    double Psidpp =  tempd1*Eqp + tempd2*Psidp;
    double Psiqpp = -tempq1*Edp - tempq2*Psiqp;
    double Psipp  = sqrt(Psidpp*Psidpp + Psiqpp*Psiqpp);

    double Vd = -Psiqpp*(1 + dw);
    double Vq =  Psidpp*(1 + dw);

    // dq Axis currents
    Id = (Vd-Vdterm)*G - (Vq-Vqterm)*B;
    Iq = (Vd-Vdterm)*B + (Vq-Vqterm)*G;

    // Electrical torque
    double Telec = Psidpp*Iq - Psiqpp*Id;

    // Field current
    double LadIfd = getFieldCurrent();

    // RESIDUAL_EVAL for state 1 to 6
    // Machine rotor angle
    values[delta_idx] = dw * OMEGA_S - ddelta;

    // Speed
    values[dw_idx] = 1 / (2 * H) * ((Pmech - D * dw) / (1 + dw) - Telec) - ddw; // Pmech can be called from Governor
    //    values[dw_idx] = 1 / (2 * H) * ((Pmech - D * dw) - Telec) - ddw; // Ignoring speed effects

    // Q-axis transient EMF
    values[Eqp_idx] = (Efd - LadIfd) / Tdop - dEqp; 

    // D-axis transient flux
    values[Psidp_idx] = (-Psidp - (Xdp - Xl) * Id + Eqp) / Tdopp - dPsidp;
    
    // Q-axis transient flux
    values[Psiqp_idx] = (-Psiqp + (Xqp - Xl) * Iq + Edp) / Tqopp - dPsiqp;

    // D-axis transient EMF
    double TempQ = (Xqp - Xqpp) / ((Xqp - Xl) * (Xqp - Xl)) * (-Psiqp + (Xqp - Xl) * Iq + Edp);
    values[Edp_idx] = (-Edp + (Xq - Xqp) * (Iq - TempQ)  + (Xq - Xl)/(Xd - Xl)*Psiqpp*Sat(Psidpp,Psiqpp)/Psipp) / Tqop - dEdp; 

  }
  
  return true;
}

/**
 * Set Jacobian values
 * @param values a 2-d array of Jacobian block for the bus
 */
bool GenrouGen::setJacobian(gridpack::ComplexType **values)
{
  // The voltage variables are first two variable in the bus array
  int VD_idx = 0; /* Row/col number for bus voltage VD variable */
  int VQ_idx = 1; /* Row/col number for bus voltage VQ variable */
  // The network current balance equations are the first two in the set
  // of equations at the bus. Note that DSim orders the reactive current (IQ)
  // first and then the real current (ID). Hence, IGQ is ordered first and then IGD
  int IGQ_idx = 0; /* Row/col location for IGQ equations */
  int IGD_idx = 1; /* Row/col location for IGD equations */

  // Offsets for variables in the bus variables array
  int delta_idx = offsetb;
  int dw_idx    = offsetb+1;
  int Eqp_idx   = offsetb+2;
  int Psidp_idx = offsetb+3;
  int Psiqp_idx = offsetb+4;
  int Edp_idx   = offsetb+5;

  // Generator equations
  double theta = delta - PI/2.0;
  double Vdterm = VD * cos(theta) + VQ * sin(theta);
  double Vqterm = VD *-sin(theta) + VQ * cos(theta);
  
  double dVdterm_dVD = cos(theta), dVdterm_dVQ = sin(theta);
  double dVqterm_dVD = -sin(theta), dVqterm_dVQ = cos(theta);
  
  double dVdterm_ddelta = -VD*sin(theta) + VQ*cos(theta);
  double dVqterm_ddelta = -VD*cos(theta) - VQ*sin(theta);
  
  double tempd1,tempd2,tempq1,tempq2;
  tempd1 = (Xdpp - Xl)/(Xdp - Xl);
  tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
  tempq1 = (Xdpp - Xl)/(Xqp - Xl);
  tempq2 = (Xqp - Xdpp)/(Xqp - Xl);
  
  double Psidpp =  tempd1*Eqp + tempd2*Psidp;
  double Psiqpp = -tempq1*Edp - tempq2*Psiqp;
  double Psipp  = sqrt(Psidpp*Psidpp + Psiqpp*Psiqpp);
  
  double dPsidpp_dEqp   =  tempd1;
  double dPsidpp_dPsidp =  tempd2;
  double dPsiqpp_dEdp   = -tempq1;
  double dPsiqpp_dPsiqp = -tempq2;
  
  double dPsipp_dPsidpp = Psidpp/Psipp;
  double dPsipp_dPsiqpp = Psiqpp/Psipp;
  
  double dPsipp_dEqp   = dPsipp_dPsidpp*dPsidpp_dEqp;
  double dPsipp_dPsidp = dPsipp_dPsidpp*dPsidpp_dPsidp;
  double dPsipp_dEdp   = dPsipp_dPsiqpp*dPsiqpp_dEdp;
  double dPsipp_dPsiqp = dPsipp_dPsiqpp*dPsiqpp_dPsiqp;
  
  double Vd = -Psiqpp*(1 + dw);
  double Vq =  Psidpp*(1 + dw);
  
  double dVd_dPsiqpp = -(1 + dw);
  double dVd_ddw     = -Psiqpp;
  double dVq_dPsidpp =  (1 + dw);
  double dVq_ddw     =  Psidpp;
  
  double dVd_dEdp    = dVd_dPsiqpp*dPsiqpp_dEdp;
  double dVd_dPsiqp  = dVd_dPsiqpp*dPsiqpp_dPsiqp;
  double dVq_dEqp    = dVq_dPsidpp*dPsidpp_dEqp;
  double dVq_dPsidp  = dVq_dPsidpp*dPsidpp_dPsidp;
  
  // dq Axis currents
  Id = (Vd-Vdterm)*G - (Vq-Vqterm)*B;
  Iq = (Vd-Vdterm)*B + (Vq-Vqterm)*G;
  
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
  
  double dId_dEdp = dId_dVd*dVd_dEdp;
  double dIq_dEdp = dIq_dVd*dVd_dEdp;
  
  double dId_dPsiqp = dId_dVd*dVd_dPsiqp;
  double dIq_dPsiqp = dIq_dVd*dVd_dPsiqp;
  
  double dId_dVD = dId_dVdterm*dVdterm_dVD + dId_dVqterm*dVqterm_dVD;
  double dId_dVQ = dId_dVdterm*dVdterm_dVQ + dId_dVqterm*dVqterm_dVQ;
  
  double dIq_dVD = dIq_dVdterm*dVdterm_dVD + dIq_dVqterm*dVqterm_dVD;
  double dIq_dVQ = dIq_dVdterm*dVdterm_dVQ + dIq_dVqterm*dVqterm_dVQ;
  
  // Electrical torque
  double Telec = Psidpp*Iq - Psiqpp*Id;
  
  double dTelec_ddelta = Psidpp*dIq_ddelta - Psiqpp*dId_ddelta;
  double dTelec_ddw     = Psidpp*dIq_ddw    - Psiqpp*dId_ddw;
  double dTelec_dEqp    = Psidpp*dIq_dEqp   + dPsidpp_dEqp*Iq   - Psiqpp*dId_dEqp;
  double dTelec_dPsidp  = Psidpp*dIq_dPsidp + dPsidpp_dPsidp*Iq - Psiqpp*dId_dPsidp;
  double dTelec_dEdp    = Psidpp*dIq_dEdp   - (Psiqpp*dId_dEdp + dPsiqpp_dEdp*Id);
  double dTelec_dPsiqp  = Psidpp*dIq_dPsiqp - (Psiqpp*dId_dPsiqp + dPsiqpp_dPsiqp*Id);
  
  double dTelec_dVD = Psidpp*dIq_dVD - Psiqpp*dId_dVD;
  double dTelec_dVQ = Psidpp*dIq_dVQ - Psiqpp*dId_dVQ;
  
  // Field current
  double LadIfd = getFieldCurrent();
  
  if(p_mode == FAULT_EVAL) {
    // Generator variables held constant
    // dF_dX
    // Set diagonal values to 1.0
    values[delta_idx][delta_idx] = 1.0;
    values[dw_idx][dw_idx] = 1.0;
    values[Eqp_idx][Eqp_idx] = 1.0;
    values[Psidp_idx][Psidp_idx] = 1.0;
    values[Psiqp_idx][Psiqp_idx] = 1.0;
    values[Edp_idx][Edp_idx] = 1.0;

    // dG_dV
  } else {
    // Partial derivatives of ddelta_dt equation
    values[delta_idx][delta_idx] = -shift;
    values[dw_idx][delta_idx]    = OMEGA_S;
    
    // Partial derivatives of dw_dt equation
    double const1 = 1/(2*H);

    values[delta_idx][dw_idx] = const1*-dTelec_ddelta;
    values[dw_idx][dw_idx] = const1*(-(Pmech-D*dw)/((1+dw)*(1+dw)) - D/(1+dw) - dTelec_ddw) - shift;
    values[Eqp_idx][dw_idx] = const1*-dTelec_dEqp;
    values[Psidp_idx][dw_idx] = const1*-dTelec_dPsidp;
    values[Psiqp_idx][dw_idx] = const1*-dTelec_dPsiqp;
    values[Edp_idx][dw_idx]   = const1*-dTelec_dEdp;

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

    // Partial derivatives of dEqp_dt equation
    // Note: No saturation considered yet
    double const2 =  (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl));
    double dLadIfd_ddelta =  (Xd - Xdp)*(dId_ddelta + const2*(-(Xdp - Xl)*dId_ddelta));
    double dLadIfd_ddw    =  (Xd - Xdp)*(dId_ddw + const2*(-(Xdp - Xl)*dId_ddw));
    double dLadIfd_dEqp   = 1.0 + (Xd - Xdp)*(dId_dEqp + const2*(1.0 - (Xdp - Xl)*dId_dEqp));
    double dLadIfd_dPsidp =  (Xd - Xdp)*(dId_dPsidp + const2*(-1.0 -(Xdp - Xl)*dId_dPsidp));
    double dLadIfd_dPsiqp =  (Xd - Xdp)*(dId_dPsiqp + const2*(-(Xdp - Xl)*dId_dPsiqp));
    double dLadIfd_dEdp   =  (Xd - Xdp)*(dId_dEdp + const2*(-(Xdp - Xl)*dId_dEdp));
    double dLadIfd_dVD   =  (Xd - Xdp)*(dId_dVD + const2*(-(Xdp - Xl)*dId_dVD));
    double dLadIfd_dVQ   =  (Xd - Xdp)*(dId_dVQ + const2*(-(Xdp - Xl)*dId_dVQ));

    values[delta_idx][Eqp_idx] = -dLadIfd_ddelta/Tdop;
    values[dw_idx][Eqp_idx]    = -dLadIfd_ddw/Tdop;
    values[Eqp_idx][Eqp_idx]   = -dLadIfd_dEqp/Tdop - shift;
    values[Psidp_idx][Eqp_idx] = -dLadIfd_dPsidp/Tdop;
    values[Psiqp_idx][Eqp_idx] = -dLadIfd_dPsiqp/Tdop;
    values[Edp_idx][Eqp_idx]   = -dLadIfd_dEdp/Tdop;

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

    // Partial derivatives of dPsidp_dt equation
    values[delta_idx][Psidp_idx] = -(Xdp - Xl)*dId_ddelta/Tdopp;
    values[dw_idx][Psidp_idx]    = -(Xdp - Xl)*dId_ddw/Tdopp;
    values[Eqp_idx][Psidp_idx]   = (-(Xdp - Xl)*dId_dEqp + 1.0)/Tdopp;
    values[Psidp_idx][Psidp_idx] = (-1.0 -(Xdp - Xl)*dId_dPsidp)/Tdopp - shift;
    values[Psiqp_idx][Psidp_idx] = -(Xdp - Xl)*dId_dPsiqp/Tdopp;
    values[Edp_idx][Psidp_idx]   = -(Xdp - Xl)*dId_dEdp/Tdopp;

    values[VD_idx][Psidp_idx]    = -(Xdp - Xl)*dId_dVD/Tdopp;
    values[VQ_idx][Psidp_idx]    = -(Xdp - Xl)*dId_dVQ/Tdopp;

    // Partial derivatives of dPsiqp_dt equation
    values[delta_idx][Psiqp_idx] = (Xqp - Xl)*dIq_ddelta/Tqopp;
    values[dw_idx][Psiqp_idx]    = (Xqp - Xl)*dIq_ddw/Tqopp;
    values[Eqp_idx][Psiqp_idx]   = (Xqp - Xl)*dIq_dEqp/Tqopp;
    values[Psidp_idx][Psiqp_idx] = (Xqp - Xl)*dIq_dPsidp/Tqopp;
    values[Psiqp_idx][Psiqp_idx] = (-1.0 + (Xqp - Xl)*dIq_dPsiqp)/Tqopp - shift;
    values[Edp_idx][Psiqp_idx]   = ((Xqp - Xl)*dIq_dEdp + 1.0)/Tqopp;

    values[VD_idx][Psiqp_idx]    = (Xqp - Xl)*dIq_dVD/Tqopp;
    values[VQ_idx][Psiqp_idx]    = (Xqp - Xl)*dIq_dVQ/Tqopp;

    // Partial derivatives of dEdp_dt equation
    double const3 = (Xqp - Xqpp) / ((Xqp - Xl) * (Xqp - Xl));
    values[delta_idx][Edp_idx] = (Xq - Xqp)*(dIq_ddelta - const3*(Xqp - Xl)*dIq_ddelta)/Tqop;
    values[dw_idx][Edp_idx]    = (Xq - Xqp)*(dIq_ddw - const3*(Xqp - Xl)*dIq_ddw)/Tqop;
    values[Eqp_idx][Edp_idx]   = (Xq - Xqp)*(dIq_dEqp - const3*(Xqp - Xl)*dIq_dEqp)/Tqop;
    values[Psidp_idx][Edp_idx] = (Xq - Xqp)*(dIq_dPsidp - const3*(Xqp - Xl)*dIq_dPsidp)/Tqop;
    values[Psiqp_idx][Edp_idx] = (Xq - Xqp)*(dIq_dPsiqp - const3*(-1.0 + (Xqp - Xl)*dIq_dPsiqp))/Tqop;
    values[Edp_idx][Edp_idx]   = (-1.0 + (Xq - Xqp)*(dIq_dEdp - const3*((Xqp - Xl)*dIq_dEdp + 1.0)))/Tqop - shift;

    values[VD_idx][Edp_idx]    = (Xq - Xqp)*(dIq_dVD - const3*(Xqp - Xl)*dIq_dVD)/Tqop;
    values[VQ_idx][Edp_idx]    = (Xq - Xqp)*(dIq_dVQ - const3*(Xqp - Xl)*dIq_dVQ)/Tqop;
  }
   
  // Partial derivatives of generator currents IGQ and IGD w.r.t x and V
  values[delta_idx][IGQ_idx]  =  dId_ddelta*sin(theta) + Id*cos(theta)
    + dIq_ddelta*cos(theta) + Iq*-sin(theta); 
  values[dw_idx][IGQ_idx]     =  dId_ddw*sin(theta) + dIq_ddw*cos(theta);
  values[Eqp_idx][IGQ_idx]    =  dId_dEqp*sin(theta) + dIq_dEqp*cos(theta);
  values[Psidp_idx][IGQ_idx]  =  dId_dPsidp*sin(theta) + dIq_dPsidp*cos(theta);
  values[Psiqp_idx][IGQ_idx]  =  dId_dPsiqp*sin(theta) + dIq_dPsiqp*cos(theta);
  values[Edp_idx][IGQ_idx]    =  dId_dEdp*sin(theta) + dIq_dEdp*cos(theta);

  // Note the values are added, since different components (bus, load, etc.) add contributions to these locations
  values[VD_idx][IGQ_idx]    +=  dId_dVD*sin(theta) + dIq_dVD*cos(theta);
  values[VQ_idx][IGQ_idx]    +=  dId_dVQ*sin(theta) + dIq_dVQ*cos(theta);

  values[delta_idx][IGD_idx] =  dId_ddelta*cos(theta) + Id*-sin(theta)
    - dIq_ddelta*sin(theta) - Iq*cos(theta); 
  values[dw_idx][IGD_idx]    = dId_ddw*cos(theta) - dIq_ddw*sin(theta);
  values[Eqp_idx][IGD_idx]   = dId_dEqp*cos(theta) - dIq_dEqp*sin(theta);
  values[Psidp_idx][IGD_idx] = dId_dPsidp*cos(theta) - dIq_dPsidp*sin(theta);
  values[Psiqp_idx][IGD_idx] = dId_dPsiqp*cos(theta) - dIq_dPsiqp*sin(theta);
  values[Edp_idx][IGD_idx]   = dId_dEdp*cos(theta) - dIq_dEdp*sin(theta);
  
  // Note the values are added, since different components (bus, load, etc.) add contributions to these locations
  values[VD_idx][IGD_idx]    += dId_dVD*cos(theta) - dIq_dVD*sin(theta);
  values[VQ_idx][IGD_idx]    += dId_dVQ*cos(theta) - dIq_dVQ*sin(theta);

  return true;
}

/**
 * Set Jacobian values
 * @param value_map standard map containing indices and values of matrix
 *                   elements
 */
bool GenrouGen::setJacobian(std::map<std::pair<int,int>,gridpack::ComplexType> &value_map)
{
  // The voltage variables are first two variable in the bus array
  int VD_idx = 0; /* Row/col number for bus voltage VD variable */
  int VQ_idx = 1; /* Row/col number for bus voltage VQ variable */
  // The network current balance equations are the first two in the set
  // of equations at the bus. Note that DSim orders the reactive current (IQ)
  // first and then the real current (ID). Hence, IGQ is ordered first and then IGD
  int IGQ_idx = 0; /* Row/col location for IGQ equations */
  int IGD_idx = 1; /* Row/col location for IGD equations */

  // Offsets for variables in the bus variables array
  int delta_idx = offsetb;
  int dw_idx    = offsetb+1;
  int Eqp_idx   = offsetb+2;
  int Psidp_idx = offsetb+3;
  int Psiqp_idx = offsetb+4;
  int Edp_idx   = offsetb+5;

  // Generator equations
  double theta = delta - PI/2.0;
  double Vdterm = VD * cos(theta) + VQ * sin(theta);
  double Vqterm = VD *-sin(theta) + VQ * cos(theta);

  double dVdterm_dVD = cos(theta), dVdterm_dVQ = sin(theta);
  double dVqterm_dVD = -sin(theta), dVqterm_dVQ = cos(theta);

  double dVdterm_ddelta = -VD*sin(theta) + VQ*cos(theta);
  double dVqterm_ddelta = -VD*cos(theta) - VQ*sin(theta);

  double tempd1,tempd2,tempq1,tempq2;
  tempd1 = (Xdpp - Xl)/(Xdp - Xl);
  tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
  tempq1 = (Xdpp - Xl)/(Xqp - Xl);
  tempq2 = (Xqp - Xdpp)/(Xqp - Xl);

  double Psidpp =  tempd1*Eqp + tempd2*Psidp;
  double Psiqpp = -tempq1*Edp - tempq2*Psiqp;
  double Psipp  = sqrt(Psidpp*Psidpp + Psiqpp*Psiqpp);

  double dPsidpp_dEqp   =  tempd1;
  double dPsidpp_dPsidp =  tempd2;
  double dPsiqpp_dEdp   = -tempq1;
  double dPsiqpp_dPsiqp = -tempq2;

  double dPsipp_dPsidpp = Psidpp/Psipp;
  double dPsipp_dPsiqpp = Psiqpp/Psipp;

  double dPsipp_dEqp   = dPsipp_dPsidpp*dPsidpp_dEqp;
  double dPsipp_dPsidp = dPsipp_dPsidpp*dPsidpp_dPsidp;
  double dPsipp_dEdp   = dPsipp_dPsiqpp*dPsiqpp_dEdp;
  double dPsipp_dPsiqp = dPsipp_dPsiqpp*dPsiqpp_dPsiqp;

  double Vd = -Psiqpp*(1 + dw);
  double Vq =  Psidpp*(1 + dw);

  double dVd_dPsiqpp = -(1 + dw);
  double dVd_ddw     = -Psiqpp;
  double dVq_dPsidpp =  (1 + dw);
  double dVq_ddw     =  Psidpp;

  double dVd_dEdp    = dVd_dPsiqpp*dPsiqpp_dEdp;
  double dVd_dPsiqp  = dVd_dPsiqpp*dPsiqpp_dPsiqp;
  double dVq_dEqp    = dVq_dPsidpp*dPsidpp_dEqp;
  double dVq_dPsidp  = dVq_dPsidpp*dPsidpp_dPsidp;

  // dq Axis currents
  Id = (Vd-Vdterm)*G - (Vq-Vqterm)*B;
  Iq = (Vd-Vdterm)*B + (Vq-Vqterm)*G;

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

  double dId_dEdp = dId_dVd*dVd_dEdp;
  double dIq_dEdp = dIq_dVd*dVd_dEdp;

  double dId_dPsiqp = dId_dVd*dVd_dPsiqp;
  double dIq_dPsiqp = dIq_dVd*dVd_dPsiqp;

  double dId_dVD = dId_dVdterm*dVdterm_dVD + dId_dVqterm*dVqterm_dVD;
  double dId_dVQ = dId_dVdterm*dVdterm_dVQ + dId_dVqterm*dVqterm_dVQ;

  double dIq_dVD = dIq_dVdterm*dVdterm_dVD + dIq_dVqterm*dVqterm_dVD;
  double dIq_dVQ = dIq_dVdterm*dVdterm_dVQ + dIq_dVqterm*dVqterm_dVQ;

  // Electrical torque
  double Telec = Psidpp*Iq - Psiqpp*Id;

  double dTelec_ddelta = Psidpp*dIq_ddelta - Psiqpp*dId_ddelta;
  double dTelec_ddw     = Psidpp*dIq_ddw    - Psiqpp*dId_ddw;
  double dTelec_dEqp    = Psidpp*dIq_dEqp   + dPsidpp_dEqp*Iq   - Psiqpp*dId_dEqp;
  double dTelec_dPsidp  = Psidpp*dIq_dPsidp + dPsidpp_dPsidp*Iq - Psiqpp*dId_dPsidp;
  double dTelec_dEdp    = Psidpp*dIq_dEdp   - (Psiqpp*dId_dEdp + dPsiqpp_dEdp*Id);
  double dTelec_dPsiqp  = Psidpp*dIq_dPsiqp - (Psiqpp*dId_dPsiqp + dPsiqpp_dPsiqp*Id);

  double dTelec_dVD = Psidpp*dIq_dVD - Psiqpp*dId_dVD;
  double dTelec_dVQ = Psidpp*dIq_dVQ - Psiqpp*dId_dVQ;

  // Field current
  double LadIfd = getFieldCurrent();

#define MAP_PAIR(_i, _j, _val) \
  std::pair<std::pair<int,int>,gridpack::ComplexType>( \
      std::pair<int,int>(_j,_i),_val)
  if(p_mode == FAULT_EVAL) {
    // Generator variables held constant
    // dF_dX
    // Set diagonal values to 1.0
#if 0
    value_map.insert(MAP_PAIR(delta_idx,delta_idx,1.0));
    value_map.insert(MAP_PAIR(dw_idx,dw_idx,1.0));
    value_map.insert(MAP_PAIR(Eqp_idx,Eqp_idx,1.0));
    value_map.insert(MAP_PAIR(Psidp_idx,Psidp_idx,1.0));
    value_map.insert(MAP_PAIR(Psiqp_idx,Psiqp_idx,1.0));
    value_map.insert(MAP_PAIR(Edp_idx,Edp_idx,1.0));
#endif
    // Partial derivatives of ddelta_dt equation
    value_map.insert(MAP_PAIR(delta_idx,delta_idx,1.0));
    value_map.insert(MAP_PAIR(dw_idx,delta_idx,0.0));

    value_map.insert(MAP_PAIR(delta_idx,dw_idx,0.0));
    value_map.insert(MAP_PAIR(dw_idx,dw_idx,1.0));
    value_map.insert(MAP_PAIR(Eqp_idx,dw_idx,0.0));
    value_map.insert(MAP_PAIR(Psidp_idx,dw_idx,0.0));
    value_map.insert(MAP_PAIR(Psiqp_idx,dw_idx,0.0));
    value_map.insert(MAP_PAIR(Edp_idx,dw_idx,0.0));

    value_map.insert(MAP_PAIR(VD_idx,dw_idx,0.0));
    value_map.insert(MAP_PAIR(VQ_idx,dw_idx,0.0));

    // Add Pmech contributions
    if(p_hasGovernor) {
      int nxgov,i;
      p_governor->vectorSize(&nxgov);
      p_governor->getMechanicalPowerPartialDerivatives(xgov_loc,dPmech_dxgov);

      /* Partials w.r.t. governor mechanical power Pmech */
      for(i=0; i < nxgov; i++) {
        value_map.insert(MAP_PAIR(xgov_loc[i],dw_idx,0.0));
      }
    }

    // Partial derivatives of dEqp_dt equation
    // Note: No saturation considered yet
    value_map.insert(MAP_PAIR(delta_idx,Eqp_idx,0.0));
    value_map.insert(MAP_PAIR(dw_idx,Eqp_idx,0.0));
    value_map.insert(MAP_PAIR(Eqp_idx,Eqp_idx,1.0));
    value_map.insert(MAP_PAIR(Psidp_idx,Eqp_idx,0.0));
    value_map.insert(MAP_PAIR(Psiqp_idx,Eqp_idx,0.0));
    value_map.insert(MAP_PAIR(Edp_idx,Eqp_idx,0.0));

    value_map.insert(MAP_PAIR(VD_idx,Eqp_idx,0.0));
    value_map.insert(MAP_PAIR(VQ_idx,Eqp_idx,0.0));

    // Add Efd contributions
    if(p_hasExciter) {
      int nexc,i;
      p_exciter->vectorSize(&nexc);
      p_exciter->getFieldVoltagePartialDerivatives(xexc_loc,dEfd_dxexc,
          dEfd_dxgen);

      /* Partials w.r.t. exciter variables */
      for(i=0; i < nexc; i++) {
        value_map.insert(MAP_PAIR(xexc_loc[i],Eqp_idx,0.0));
      }

      /* Partials contributions for Efd w.r.t. generator variables. Note
       * that this may be because the exciter Efd calculation uses the
       * field current LadIfd (which is a function of generator variables)
       */
      for(i=0; i < nxgen; i++) {
        std::map<std::pair<int,int>,gridpack::ComplexType>::iterator it;
        it = value_map.find(std::pair<int,int>(offsetb+i,Eqp_idx));
        if (it == value_map.end()) {
        } else {
          value_map.insert(MAP_PAIR(VQ_idx,Eqp_idx,0.0));
        }
      }
    }

    // Partial derivatives of dPsidp_dt equation
    value_map.insert(MAP_PAIR(delta_idx,Psidp_idx,0.0));
    value_map.insert(MAP_PAIR(dw_idx,Psidp_idx,0.0));
    value_map.insert(MAP_PAIR(Eqp_idx,Psidp_idx,0.0));
    value_map.insert(MAP_PAIR(Psidp_idx,Psidp_idx,1.0));
    value_map.insert(MAP_PAIR(Psiqp_idx,Psidp_idx,0.0));
    value_map.insert(MAP_PAIR(Edp_idx,Psidp_idx,0.0));

    value_map.insert(MAP_PAIR(VD_idx,Psidp_idx,0.0));
    value_map.insert(MAP_PAIR(VQ_idx,Psidp_idx,0.0));

    // Partial derivatives of dPsiqp_dt equation
    value_map.insert(MAP_PAIR(delta_idx,Psiqp_idx,0.0));
    value_map.insert(MAP_PAIR(dw_idx,Psiqp_idx,0.0));
    value_map.insert(MAP_PAIR(Eqp_idx,Psiqp_idx,0.0));
    value_map.insert(MAP_PAIR(Psidp_idx,Psiqp_idx,0.0));
    value_map.insert(MAP_PAIR(Psiqp_idx,Psiqp_idx,1.0));
    value_map.insert(MAP_PAIR(Edp_idx,Psiqp_idx,0.0));

    value_map.insert(MAP_PAIR(VD_idx,Psiqp_idx,0.0));
    value_map.insert(MAP_PAIR(VQ_idx,Psiqp_idx,0.0));

    // Partial derivatives of dEdp_dt equation
    value_map.insert(MAP_PAIR(delta_idx,Edp_idx,0.0));
    value_map.insert(MAP_PAIR(dw_idx,Edp_idx,0.0));
    value_map.insert(MAP_PAIR(Eqp_idx,Edp_idx,0.0));
    value_map.insert(MAP_PAIR(Psidp_idx,Edp_idx,0.0));
    value_map.insert(MAP_PAIR(Psiqp_idx,Edp_idx,0.0));
    value_map.insert(MAP_PAIR(Edp_idx,Edp_idx,1.0));

    value_map.insert(MAP_PAIR(VD_idx,Edp_idx,0.0));
    value_map.insert(MAP_PAIR(VQ_idx,Edp_idx,0.0));

    // Partial derivatives of generator currents IGQ and IGD w.r.t x and V
    value_map.insert(MAP_PAIR(delta_idx,IGQ_idx,0.0));
    value_map.insert(MAP_PAIR(dw_idx,IGQ_idx,0.0));
    value_map.insert(MAP_PAIR(Eqp_idx,IGQ_idx,0.0));
    value_map.insert(MAP_PAIR(Psidp_idx,IGQ_idx,0.0));
    value_map.insert(MAP_PAIR(Psiqp_idx,IGQ_idx,0.0));
    value_map.insert(MAP_PAIR(Edp_idx,IGQ_idx,0.0));

    value_map.insert(MAP_PAIR(delta_idx,IGD_idx,0.0));
    value_map.insert(MAP_PAIR(dw_idx,IGD_idx,0.0));
    value_map.insert(MAP_PAIR(Eqp_idx,IGD_idx,0.0));
    value_map.insert(MAP_PAIR(Psidp_idx,IGD_idx,0.0));
    value_map.insert(MAP_PAIR(Psiqp_idx,IGD_idx,0.0));
    value_map.insert(MAP_PAIR(Edp_idx,IGD_idx,0.0));
    // dG_dV
  } else {
    gridpack::ComplexType zval;
    // Partial derivatives of ddelta_dt equation
    value_map.insert(MAP_PAIR(delta_idx,delta_idx,-shift));
    value_map.insert(MAP_PAIR(dw_idx,delta_idx,OMEGA_S));

    // Partial derivatives of dw_dt equation
    double const1 = 1/(2*H);

    zval = const1*-dTelec_ddelta;
    value_map.insert(MAP_PAIR(delta_idx,dw_idx,zval));
    zval = const1*(-(Pmech-D*dw)/((1+dw)*(1+dw)) - D/(1+dw) - dTelec_ddw)
      - shift;
    value_map.insert(MAP_PAIR(dw_idx,dw_idx,zval));
    zval = const1*-dTelec_dEqp;
    value_map.insert(MAP_PAIR(Eqp_idx,dw_idx,zval));
    zval = const1*-dTelec_dPsidp;
    value_map.insert(MAP_PAIR(Psidp_idx,dw_idx,zval));
    zval = const1*-dTelec_dPsiqp;
    value_map.insert(MAP_PAIR(Psiqp_idx,dw_idx,zval));
    zval   = const1*-dTelec_dEdp;
    value_map.insert(MAP_PAIR(Edp_idx,dw_idx,zval));

    zval = const1*-dTelec_dVD;
    value_map.insert(MAP_PAIR(VD_idx,dw_idx,zval));
    zval = const1*-dTelec_dVQ;
    value_map.insert(MAP_PAIR(VQ_idx,dw_idx,zval));

    // Add Pmech contributions
    if(p_hasGovernor) {
      int nxgov,i;
      p_governor->vectorSize(&nxgov);
      p_governor->getMechanicalPowerPartialDerivatives(xgov_loc,dPmech_dxgov);

      /* Partials w.r.t. governor mechanical power Pmech */
      for(i=0; i < nxgov; i++) {
        zval = const1*dPmech_dxgov[i]/(1+dw);
        value_map.insert(MAP_PAIR(xgov_loc[i],dw_idx,zval));
      }
    }

    // Partial derivatives of dEqp_dt equation
    // Note: No saturation considered yet
    double const2 =  (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl));
    double dLadIfd_ddelta =  (Xd - Xdp)*(dId_ddelta + const2*(-(Xdp - Xl)*dId_ddelta));
    double dLadIfd_ddw    =  (Xd - Xdp)*(dId_ddw + const2*(-(Xdp - Xl)*dId_ddw));
    double dLadIfd_dEqp   = 1.0 + (Xd - Xdp)*(dId_dEqp + const2*(1.0 - (Xdp - Xl)*dId_dEqp));
    double dLadIfd_dPsidp =  (Xd - Xdp)*(dId_dPsidp + const2*(-1.0 -(Xdp - Xl)*dId_dPsidp));
    double dLadIfd_dPsiqp =  (Xd - Xdp)*(dId_dPsiqp + const2*(-(Xdp - Xl)*dId_dPsiqp));
    double dLadIfd_dEdp   =  (Xd - Xdp)*(dId_dEdp + const2*(-(Xdp - Xl)*dId_dEdp));
    double dLadIfd_dVD   =  (Xd - Xdp)*(dId_dVD + const2*(-(Xdp - Xl)*dId_dVD));
    double dLadIfd_dVQ   =  (Xd - Xdp)*(dId_dVQ + const2*(-(Xdp - Xl)*dId_dVQ));

    zval = -dLadIfd_ddelta/Tdop;
    value_map.insert(MAP_PAIR(delta_idx,Eqp_idx,zval));
    zval    = -dLadIfd_ddw/Tdop;
    value_map.insert(MAP_PAIR(dw_idx,Eqp_idx,zval));
    zval   = -dLadIfd_dEqp/Tdop - shift;
    value_map.insert(MAP_PAIR(Eqp_idx,Eqp_idx,zval));
    zval = -dLadIfd_dPsidp/Tdop;
    value_map.insert(MAP_PAIR(Psidp_idx,Eqp_idx,zval));
    zval = -dLadIfd_dPsiqp/Tdop;
    value_map.insert(MAP_PAIR(Psiqp_idx,Eqp_idx,zval));
    zval   = -dLadIfd_dEdp/Tdop;
    value_map.insert(MAP_PAIR(Edp_idx,Eqp_idx,zval));

    zval    = -dLadIfd_dVD/Tdop;
    value_map.insert(MAP_PAIR(VD_idx,Eqp_idx,zval));
    zval    = -dLadIfd_dVQ/Tdop;
    value_map.insert(MAP_PAIR(VQ_idx,Eqp_idx,zval));

    // Add Efd contributions
    if(p_hasExciter) {
      int nexc,i;
      p_exciter->vectorSize(&nexc);
      p_exciter->getFieldVoltagePartialDerivatives(xexc_loc,dEfd_dxexc,
          dEfd_dxgen);

      /* Partials w.r.t. exciter variables */
      for(i=0; i < nexc; i++) {
        zval = dEfd_dxexc[i]/Tdop;
        value_map.insert(MAP_PAIR(xexc_loc[i],Eqp_idx,zval));
      }

      /* Partials contributions for Efd w.r.t. generator variables. Note
       * that this may be because the exciter Efd calculation uses the
       * field current LadIfd (which is a function of generator variables)
       */
      for(i=0; i < nxgen; i++) {
        zval = dEfd_dxgen[i]/Tdop;
        std::map<std::pair<int,int>,gridpack::ComplexType>::iterator it;
        it = value_map.find(std::pair<int,int>(offsetb+i,Eqp_idx));
        if (it != value_map.end()) {
          it->second += zval;
        } else {
          value_map.insert(MAP_PAIR(offsetb+i,Eqp_idx,zval));
        }
      }
    }

    // Partial derivatives of dPsidp_dt equation
    zval = -(Xdp - Xl)*dId_ddelta/Tdopp;
    value_map.insert(MAP_PAIR(delta_idx,Psidp_idx,zval));
    zval    = -(Xdp - Xl)*dId_ddw/Tdopp;
    value_map.insert(MAP_PAIR(dw_idx,Psidp_idx,zval));
    zval   = (-(Xdp - Xl)*dId_dEqp + 1.0)/Tdopp;
    value_map.insert(MAP_PAIR(Eqp_idx,Psidp_idx,zval));
    zval = (-1.0 -(Xdp - Xl)*dId_dPsidp)/Tdopp - shift;
    value_map.insert(MAP_PAIR(Psidp_idx,Psidp_idx,zval));
    zval = -(Xdp - Xl)*dId_dPsiqp/Tdopp;
    value_map.insert(MAP_PAIR(Psiqp_idx,Psidp_idx,zval));
    zval   = -(Xdp - Xl)*dId_dEdp/Tdopp;
    value_map.insert(MAP_PAIR(Edp_idx,Psidp_idx,zval));

    zval    = -(Xdp - Xl)*dId_dVD/Tdopp;
    value_map.insert(MAP_PAIR(VD_idx,Psidp_idx,zval));
    zval    = -(Xdp - Xl)*dId_dVQ/Tdopp;
    value_map.insert(MAP_PAIR(VQ_idx,Psidp_idx,zval));

    // Partial derivatives of dPsiqp_dt equation
    zval = (Xqp - Xl)*dIq_ddelta/Tqopp;
    value_map.insert(MAP_PAIR(delta_idx,Psiqp_idx,zval));
    zval = (Xqp - Xl)*dIq_ddw/Tqopp;
    value_map.insert(MAP_PAIR(dw_idx,Psiqp_idx,zval));
    zval   = (Xqp - Xl)*dIq_dEqp/Tqopp;
    value_map.insert(MAP_PAIR(Eqp_idx,Psiqp_idx,zval));
    zval = (Xqp - Xl)*dIq_dPsidp/Tqopp;
    value_map.insert(MAP_PAIR(Psidp_idx,Psiqp_idx,zval));
    zval = (-1.0 + (Xqp - Xl)*dIq_dPsiqp)/Tqopp - shift;
    value_map.insert(MAP_PAIR(Psiqp_idx,Psiqp_idx,zval));
    zval   = ((Xqp - Xl)*dIq_dEdp + 1.0)/Tqopp;
    value_map.insert(MAP_PAIR(Edp_idx,Psiqp_idx,zval));

    zval = (Xqp - Xl)*dIq_dVD/Tqopp;
    value_map.insert(MAP_PAIR(VD_idx,Psiqp_idx,zval));
    zval = (Xqp - Xl)*dIq_dVQ/Tqopp;
    value_map.insert(MAP_PAIR(VQ_idx,Psiqp_idx,zval));

    // Partial derivatives of dEdp_dt equation
    double const3 = (Xqp - Xqpp) / ((Xqp - Xl) * (Xqp - Xl));
    zval = (Xq - Xqp)*(dIq_ddelta - const3*(Xqp - Xl)*dIq_ddelta)/Tqop;
    value_map.insert(MAP_PAIR(delta_idx,Edp_idx,zval));
    zval = (Xq - Xqp)*(dIq_ddw - const3*(Xqp - Xl)*dIq_ddw)/Tqop;
    value_map.insert(MAP_PAIR(dw_idx,Edp_idx,zval));
    zval = (Xq - Xqp)*(dIq_dEqp - const3*(Xqp - Xl)*dIq_dEqp)/Tqop;
    value_map.insert(MAP_PAIR(Eqp_idx,Edp_idx,zval));
    zval = (Xq - Xqp)*(dIq_dPsidp - const3*(Xqp - Xl)*dIq_dPsidp)/Tqop;
    value_map.insert(MAP_PAIR(Psidp_idx,Edp_idx,zval));
    zval = (Xq - Xqp)*(dIq_dPsiqp - const3*(-1.0 + (Xqp - Xl)*dIq_dPsiqp))/Tqop;
    value_map.insert(MAP_PAIR(Psiqp_idx,Edp_idx,zval));
    zval = (-1.0 + (Xq - Xqp)*(dIq_dEdp - const3*((Xqp - Xl)*dIq_dEdp
            + 1.0)))/Tqop - shift;
    value_map.insert(MAP_PAIR(Edp_idx,Edp_idx,zval));

    zval = (Xq - Xqp)*(dIq_dVD - const3*(Xqp - Xl)*dIq_dVD)/Tqop;
    value_map.insert(MAP_PAIR(VD_idx,Edp_idx,zval));
    zval = (Xq - Xqp)*(dIq_dVQ - const3*(Xqp - Xl)*dIq_dVQ)/Tqop;
    value_map.insert(MAP_PAIR(VQ_idx,Edp_idx,zval));

    // Partial derivatives of generator currents IGQ and IGD w.r.t x and V
    zval = dId_ddelta*sin(theta) + Id*cos(theta)
      + dIq_ddelta*cos(theta) + Iq*-sin(theta); 
    value_map.insert(MAP_PAIR(delta_idx,IGQ_idx,zval));
    zval = dId_ddw*sin(theta) + dIq_ddw*cos(theta);
    value_map.insert(MAP_PAIR(dw_idx,IGQ_idx,zval));
    zval = dId_dEqp*sin(theta) + dIq_dEqp*cos(theta);
    value_map.insert(MAP_PAIR(Eqp_idx,IGQ_idx,zval));
    zval = dId_dPsidp*sin(theta) + dIq_dPsidp*cos(theta);
    value_map.insert(MAP_PAIR(Psidp_idx,IGQ_idx,zval));
    zval = dId_dPsiqp*sin(theta) + dIq_dPsiqp*cos(theta);
    value_map.insert(MAP_PAIR(Psiqp_idx,IGQ_idx,zval));
    zval = dId_dEdp*sin(theta) + dIq_dEdp*cos(theta);
    value_map.insert(MAP_PAIR(Edp_idx,IGQ_idx,zval));

    // Note the values are added, since different components (bus, load, etc.) add contributions to these locations
    std::map<std::pair<int,int>,gridpack::ComplexType>::iterator it;
    it = value_map.find(std::pair<int,int>(VD_idx,IGQ_idx));
    if (it != value_map.end())  {
      zval = dId_dVD*sin(theta) + dIq_dVD*cos(theta);
      it->second += zval;
    }
    it = value_map.find(std::pair<int,int>(VQ_idx,IGQ_idx));
    if (it != value_map.end())  {
      zval = dId_dVQ*sin(theta) + dIq_dVQ*cos(theta);
      it->second += zval;
    }

    zval = dId_ddelta*cos(theta) + Id*-sin(theta)
      - dIq_ddelta*sin(theta) - Iq*cos(theta); 
    value_map.insert(MAP_PAIR(delta_idx,IGD_idx,zval));
    zval = dId_ddw*cos(theta) - dIq_ddw*sin(theta);
    value_map.insert(MAP_PAIR(dw_idx,IGD_idx,zval));
    zval = dId_dEqp*cos(theta) - dIq_dEqp*sin(theta);
    value_map.insert(MAP_PAIR(Eqp_idx,IGD_idx,zval));
    zval = dId_dPsidp*cos(theta) - dIq_dPsidp*sin(theta);
    value_map.insert(MAP_PAIR(Psidp_idx,IGD_idx,zval));
    zval = dId_dPsiqp*cos(theta) - dIq_dPsiqp*sin(theta);
    value_map.insert(MAP_PAIR(Psiqp_idx,IGD_idx,zval));
    zval = dId_dEdp*cos(theta) - dIq_dEdp*sin(theta);
    value_map.insert(MAP_PAIR(Edp_idx,IGD_idx,zval));

    // Note the values are added, since different components
    // (bus, load, etc.) add contributions to these locations
    it = value_map.find(std::pair<int,int>(VD_idx,IGD_idx));
    if (it != value_map.end())  {
      zval = dId_dVD*cos(theta) - dIq_dVD*sin(theta);
      it->second += zval;
    }
    it = value_map.find(std::pair<int,int>(VQ_idx,IGD_idx));
    if (it != value_map.end())  {
      zval = dId_dVQ*cos(theta) - dIq_dVQ*sin(theta);
      it->second += zval;
    }
  }
#undef MAP_PAIR
  return true;
}

/**
 * Number of rows (equations) contributed to by generator
 * @return number of rows
 */
int GenrouGen::matrixNumRows()
{
  return 6;
}

/**
 * Number of rows (equations) contributed to by generator
 * @return number of rows
 */
int GenrouGen::matrixNumCols()
{
  return 6;
}

/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values
 */
int GenrouGen::matrixNumValues()
{
  // The voltage variables are first two variable in the bus array
  int VD_idx = 0; /* Row/col number for bus voltage VD variable */
  int VQ_idx = 1; /* Row/col number for bus voltage VQ variable */
  // The network current balance equations are the first two in the set
  // of equations at the bus. Note that DSim orders the reactive current (IQ)
  // first and then the real current (ID). Hence, IGQ is ordered first and then IGD
  int IGQ_idx = 0; /* Row/col location for IGQ equations */
  int IGD_idx = 1; /* Row/col location for IGD equations */

  std::set<std::pair<int,int> > set;
  // Offsets for variables in the bus variables array
  int delta_idx = offsetb;
  int dw_idx    = offsetb+1;
  int Eqp_idx   = offsetb+2;
  int Psidp_idx = offsetb+3;
  int Psiqp_idx = offsetb+4;
  int Edp_idx   = offsetb+5;

  if(p_mode == FAULT_EVAL) {
    set.insert(std::pair<int,int>(delta_idx,delta_idx));
    set.insert(std::pair<int,int>(dw_idx,dw_idx));
    set.insert(std::pair<int,int>(Eqp_idx,Eqp_idx));
    set.insert(std::pair<int,int>(Psidp_idx,Psidp_idx));
    set.insert(std::pair<int,int>(Psiqp_idx,Psiqp_idx));
    set.insert(std::pair<int,int>(Edp_idx,Edp_idx));

    // dG_dV
  } else {
    // Partial derivatives of ddelta_dt equation
    set.insert(std::pair<int,int>(delta_idx,delta_idx));
    set.insert(std::pair<int,int>(dw_idx,delta_idx));

    // Partial derivatives of dw_dt equation
    double const1 = 1/(2*H);

    set.insert(std::pair<int,int>(delta_idx,dw_idx));
    set.insert(std::pair<int,int>(dw_idx,dw_idx));
    set.insert(std::pair<int,int>(Eqp_idx,dw_idx));
    set.insert(std::pair<int,int>(Psidp_idx,dw_idx));
    set.insert(std::pair<int,int>(Psiqp_idx,dw_idx));
    set.insert(std::pair<int,int>(Edp_idx,dw_idx));

    set.insert(std::pair<int,int>(VD_idx,dw_idx));
    set.insert(std::pair<int,int>(VQ_idx,dw_idx));

    // Add Pmech contributions
    if(p_hasGovernor) {
      int nxgov,i;
      p_governor->vectorSize(&nxgov);
      /* Partials w.r.t. governor mechanical power Pmech */
      for(i=0; i < nxgov; i++) {
        set.insert(std::pair<int,int>(xgov_loc[i],dw_idx));
      }
    }

    // Partial derivatives of dEqp_dt equation
    // Note: No saturation considered yet
    set.insert(std::pair<int,int>(delta_idx,Eqp_idx));
    set.insert(std::pair<int,int>(dw_idx,Eqp_idx));
    set.insert(std::pair<int,int>(Eqp_idx,Eqp_idx));
    set.insert(std::pair<int,int>(Psidp_idx,Eqp_idx));
    set.insert(std::pair<int,int>(Psiqp_idx,Eqp_idx));
    set.insert(std::pair<int,int>(Edp_idx,Eqp_idx));

    set.insert(std::pair<int,int>(VD_idx,Eqp_idx));
    set.insert(std::pair<int,int>(VQ_idx,Eqp_idx));

    // Add Efd contributions
    if(p_hasExciter) {
      int nexc,i;
      p_exciter->vectorSize(&nexc);

      /* Partials w.r.t. exciter variables */
      for(i=0; i < nexc; i++) {
        set.insert(std::pair<int,int>(xexc_loc[i],Eqp_idx));
      }

      /* Partials contributions for Efd w.r.t. generator variables.
       * Note that this may be because the exciter Efd calculation
       * uses the field current LadIfd (which is a function of
       * generator variables)
       */
      for(i=0; i < nxgen; i++) {
        if (set.find(std::pair<int,int>(offsetb+i,Eqp_idx)) == set.end()) {
          set.insert(std::pair<int,int>(offsetb+i,Eqp_idx));
        } 

      }
    }

    // Partial derivatives of dPsidp_dt equation
    set.insert(std::pair<int,int>(delta_idx,Psidp_idx));
    set.insert(std::pair<int,int>(dw_idx,Psidp_idx));
    set.insert(std::pair<int,int>(Eqp_idx,Psidp_idx));
    set.insert(std::pair<int,int>(Psidp_idx,Psidp_idx));
    set.insert(std::pair<int,int>(Psiqp_idx,Psidp_idx));
    set.insert(std::pair<int,int>(Edp_idx,Psidp_idx));

    set.insert(std::pair<int,int>(VD_idx,Psidp_idx));
    set.insert(std::pair<int,int>(VQ_idx,Psidp_idx));

    // Partial derivatives of dPsiqp_dt equation
    set.insert(std::pair<int,int>(delta_idx,Psiqp_idx));
    set.insert(std::pair<int,int>(dw_idx,Psiqp_idx));
    set.insert(std::pair<int,int>(Eqp_idx,Psiqp_idx));
    set.insert(std::pair<int,int>(Psidp_idx,Psiqp_idx));
    set.insert(std::pair<int,int>(Psiqp_idx,Psiqp_idx));
    set.insert(std::pair<int,int>(Edp_idx,Psiqp_idx));

    set.insert(std::pair<int,int>(VD_idx,Psiqp_idx));
    set.insert(std::pair<int,int>(VQ_idx,Psiqp_idx));

    // Partial derivatives of dEdp_dt equation
    set.insert(std::pair<int,int>(delta_idx,Edp_idx));
    set.insert(std::pair<int,int>(dw_idx,Edp_idx));
    set.insert(std::pair<int,int>(Eqp_idx,Edp_idx));
    set.insert(std::pair<int,int>(Psidp_idx,Edp_idx));
    set.insert(std::pair<int,int>(Psiqp_idx,Edp_idx));
    set.insert(std::pair<int,int>(Edp_idx,Edp_idx));

    set.insert(std::pair<int,int>(VD_idx,Edp_idx));
    set.insert(std::pair<int,int>(VQ_idx,Edp_idx));

    // Partial derivatives of generator currents IGQ and IGD w.r.t x and V
    set.insert(std::pair<int,int>(delta_idx,IGQ_idx));
    set.insert(std::pair<int,int>(dw_idx,IGQ_idx));
    set.insert(std::pair<int,int>(Eqp_idx,IGQ_idx));
    set.insert(std::pair<int,int>(Psidp_idx,IGQ_idx));
    set.insert(std::pair<int,int>(Psiqp_idx,IGQ_idx));
    set.insert(std::pair<int,int>(Edp_idx,IGQ_idx));

    // Note the matrix elements at [VD_idx,IGQ_idx],[VQ_idx,IGQ_idx] are added by
    // bus, no need to add them here. Different components (bus, load, etc.)
    // add contributions to these locations

    set.insert(std::pair<int,int>(delta_idx,IGD_idx));
    set.insert(std::pair<int,int>(dw_idx,IGD_idx));
    set.insert(std::pair<int,int>(Eqp_idx,IGD_idx));
    set.insert(std::pair<int,int>(Psidp_idx,IGD_idx));
    set.insert(std::pair<int,int>(Psiqp_idx,IGD_idx));
    set.insert(std::pair<int,int>(Edp_idx,IGD_idx));

    // Note the matrix elements at [VD_idx,IGD_idx],[VQ_idx,IGD_idx] are added by
    // bus, no need to add them here. Different components (bus, load, etc.)
    // add contributions to these locations
  }

  return set.size();
}

/**
 * Get list of matrix values contributed by generator
 * @params values list of matrix values
 * @params rows list of local row indices
 */
void GenrouGen::matrixGetValues(gridpack::ComplexType *values,
    int *rows, int *cols)
{
}

/**
 * Return the generator current injection (in rectangular form) 
 * @param [output] IGD - real part of the generator current
 * @param [output] IGQ - imaginary part of the generator current
*/
void GenrouGen::getCurrent(double *IGD, double *IGQ)
{
    double theta = delta - PI/2.0;
    double Vdterm = VD * cos(theta) + VQ * sin(theta);
    double Vqterm = VD *-sin(theta) + VQ * cos(theta);

    double tempd1,tempd2,tempq1,tempq2;
    tempd1 = (Xdpp - Xl)/(Xdp - Xl);
    tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
    tempq1 = (Xdpp - Xl)/(Xqp - Xl);
    tempq2 = (Xqp - Xdpp)/(Xqp - Xl);

    double Psidpp =  tempd1*Eqp + tempd2*Psidp;
    double Psiqpp = -tempq1*Edp - tempq2*Psiqp;

    double Vd = -Psiqpp*(1 + dw);
    double Vq =  Psidpp*(1 + dw);

    // dq Axis currents
    Id = (Vd-Vdterm)*G - (Vq-Vqterm)*B;
    Iq = (Vd-Vdterm)*B + (Vq-Vqterm)*G;

    // Generator current injections in the network
    *IGD = + Id * cos(theta) - Iq * sin(theta);
    *IGQ =   Id * sin(theta) + Iq * cos(theta);

}

double GenrouGen::getRotorSpeedDeviation()
{
  return dw;
}

int GenrouGen::getRotorSpeedDeviationLocation()
{
  return offsetb+1;
}

double GenrouGen::getFieldCurrent()
{
  double theta = delta - PI/2.0;
  double Vdterm = VD * cos(theta) + VQ * sin(theta);
  double Vqterm = VD *-sin(theta) + VQ * cos(theta);
  
  double tempd1,tempd2,tempq1,tempq2;
  tempd1 = (Xdpp - Xl)/(Xdp - Xl);
  tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
  tempq1 = (Xdpp - Xl)/(Xqp - Xl);
  tempq2 = (Xqp - Xdpp)/(Xqp - Xl);
  
  double Psidpp =  tempd1*Eqp + tempd2*Psidp;
  double Psiqpp = -tempq1*Edp - tempq2*Psiqp;
  double Psipp  = sqrt(Psidpp*Psidpp + Psiqpp*Psiqpp);
  
  double Vd = -Psiqpp*(1 + dw);
  double Vq =  Psidpp*(1 + dw);

  // dq Axis currents
  Id = (Vd-Vdterm)*G - (Vq-Vqterm)*B;
  Iq = (Vd-Vdterm)*B + (Vq-Vqterm)*G;

  // Field current
  double temp = (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl))*(Eqp - Psidp - (Xdp - Xl)*Id);
  double LadIfd = Eqp + (Xd - Xdp)*(Id + temp) +  Psidpp*Sat(Psidpp,Psiqpp)/Psipp;

  return LadIfd;
}
