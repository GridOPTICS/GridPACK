/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   genrou.cpp
 * @author Shuangshuang Jin 
 * @author Shrirang Abhyankar
 * @Last modified:   12/30/19
 *  
 * @brief  
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

}

GenrouGen::~GenrouGen(void)
{
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
  data->getValue(BUS_NUMBER, &bid);
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
  if (!data->getValue(GENERATOR_XDPP, &Xqpp, idx)) Xqpp=0.0; // Xqpp // SJin: no GENERATOR_XQPP in dictionary.hpp, read XDPP instead (Xqpp = Xdpp)
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

  P = pg / mbase;
  Q = qg / mbase;

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
  p_hasExciter = getphasExciter();
  if (p_hasExciter) {
    p_exciter = getExciter();
    p_exciter->setInitialFieldVoltage(Efd);
  }
    
  // Initialize governor mechanical power and speed deviation
  p_hasGovernor = getphasGovernor();
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
  /*if (!strcmp(signal,"standard")) {
       sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f	%12.6f	%12.6f\n",
          p_bus_id, p_ckt.c_str(), delta_1, dw_1, Eqp_1, Psidp_1, Psiqp_1, Edp_1);
  }*/
  //printf("...........bid=%d: %f\t%f\t%f\t%f\t%f\t%f\n", bid, delta, dw, Eqp, Psidp, Psiqp, Edp);
}

/**
 *  Set the number of variables for this generator model
 *  @param [output] number of variables for this model
 */
bool GenrouGen::vectorSize(int *nvar) const
{
  *nvar = 6;
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

/*void GenrouGen::setExciter(boost::shared_ptr<BaseExcModel> &exciter)
{
  p_exciter = exciter;
}

boost::shared_ptr<BaseExcModel> GenrouGen::getExciter()
{
  return p_exciter;
}*/

