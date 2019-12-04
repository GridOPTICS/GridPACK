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
  x1d = 0.0; // Rotor angle
  x2w = 0.0; // Rotor speed
  x3Eqp = 0.0; // Transient Q axis Eq 
  x4Psidp = 0.0; // Transient D axis flux
  x5Psiqp = 0.0; // Transient Q axis flux
  dx1d = 0.0;
  dx2w = 0.0;
  dx3Eqp = 0.0;
  dx4Psidp = 0.0;
  dx5Psiqp = 0.0;
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
  Tqop = 0.0;

  B = 0.0;
  G = 0.0;

  Vd = 0.0;
  Vq = 0.0;

  p_hasExciter = false;
  p_hasGovernor = false;
}

GensalGen::~GensalGen(void)
{
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
  if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.55; // S12  // SJin: no GENERATOR_XQPP in dictionary.hpp, read XDPP instead (Xqpp = Xdpp)
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
  //printf("H=%f,D=%f,Ra=%f,Xd=%f,Xq=%f,Xdp=%f,Xdpp=%f,Xl=%f,Tdop=%f,Tdopp=%f,Tqopp=%f,S10=%f,S12=%f,Xqp=%f,Xqpp=%f,Tqop=%f\n", H,D,Ra,Xd,Xq,Xdp,Xdpp,Xl,Tdop,Tdopp,Tqopp,S10,S12,Xqp,Xqpp,Tqop);


}

/**
 * Saturation function
 * @ param x
 */
double GensalGen::Sat(double x)
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
 * @param [output] values - array where initialized generator variables should be set
 */
void GensalGen::init(gridpack::ComplexType* values) 
{
  double Vterm = sqrt(VD*VD + VQ*VQ); // SJin: voltage VD and VQ come from base_gen_model.hpp
  double P, Q; // Generator real and reactive power
  double Vrterm = VD;
  double Viterm = VQ;

  // Admittance
  B = -Xdpp / (Ra * Ra + Xdpp * Xdpp);
  G = Ra / (Ra * Ra + Xdpp * Xdpp);

  P = pg / sbase;
  Q = qg / sbase;

  Ir = (P * Vrterm + Q * Viterm) / (Vterm * Vterm);
  Ii = (P * Viterm - Q * Vrterm) / (Vterm * Vterm);
  x2w = 0;
  x1d = atan2(Viterm + Ir * Xq + Ii * Ra, Vrterm + Ir * Ra - Ii * Xq);
  Id = Ir * sin(x1d) - Ii * cos(x1d); // convert values to the dq axis
  Iq = Ir * cos(x1d) + Ii * sin(x1d); // convert values to the dq axis
  double Vr = Vrterm + Ra * Ir - Xdpp * Ii; // internal voltage on network reference
  double Vi = Viterm + Ra * Ii + Xdpp * Ir; // internal voltage on network reference
  double Vd = Vr * sin(x1d) - Vi * cos(x1d); // convert to dq reference
  double Vq = Vr * cos(x1d) + Vi * sin(x1d); // convert to dq reference
  double Psiqpp = -Vd / (1 + x2w);
  double Psidpp = +Vq / (1 + x2w);

  //printf("Psiqpp=%f\n", Psiqpp);
  x4Psidp = Psidpp - Id * (Xdpp - Xl);
  x3Eqp = x4Psidp + Id * (Xdp - Xl);
  x5Psiqp = (Xdpp - Xq) * Iq;
  Efd = x3Eqp + Id * (Xd - Xdp); 
  LadIfd = Efd;
  Pmech = Psidpp * Iq - Psiqpp * Id;

  values[0] = x1d;
  values[1] = x2w;
  values[2] = x3Eqp;
  values[3] = x4Psidp;
  values[4] = x5Psiqp;
  
  // Initialize exciters
  p_hasExciter = getphasExciter();
  if (p_hasExciter) {
    p_exciter = getExciter();
    p_exciter->setInitialFieldVoltage(Efd);
    //p_exciter->setTimestep(0.0); // SJin: to be read from input file
    //p_exciter->setTimeincrement(0.01); // SJin: to be read from input file
    p_exciter->init(values);
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
  *nvar = 5;
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
    x1d = real(values[0]);
    x2w = real(values[1]);
    x3Eqp = real(values[2]);
    x4Psidp = real(values[3]);
    x5Psiqp = real(values[4]);
  } else if(p_mode == XDOTVECTOBUS) {
    dx1d = real(values[0]);
    dx2w = real(values[1]);
    dx3Eqp = real(values[2]);
    dx4Psidp = real(values[3]);
    dx5Psiqp = real(values[4]);
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
  double Vterm = sqrt(VD*VD + VQ*VQ);
  int x1d_idx = 0;
  int x2w_idx = 1;
  int x3Eqp_idx = 2;
  int x4Psidp_idx = 3;
  int x5Psiqp_idx = 4;
  // On fault (p_mode == FAULT_EVAL flag), the generator variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    values[x1d_idx] = 0.0;
    values[x2w_idx] = 0.0;
    values[x3Eqp_idx] = 0.0;
    values[x4Psidp_idx] = 0.0;
    values[x5Psiqp_idx] = 0.0;
  } else if(p_mode == RESIDUAL_EVAL) {
    if (p_hasExciter) {
      p_exciter = getExciter();
      Efd = p_exciter->getFieldVoltage(); // Efd are called from Exciter
    }
    //    if (p_hasGovernor) {
    //      p_governor = getGovernor();
    //      Pmech = p_governor->getMechanicalPower();
    //    }

    // Generator equations
    double Psidpp = + x3Eqp * (Xdpp - Xl) / (Xdp - Xl) + x4Psidp * (Xdp - Xdpp) / (Xdp - Xl);
  
    double Telec = Psidpp * Iq - Psiqpp * Id;
    double TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl)) * (-x4Psidp - (Xdp - Xl) * Id + x3Eqp);
    //LadIfd = x3Eqp * (1 + Sat(x3Eqp)) + (Xd - Xdp) * (Id + TempD); // update Ifd later
    LadIfd = x3Eqp + (Xd - Xdp) * (Id + TempD); // update Ifd later

    // CHECK EQUATIONS !!!! 
    // RESIDUAL_EVAL for state 1 to 5
    values[x1d_idx] = x2w * OMEGA_S - dx1d;
    values[x2w_idx] = 1 / (2 * H) * ((Pmech - D * x2w) / (1 + x2w) - Telec) - dx2w; // Pmech can be called from Governor
    //printf("Efd = %f, LadIfd = %f\n", Efd, LadIfd);
    values[x3Eqp_idx] = (Efd - LadIfd) / Tdop - dx3Eqp; 
    values[x4Psidp_idx] = (-x4Psidp - (Xdp - Xl) * Id + x3Eqp) / Tdopp; //Zakaria: I removed (-dx4Psidp) from the end of the equation;
    values[x5Psiqp_idx] = (-x5Psiqp - (Xq - Xdpp) * Iq ) / Tqopp;//Zakaria: I change this equation: (-x5Psiqp + (Xqp - Xl) * Iq + x6Edp) / Tqopp - dx5Psiqp, with a one similiar to that of genrou.cpp in the dynamic full application since x6Edp is not defined here.

  }
  
  return true;
}

/**
 * Return the generator current injection (in rectangular form) 
 * @param [output] IGD - real part of the generator current // SJin: match to Ir
 * @param [output] IGQ - imaginary part of the generator current // SJin: match to Ii 
*/
void GensalGen::getCurrent(double *IGD, double *IGQ)
{
  double Vrterm = VD;
  double Viterm = VQ;

  // Setup
  double Psiqpp = x5Psiqp; 
  double Psidpp = + x3Eqp * (Xdpp - Xl) / (Xdp - Xl) + x4Psidp * (Xdp - Xdpp) / (Xdp - Xl); 
  double Vd = -Psiqpp * (1 + x2w);
  double Vq = +Psidpp * (1 + x2w);

  double Vdterm = Vrterm * sin(x1d) - Viterm * cos(x1d);
  double Vqterm = Vrterm * cos(x1d) + Viterm * sin(x1d);
  // DQ Axis
  Id = (Vd - Vdterm) * G - (Vq - Vqterm) * B;
  Iq = (Vd - Vdterm) * B + (Vq - Vqterm) * G;
  // Generator current injections in the network
  Ir = + Id * sin(x1d) + Iq * cos(x1d);
  Ii = - Id * cos(x1d) + Iq * sin(x1d);

  *IGD = Ir;
  *IGQ = Ii;
}

/**
 * Return the matrix entries
 * @param [output] nval - number of values set
 * @param [output] row - row indics for matrix entries
 * @param [output] col - col indices for matrix entries
 * @param [output] values - matrix entries
 * return true when matrix entries set
 */


bool GensalGen::matrixDiagEntries(int *nval,int *row, int *col, gridpack::ComplexType *values)
{
  // Not implemented yet
  int idx = 0;
  if(p_mode == FAULT_EVAL) {
    *nval = idx;
  } else if(p_mode == DIG_DV) {

    *nval = idx;
  } else if(p_mode == DFG_DV) {

    *nval = idx;
  } else if(p_mode == DIG_DX) {
 
    *nval = idx;
  }
  return true;
}

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
