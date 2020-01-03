/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   genrou.cpp
 * @author Shuangshuang Jin 
 *
   @author Shrirang Abhyankar
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
  if (!data->getValue(GENERATOR_S1, &S10, idx)) S10=0.17; // S10 TBD: check parser
  if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.55; // S12 TBD: check parser
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
}

/**
 * Saturation function
 * @ param x
 */
double GenrouGen::Sat(double x)
{
    double a_ = S12 / S10 - 1.0 / 1.2;
    double b_ = -2 * S12 / S10 + 2;
    double c_ = S12 / S10 - 1.2;
    double A = (-b_ - sqrt(b_ * b_ - 4 * a_ * c_)) / (2 * a_);
    double B_ = S10 / ((1.0 - A) * (1.0 - A));
    double result = B_ * (x - A) * (x - A) / x;
    return result; // Scaled Quadratic with 1.7.1 equations

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
  P = pg / mbase;
  Q = qg / mbase;

  // Generator currents in network reference frame
  Ir = (P * VD + Q * VQ) / (Vm * Vm);
  Ii = (P * VQ - Q * VD) / (Vm * Vm);

  // Machine angle and speed deviation
  delta = atan2(VQ + Ir * Xq + Ii * Ra, VD + Ir * Ra - Ii * Xq);
  dw = 0.0;

  theta = PI/2.0 - delta; // Rotor angle lags quadrature axis by 90 deg.

  // Generator currents in dq axis reference frame
  Id = Ir * cos(theta) - Ii * sin(theta); // convert values to the dq axis
  Iq = Ir * sin(theta) + Ii * cos(theta); // convert values to the dq axis

  Vd = VD * cos(theta) - VQ * sin(theta); // convert to dq reference
  Vq = VD * sin(theta) + VQ * cos(theta); // convert to dq reference

  Eqp = Vq + Xdp*Id;
  Edp = Vd - Xqp*Iq;

  Psidp =  Eqp - Id * (Xdp - Xl);
  Psiqp = -Edp - Iq * (Xqp - Xl);

  Efd = Eqp * (1 + Sat(Eqp)) + Id * (Xd - Xdp);
  LadIfd = Efd;
  Pmech = P;

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
    p_governor->setRotorSpeedDeviation(dw);
    p_governor->setInitialTimestep(0.01); // SJin: to be read from input file
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

    if (p_hasGovernor) {
      p_governor->setRotorSpeedDeviation(dw);
    }
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

    double Vdterm = VD * sin(delta) - VQ * cos(delta);
    double Vqterm = VD * cos(delta) + VQ * sin(delta);

    double tempd1,tempd2,tempq1,tempq2;
    tempd1 = (Xdpp - Xl)/(Xdp - Xl);
    tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
    tempq1 = (Xdpp - Xl)/(Xqp - Xl);
    tempq2 = (Xqp - Xdpp)/(Xqp - Xl);

    // dq Axis currents
    Id = (tempd1*Eqp + tempd2*Psidp - Vqterm)/Xdpp;
    Iq = (tempq1*Edp - tempq2*Psiqp - Vdterm)/-Xdpp;

    // Electrical torque
    double Telec = tempd1*Eqp*Iq + tempd2*Psidp*Iq + tempq1*Edp*Id - tempq2*Psiqp*Id; 

    double TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl)) * (-Psidp - (Xdp - Xl) * Id + Eqp);

    // Field current 
    LadIfd = Eqp * (1 + Sat(Eqp)) + (Xd - Xdp) * (Id + TempD);

    // RESIDUAL_EVAL for state 1 to 6
    // Machine rotor angle
    values[delta_idx] = dw * OMEGA_S - ddelta;

    // Speed
    values[dw_idx] = 1 / (2 * H) * ((Pmech - D * dw) / (1 + dw) - Telec) - ddw; // Pmech can be called from Governor

    // Q-axis transient EMF
    values[Eqp_idx] = (Efd - LadIfd) / Tdop - dEqp; 

    // D-axis transient flux
    values[Psidp_idx] = (-Psidp - (Xdp - Xl) * Id + Eqp) / Tdopp - dPsidp;
    
    // Q-axis transient flux
    values[Psiqp_idx] = (-Psiqp - (Xqp - Xl) * Iq - Edp) / Tqopp - dPsiqp;

    // D-axis transient EMF
    double TempQ = (Xqp - Xqpp) / ((Xqp - Xl) * (Xqp - Xl)) * (Psiqp + (Xqp - Xl) * Iq + Edp);
    values[Edp_idx] = (-Edp + (Xq - Xqp) * (Iq - TempQ)) / Tqop - dEdp; 

    if (p_hasGovernor) {
      p_governor->setRotorSpeedDeviation(dw);
    }
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
  double Vm = VD*VD + VQ*VQ;

  double Vdterm = VD * sin(delta) - VQ * cos(delta);
  double Vqterm = VD * cos(delta) + VQ * sin(delta);

  double tempd1,tempd2,tempq1,tempq2;
  tempd1 = (Xdpp - Xl)/(Xdp - Xl);
  tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
  tempq1 = (Xdpp - Xl)/(Xqp - Xl);
  tempq2 = (Xqp - Xdpp)/(Xqp - Xl);
  // dq Axis currents
  Id = (tempd1*Eqp + tempd2*Psidp - Vqterm)/Xdpp;
  Iq = (tempq1*Edp - tempq2*Psiqp - Vdterm)/-Xdpp;

  // Generator current injections in the network
  *IGD = + Id * sin(delta) + Iq * cos(delta);
  *IGQ = - Id * cos(delta) + Iq * sin(delta);

}

double GenrouGen::getFieldCurrent()
{
    double Vdterm = VD * sin(delta) - VQ * cos(delta);
    double Vqterm = VD * cos(delta) + VQ * sin(delta);

    double tempd1,tempd2,tempq1,tempq2;
    tempd1 = (Xdpp - Xl)/(Xdp - Xl);
    tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
    tempq1 = (Xdpp - Xl)/(Xqp - Xl);
    tempq2 = (Xqp - Xdpp)/(Xqp - Xl);

    // dq Axis currents
    Id = (tempd1*Eqp + tempd2*Psidp - Vqterm)/Xdpp;
    Iq = (tempq1*Edp - tempq2*Psiqp - Vdterm)/-Xdpp;

    double TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl)) * (-Psidp - (Xdp - Xl) * Id + Eqp);

    // Field current 
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
bool GenrouGen::matrixDiagEntries(int *nval,int *row, int *col, gridpack::ComplexType *values)
{
  int idx = 0;
  if(p_mode == FAULT_EVAL) { // SJin: put values 1 along diagonals, 0 along off diagonals
  // On fault (p_mode == FAULT_EVAL flag), the generator variables are held constant. This is done by setting the diagonal matrix entries to 1.0 and all other entries to 0. The residual function values are already set to 0.0 in the vector values function. This results in the equation 1*dx = 0.0 such that dx = 0.0 and hence x does not get changed.
    row[idx] = 0; col[idx] = 0;
    values[idx] = 1.0;
    idx++;
    row[idx] = 1; col[idx] = 1;
    values[idx] = 1.0;
    idx++;
    row[idx] = 2; col[idx] = 2;
    values[idx] = 1.0;
    idx++;
    row[idx] = 3; col[idx] = 3;
    values[idx] = 1.0;
    idx++;
    row[idx] = 4; col[idx] = 4;
    values[idx] = 1.0;
    idx++;
    row[idx] = 5; col[idx] = 5;
    values[idx] = 1.0;
    idx++;
    *nval = idx;
  } else if(p_mode == DIG_DV) { // SJin: Jacobian matrix block Jgy 
    //printf("B = %f, G = %f, VD = %f, VQ = %f, delta = %f, sin(delta) = %f, cos(delta) = %f\n", B, G, VD, VQ, delta, sin(delta), cos(delta));
    // These are the partial derivatives of the generator currents (see getCurrent function) w.r.t to the voltage variables VD and VQ    
    /*row[idx] = 0; col[idx] = 0;
    values[idx] = sin(delta)*(B*((Vd*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) + 1) - G*((Vd*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq))) - cos(delta)*(B*((Vd*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vd*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) + 1));
    //values[idx] = -cos(delta);//B*cos(delta) + G*sin(delta);
    //values[idx] = B*sin(delta) - G*cos(delta);
    //printf("values[%d] = %f\n", idx, values[idx]);
    idx++;
    row[idx] = 0; col[idx] = 1; 
    values[idx] = - cos(delta)*(B*((Vq*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - 1) - G*((Vq*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq))) - sin(delta)*(B*((Vq*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vq*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - 1));
    //values[idx] = sin(delta); //G*cos(delta) - B*sin(delta);
    //values[idx] = B*cos(delta) + G*sin(delta);
    //printf("values[%d] = %f\n", idx, values[idx]);
    idx++;
    row[idx] = 1; col[idx] = 0;
    values[idx] = cos(delta)*(B*((Vd*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) + 1) - G*((Vd*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq))) + sin(delta)*(B*((Vd*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vd*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) + 1));
    //values[idx] = sin(delta); //B*sin(delta) - G*cos(delta);
    //values[idx] = B*cos(delta) + G*sin(delta);
    //printf("values[%d] = %f\n", idx, values[idx]);
    idx++;
    row[idx] = 1; col[idx] = 1;
    values[idx] = sin(delta)*(B*((Vq*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - 1) - G*((Vq*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq))) - cos(delta)*(B*((Vq*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vq*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - 1));
    //values[idx] = cos(delta); //B*cos(delta) + G*sin(delta);
    //values[idx] = G*cos(delta) - B*sin(delta);
    //printf("values[%d] = %f\n", idx, values[idx]);
    idx++;*/

    *nval = idx;
  } else if(p_mode == DFG_DV) {  // SJin: Jacobian matrix block Jfyi 
    // These are the partial derivatives of the generator equations w.r.t variables VD and VQ
    /*row[idx] = 1; col[idx] = 0;
    values[idx] = -((B*((Vd*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) + 1) - G*((Vd*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq)))*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl)) + (B*((Vd*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vd*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) + 1))*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp)))/(2*H);
    idx++;
    row[idx] = 1; col[idx] = 1;
    values[idx] = ((B*((Vq*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vq*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - 1))*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl)) - (B*((Vq*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - 1) - G*((Vq*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq)))*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp)))/(2*H);
    idx++;
    row[idx] = 2; col[idx] = 0;
    values[idx] = -((Xd - Xdp)*(B*((Vd*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vd*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) + 1) - ((Xdp - Xdpp)*(B*((Vd*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vd*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) + 1)))/(Xdp - Xl)))/Tdop;
    idx++;
    row[idx] = 2; col[idx] = 1;
    values[idx] = ((Xd - Xdp)*(G*((Vq*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq)) - B*((Vq*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - 1) + ((Xdp - Xdpp)*(B*((Vq*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - 1) - G*((Vq*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq))))/(Xdp - Xl)))/Tdop;
    idx++;
    row[idx] = 3; col[idx] = 0;
    values[idx] = -((Xdp - Xl)*(B*((Vd*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vd*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) + 1)))/Tdopp;
    idx++;
    row[idx] = 3; col[idx] = 1;
    values[idx] = -((Xdp - Xl)*(B*((Vq*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - 1) - G*((Vq*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq))))/Tdopp;
    idx++;
    row[idx] = 4; col[idx] = 0;
    values[idx] = -((Xl - Xqp)*(B*((Vd*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) + 1) - G*((Vd*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq))))/Tqopp;
    idx++;
    row[idx] = 4; col[idx] = 1;
    values[idx] = ((Xl - Xqp)*(B*((Vq*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vq*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - 1)))/Tqopp;
    idx++;
    row[idx] = 5; col[idx] = 0;
    values[idx] = ((Xq - Xqp)*(B*((Vd*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) + 1) - G*((Vd*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq)) + ((Xqp - Xqpp)*(B*((Vd*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) + 1) - G*((Vd*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq))))/(Xl - Xqp)))/Tqop;
    idx++;
    row[idx] = 5; col[idx] = 1;
    values[idx] = -((Xq - Xqp)*(B*((Vq*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vq*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - 1) + ((Xqp - Xqpp)*(B*((Vq*cos(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vq*cos(Theta)*cos(delta))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(delta))/sqrt(Vd*Vd+Vq*Vq) - 1)))/(Xl - Xqp)))/Tqop;
    idx++;*/

    *nval = idx;
  } else if(p_mode == DIG_DX) { // SJin: Jacobian matrix block Jgx (CORRECT! Output is consistent with Finite difference Jacobian)
    // These are the partial derivatives of the generator currents (see getCurrent) w.r.t generator variables
    row[idx] = 0; col[idx] = 0;
    values[idx] = cos(delta)*(B*(VQ*cos(delta) - VD*sin(delta) + ((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(dw + 1)) - G*(VD*cos(delta) + VQ*sin(delta) - ((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(dw + 1))) - sin(delta)*(B*(VD*cos(delta) + VQ*sin(delta)) + G*(VQ*cos(delta) - VD*sin(delta))) - cos(delta)*(B*(VQ*cos(delta) - VD*sin(delta)) - G*(VD*cos(delta) + VQ*sin(delta))) + sin(delta)*(B*(VD*cos(delta) + VQ*sin(delta) - ((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(dw + 1)) + G*(VQ*cos(delta) - VD*sin(delta) + ((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(dw + 1))); 
    idx++;
    row[idx] = 0; col[idx] = 1;
    values[idx] = sin(delta)*(B*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) + G*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))) + cos(delta)*(B*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl)) - G*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))); 
    idx++;
    row[idx] = 0; col[idx] = 2;
    values[idx] = (B*cos(delta)*(Xdpp - Xl)*(dw + 1))/(Xdp - Xl) + (G*sin(delta)*(Xdpp - Xl)*(dw + 1))/(Xdp - Xl); 
    idx++;
    row[idx] = 0; col[idx] = 3;
    values[idx] = (B*cos(delta)*(Xdp - Xdpp)*(dw + 1))/(Xdp - Xl) + (G*sin(delta)*(Xdp - Xdpp)*(dw + 1))/(Xdp - Xl); 
    idx++;
    row[idx] = 0; col[idx] = 4;
    values[idx] = (G*cos(delta)*(Xqp - Xqpp)*(dw + 1))/(Xl - Xqp) - (B*sin(delta)*(Xqp - Xqpp)*(dw + 1))/(Xl - Xqp); 
    idx++;
    row[idx] = 0; col[idx] = 5;
    values[idx] = (B*sin(delta)*(Xl - Xqpp)*(dw + 1))/(Xl - Xqp) - (G*cos(delta)*(Xl - Xqpp)*(dw + 1))/(Xl - Xqp); 
    idx++;

    row[idx] = 1; col[idx] = 0;
    values[idx] = sin(delta)*(B*(VQ*cos(delta) - VD*sin(delta)) - G*(VD*cos(delta) + VQ*sin(delta))) - cos(delta)*(B*(VD*cos(delta) + VQ*sin(delta)) + G*(VQ*cos(delta) - VD*sin(delta))) + cos(delta)*(B*(VD*cos(delta) + VQ*sin(delta) - ((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(dw + 1)) + G*(VQ*cos(delta) - VD*sin(delta) + ((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(dw + 1))) - sin(delta)*(B*(VQ*cos(delta) - VD*sin(delta) + ((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(dw + 1)) - G*(VD*cos(delta) + VQ*sin(delta) - ((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(dw + 1))); 
    idx++;
    row[idx] = 1; col[idx] = 1;
    values[idx] = cos(delta)*(B*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) + G*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))) - sin(delta)*(B*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl)) - G*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))); 
    idx++;
    row[idx] = 1; col[idx] = 2;
    values[idx] = (G*cos(delta)*(Xdpp - Xl)*(dw + 1))/(Xdp - Xl) - (B*sin(delta)*(Xdpp - Xl)*(dw + 1))/(Xdp - Xl); 
    idx++;
    row[idx] = 1; col[idx] = 3;
    values[idx] = (G*cos(delta)*(Xdp - Xdpp)*(dw + 1))/(Xdp - Xl) - (B*sin(delta)*(Xdp - Xdpp)*(dw + 1))/(Xdp - Xl); 
    idx++;
    row[idx] = 1; col[idx] = 4;
    values[idx] = - (B*cos(delta)*(Xqp - Xqpp)*(dw + 1))/(Xl - Xqp) - (G*sin(delta)*(Xqp - Xqpp)*(dw + 1))/(Xl - Xqp); 
    idx++;
    row[idx] = 1; col[idx] = 5;
    values[idx] = (B*cos(delta)*(Xl - Xqpp)*(dw + 1))/(Xl - Xqp) + (G*sin(delta)*(Xl - Xqpp)*(dw + 1))/(Xl - Xqp); 
    idx++;
  
    *nval = idx;
  } else { // SJin: Jacobian matrix block Jfxi (CORRECT! Output is consistent with Finite difference Jacobian)
    // Partials of generator equations w.r.t generator variables
    row[idx] = 0; col[idx] = 0; // row2 col2 for gen1
    values[idx] = -shift;
    idx++;
    row[idx] = 0; col[idx] = 1; // row2 col3 for gen1
    values[idx] = OMEGA_S;
    idx++;

    row[idx] = 1; col[idx] = 0;
    values[idx] = (((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(B*(VD*cos(delta) + VQ*sin(delta)) + G*(VQ*cos(delta) - VD*sin(delta))) - ((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(B*(VQ*cos(delta) - VD*sin(delta)) - G*(VD*cos(delta) + VQ*sin(delta))))/(2*H);
    idx++;
    row[idx] = 1; col[idx] = 1; // row3 col3 for gen1
    values[idx] = -shift -((Pmech - D*dw)/pow(dw + 1, 2) + ((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(B*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) + G*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))) - ((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(B*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl)) - G*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))) + D/(dw + 1))/(2*H);
    idx++;
    row[idx] = 1; col[idx] = 2;
    values[idx] = -(((B*(VQ*cos(delta) - VD*sin(delta) + ((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(dw + 1)) - G*(VD*cos(delta) + VQ*sin(delta) - ((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(dw + 1)))*(Xdpp - Xl))/(Xdp - Xl) - (B*(Xdpp - Xl)*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(dw + 1))/(Xdp - Xl) + (G*(Xdpp - Xl)*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(dw + 1))/(Xdp - Xl))/(2*H);
    idx++;
    row[idx] = 1; col[idx] = 3;
    values[idx] = -(((B*(VQ*cos(delta) - VD*sin(delta) + ((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(dw + 1)) - G*(VD*cos(delta) + VQ*sin(delta) - ((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(dw + 1)))*(Xdp - Xdpp))/(Xdp - Xl) - (B*(Xdp - Xdpp)*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(dw + 1))/(Xdp - Xl) + (G*(Xdp - Xdpp)*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(dw + 1))/(Xdp - Xl))/(2*H);
    idx++;
    row[idx] = 1; col[idx] = 4;
    values[idx] = (((B*(VD*cos(delta) + VQ*sin(delta) - ((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(dw + 1)) + G*(VQ*cos(delta) - VD*sin(delta) + ((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(dw + 1)))*(Xqp - Xqpp))/(Xl - Xqp) + (B*(Xqp - Xqpp)*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(dw + 1))/(Xl - Xqp) + (G*(Xqp - Xqpp)*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(dw + 1))/(Xl - Xqp))/(2*H);
    idx++;
    row[idx] = 1; col[idx] = 5;
    values[idx] = -(((B*(VD*cos(delta) + VQ*sin(delta) - ((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(dw + 1)) + G*(VQ*cos(delta) - VD*sin(delta) + ((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(dw + 1)))*(Xl - Xqpp))/(Xl - Xqp) + (B*(Xl - Xqpp)*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))*(dw + 1))/(Xl - Xqp) + (G*(Xl - Xqpp)*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(dw + 1))/(Xl - Xqp))/(2*H);
    idx++;

    row[idx] = 2; col[idx] = 0;
    values[idx] = ((Xd - Xdp)*(G*(VD*cos(delta) + VQ*sin(delta)) - B*(VQ*cos(delta) - VD*sin(delta)) + ((Xdp - Xdpp)*(B*(VQ*cos(delta) - VD*sin(delta)) - G*(VD*cos(delta) + VQ*sin(delta))))/(Xdp - Xl)))/Tdop;
    idx++;
    row[idx] = 2; col[idx] = 1;
    values[idx] = -((Xd - Xdp)*(G*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) - B*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl)) + ((Xdp - Xdpp)*(B*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl)) - G*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))))/(Xdp - Xl)))/Tdop;
    idx++;
    row[idx] = 2; col[idx] = 2;
    values[idx] = -shift -((Xd - Xdp)*(((B*(Xdpp - Xl)*(dw + 1) + 1)*(Xdp - Xdpp))/pow(Xdp - Xl, 2) - (B*(Xdpp - Xl)*(dw + 1))/(Xdp - Xl)) + 1)/Tdop;
    idx++;
    row[idx] = 2; col[idx] = 3;
    values[idx] = -((Xd - Xdp)*(((B*(Xdp - Xdpp)*(dw + 1) - 1)*(Xdp - Xdpp))/pow(Xdp - Xl, 2) - (B*(Xdp - Xdpp)*(dw + 1))/(Xdp - Xl)))/Tdop;
    idx++;
    row[idx] = 2; col[idx] = 4;
    values[idx] = (((G*(Xqp - Xqpp)*(dw + 1))/(Xl - Xqp) - (G*(Xdp - Xdpp)*(Xqp - Xqpp)*(dw + 1))/((Xdp - Xl)*(Xl - Xqp)))*(Xd - Xdp))/Tdop;
    idx++;
    row[idx] = 2; col[idx] = 5;
    values[idx] = -(((G*(Xl - Xqpp)*(dw + 1))/(Xl - Xqp) - (G*(Xdp - Xdpp)*(Xl - Xqpp)*(dw + 1))/((Xdp - Xl)*(Xl - Xqp)))*(Xd - Xdp))/Tdop;
    idx++;

    row[idx] = 3; col[idx] = 0;
    values[idx] = -((Xdp - Xl)*(B*(VQ*cos(delta) - VD*sin(delta)) - G*(VD*cos(delta) + VQ*sin(delta))))/Tdopp;
    idx++;
    row[idx] = 3; col[idx] = 1;
    values[idx] = ((Xdp - Xl)*(B*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl)) - G*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp))))/Tdopp;
    idx++;
    row[idx] = 3; col[idx] = 2;
    values[idx] = (B*(Xdpp - Xl)*(dw + 1) + 1)/Tdopp;
    idx++;
    row[idx] = 3; col[idx] = 3;
    values[idx] = -shift + (B*(Xdp - Xdpp)*(dw + 1) - 1)/Tdopp;
    idx++;
    row[idx] = 3; col[idx] = 4;
    values[idx] = (G*(Xdp - Xl)*(Xqp - Xqpp)*(dw + 1))/(Tdopp*(Xl - Xqp));
    idx++;
    row[idx] = 3; col[idx] = 5;
    values[idx] = -(G*(Xdp - Xl)*(Xl - Xqpp)*(dw + 1))/(Tdopp*(Xl - Xqp));
    idx++;

    row[idx] = 4; col[idx] = 0;
    values[idx] = ((Xl - Xqp)*(B*(VD*cos(delta) + VQ*sin(delta)) + G*(VQ*cos(delta) - VD*sin(delta))))/Tqopp;
    idx++;
    row[idx] = 4; col[idx] = 1;
    values[idx] = -((Xl - Xqp)*(B*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) + G*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))))/Tqopp;
    idx++;
    row[idx] = 4; col[idx] = 2;
    values[idx] = -(G*(Xdpp - Xl)*(Xl - Xqp)*(dw + 1))/(Tqopp*(Xdp - Xl));
    idx++;
    row[idx] = 4; col[idx] = 3;
    values[idx] = -(G*(Xdp - Xdpp)*(Xl - Xqp)*(dw + 1))/(Tqopp*(Xdp - Xl));
    idx++;
    row[idx] = 4; col[idx] = 4;
    values[idx] = -shift + (B*(Xqp - Xqpp)*(dw + 1) - 1)/Tqopp;
    idx++;
    row[idx] = 4; col[idx] = 5;
    values[idx] = -(B*(Xl - Xqpp)*(dw + 1) - 1)/Tqopp;
    idx++;

    row[idx] = 5; col[idx] = 0;
    values[idx] = -((Xq - Xqp)*(B*(VD*cos(delta) + VQ*sin(delta)) + G*(VQ*cos(delta) - VD*sin(delta)) + ((Xqp - Xqpp)*(B*(VD*cos(delta) + VQ*sin(delta)) + G*(VQ*cos(delta) - VD*sin(delta))))/(Xl - Xqp)))/Tqop;
    idx++;
    row[idx] = 5; col[idx] = 1;
    values[idx] = ((Xq - Xqp)*(B*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) + G*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl)) + ((Xqp - Xqpp)*(B*((Edp*(Xl - Xqpp))/(Xl - Xqp) - (Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) + G*((Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (Eqp*(Xdpp - Xl))/(Xdp - Xl))))/(Xl - Xqp)))/Tqop;
    idx++;
    row[idx] = 5; col[idx] = 2;
    values[idx] = (((G*(Xdpp - Xl)*(dw + 1))/(Xdp - Xl) + (G*(Xdpp - Xl)*(Xqp - Xqpp)*(dw + 1))/((Xdp - Xl)*(Xl - Xqp)))*(Xq - Xqp))/Tqop;
    idx++;
    row[idx] = 5; col[idx] = 3;
    values[idx] = (((G*(Xdp - Xdpp)*(dw + 1))/(Xdp - Xl) + (G*(Xdp - Xdpp)*(Xqp - Xqpp)*(dw + 1))/((Xdp - Xl)*(Xl - Xqp)))*(Xq - Xqp))/Tqop;
    idx++;
    row[idx] = 5; col[idx] = 4;
    values[idx] = -((Xq - Xqp)*(((B*(Xqp - Xqpp)*(dw + 1) - 1)*(Xqp - Xqpp))/pow(Xl - Xqp, 2) + (B*(Xqp - Xqpp)*(dw + 1))/(Xl - Xqp)))/Tqop;
    idx++;
    row[idx] = 5; col[idx] = 5;
    values[idx] = -shift + ((Xq - Xqp)*(((B*(Xl - Xqpp)*(dw + 1) - 1)*(Xqp - Xqpp))/pow(Xl - Xqp, 2) + (B*(Xl - Xqpp)*(dw + 1))/(Xl - Xqp)) - 1)/Tqop;
    idx++;
 
    *nval = idx;
  }
  return true;
}

/*void GenrouGen::setExciter(boost::shared_ptr<BaseExcModel> &exciter)
{
  p_exciter = exciter;
}

boost::shared_ptr<BaseExcModel> GenrouGen::getExciter()
{
  return p_exciter;
}*/

