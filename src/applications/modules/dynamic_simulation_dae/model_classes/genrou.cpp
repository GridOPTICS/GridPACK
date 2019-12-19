/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   genrou.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   10/19/19
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
  /*p_delta = 0.0;
  p_dw    = 0.0;
  p_deltadot = 0.0;
  p_dwdot    = 0.0;
  p_Rs    = 0.0;
  p_H     = 0.0;
  p_D     = 0.0;
  p_Ep    = 0.0;
  p_Pm    = 0.0;
  p_Xdp   = 0.0;*/
  x1d = 0.0; // Rotor angle
  x2w = 0.0; // Rotor speed
  x3Eqp = 0.0; // Transient Q axis Eq 
  x4Psidp = 0.0; // Transient D axis flux
  x5Psiqp = 0.0; // Transient Q axis flux
  x6Edp = 0.0; // Transient D axis Ed
  dx1d = 0.0;
  dx2w = 0.0;
  dx3Eqp = 0.0;
  dx4Psidp = 0.0;
  dx5Psiqp = 0.0;
  dx6Edp = 0.0;
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

  Vd = 0.0;
  Vq = 0.0;

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
  /*data->getValue(GENERATOR_RESISTANCE,&p_Rs,idx);
  data->getValue(GENERATOR_TRANSIENT_REACTANCE,&p_Xdp,idx);
  data->getValue(GENERATOR_INERTIA_CONSTANT_H,&p_H,idx);
  data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&p_D,idx);*/

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

  //printf("before: H=%f,D=%f,Ra=%f,Xd=%f,Xq=%f,Xdp=%f,Xdpp=%f,Xl=%f,Tdop=%f,Tdopp=%f,Tqopp=%f,S10=%f,S12=%f,Xqp=%f,Xqpp=%f,Tqop=%f\n", H,D,Ra,Xd,Xq,Xdp,Xdpp,Xl,Tdop,Tdopp,Tqopp,S10,S12,Xqp,Xqpp,Tqop);

  printf("mbase=%f,sbase=%f\n",mbase,sbase);
  // Convert generator parameters from machine base to MVA base
  /*H *= mbase/sbase;
  D *= mbase/sbase;
  Xdp /= mbase/sbase;

  Xd /= mbase/sbase;
  Xq /= mbase/sbase;
  Xdpp /= mbase/sbase;
  Xl /= mbase/sbase;
  Xqp /= mbase/sbase;
  Xqpp /= mbase/sbase;*/

  //printf("after: H=%f,D=%f,Ra=%f,Xd=%f,Xq=%f,Xdp=%f,Xdpp=%f,Xl=%f,Tdop=%f,Tdopp=%f,Tqopp=%f,S10=%f,S12=%f,Xqp=%f,Xqpp=%f,Tqop=%f\n", H,D,Ra,Xd,Xq,Xdp,Xdpp,Xl,Tdop,Tdopp,Tqopp,S10,S12,Xqp,Xqpp,Tqop);
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
    double B = S10 / ((1.0 - A) * (1.0 - A));
    double result = B * (x - A) * (x - A) / x;
    //printf("a = %f, b = %f, c = %f, A = %f, B = %f, S12 = %f, S10 = %f\n", a_, b_, c_, A, B, S12, S10);
    //printf("Sat result = %f\n", result); 
    //    return result; // Scaled Quadratic with 1.7.1 equations
    return 0.0;
}

/**
 * Initialize generator model before calculation
 * @param [output] values - array where initialized generator variables should be set
 */
void GenrouGen::init(gridpack::ComplexType* values) 
{
  double Vm = sqrt(VD*VD + VQ*VQ); // SJin: voltage VD and VQ come from base_gen_model.hpp
  double Theta = atan2(VQ, VD); // SJin: ang = arctang(VQ/VD); ?
  double P, Q; // Generator real and reactive power
  P = pg / mbase;
  Q = qg / mbase;
  //printf("Vm = %f, Theta = %f, P = %f, Q = %f\n", Vm, Theta, P, Q);
  //printf("Vm=%f, Theta=%f, VD=%f, VQ=%f\n", Vm, Theta, VD, VQ);
  Ir = (P * VD + Q * VQ) / (Vm * Vm);
  Ii = (P * VQ - Q * VD) / (Vm * Vm);
  //printf("111 Ir = %f, Ii = %f\n", Ir, Ii);
  //printf("P=%f,VD=%f,Q=%f,VQ=%f,Vm=%f\n",P,VD,Q,VQ,Vm);
  //printf("Ir = %f, Ii = %f\n", Ir, Ii);
  x2w = 0;
  x1d = atan2(VQ + Ir * Xq + Ii * Ra, VD + Ir * Ra - Ii * Xq);
  Id = Ir * sin(x1d) - Ii * cos(x1d); // convert values to the dq axis
  //printf("Ir = %f, Ii=%f, Iq = %f, Xq=%f,Xqp=%f\n", Ir, Ii, Iq,Xq,Xqp);
  Iq = Ir * cos(x1d) + Ii * sin(x1d); // convert values to the dq axis
  double Vr = VD + Ra * Ir - Xdpp * Ii; // internal voltage on network reference
  double Vi = VQ + Ra * Ii + Xdpp * Ir; // internal voltage on network reference
  //double Vd = Vr * sin(x1d) - Vi * sin(x1d); // convert to dq reference
  //double Vq = Vr * cos(x1d) + Vi * cos(x1d); // convert to dq reference
  Vd = Vr * sin(x1d) - Vi * cos(x1d); // convert to dq reference
  Vq = Vr * cos(x1d) + Vi * sin(x1d); // convert to dq reference
  double Psiqpp = -Vd / (1 + x2w);
  double Psidpp = +Vq / (1 + x2w);
  //printf("Psiqpp=%f\n", Psiqpp);
  x4Psidp = Psidpp - Id * (Xdpp - Xl);
  x3Eqp = x4Psidp + Id * (Xdp - Xl);
  x6Edp = Iq * (Xq - Xqp);
  x5Psiqp = x6Edp + Iq * (Xqp - Xl);
  //Efd = x3Eqp * (1 + Sat(x3Eqp)) + Id * (Xd - Xdp); 
  //Efd = x3Eqp + Id * (Xd - Xdp);
  Efd = x3Eqp * (1 + Sat(x3Eqp)) + Id * (Xd - Xdp);
  LadIfd = Efd;
  Pmech = Psidpp * Iq - Psiqpp * Id;

  values[0] = x1d;
  values[1] = x2w;
  values[2] = x3Eqp;
  values[3] = x4Psidp;
  values[4] = x5Psiqp;
  values[5] = x6Edp;

  //printf("VD=%f, VQ=%f, Vm=%f, mag=%f, ang=%f, Pmw=%f, Qmvar=%f, Ir=%f, Ii=%f, Id=%f, Iq=%f, I=%f\n", VD, VQ, Vm, mag, ang, pg, qg, Ir, Ii, Id, Iq, sqrt(Id*Id+Iq*Iq));
  printf("\ngenrou init: x1d = %f, x2w = %f, x3Eqp = %f, x4Psidp = %f, x5Psiqp = %f, x6Edp = %f\n", x1d, x2w, x3Eqp, x4Psidp, x5Psiqp, x6Edp);
  //printf("genrou init: Efd = %f, LadIfd = %f, Pmech = %f\n", Efd, LadIfd, Pmech);
  //printf("VD = %f, VQ=%f\n", VD, VQ);
  
  // Initialize exciters
  p_hasExciter = getphasExciter();
  if (p_hasExciter) {
    p_exciter = getExciter();
    p_exciter->setInitialFieldVoltage(Efd);
    p_exciter->setFieldCurrent(LadIfd);
    p_exciter->setInitialTimeStep(0.01); // SJin: to be read from input file
  }
    
  p_hasGovernor = getphasGovernor();
  if (p_hasGovernor) {
    p_governor = getGovernor();
    p_governor->setInitialMechanicalPower(Pmech);
    p_governor->setRotorSpeedDeviation(x2w);
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
  return x1d;
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
          p_bus_id, p_ckt.c_str(), x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqp_1, x6Edp_1);
  }*/
  //printf("...........bid=%d: %f\t%f\t%f\t%f\t%f\t%f\n", bid, x1d, x2w, x3Eqp, x4Psidp, x5Psiqp, x6Edp);
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
    x1d = real(values[0]);
    x2w = real(values[1]);
    x3Eqp = real(values[2]);
    x4Psidp = real(values[3]);
    x5Psiqp = real(values[4]);
    x6Edp = real(values[5]);
  } else if(p_mode == XDOTVECTOBUS) {
    dx1d = real(values[0]);
    dx2w = real(values[1]);
    dx3Eqp = real(values[2]);
    dx4Psidp = real(values[3]);
    dx5Psiqp = real(values[4]);
    dx6Edp = real(values[5]);
  } else if(p_mode == XVECPRETOBUS) {
    /* This state is only called once before fault on/off condition */
    x1dprev = real(values[0]);
    x2wprev = real(values[1]);
    x3Eqpprev = real(values[2]);
    x4Psidpprev = real(values[3]);
    x5Psiqpprev = real(values[4]);
    x6Edpprev = real(values[5]);
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
  int x1d_idx = 0;
  int x2w_idx = 1;
  int x3Eqp_idx = 2;
  int x4Psidp_idx = 3;
  int x5Psiqp_idx = 4;
  int x6Edp_idx = 5;
  // On fault (p_mode == FAULT_EVAL flag), the generator variables are held constant. This is done by setting the vector values of residual function to 0.0.
  if(p_mode == FAULT_EVAL) {
    values[x1d_idx] = x1d - x1dprev;
    values[x2w_idx] = x2w - x2wprev;
    values[x3Eqp_idx] = x3Eqp - x3Eqpprev;
    values[x4Psidp_idx] = x4Psidp - x4Psidpprev;
    values[x5Psiqp_idx] = x5Psiqp - x5Psiqpprev;
    values[x6Edp_idx] = x6Edp - x6Edpprev;

    if (p_hasGovernor) {
      p_governor->setRotorSpeedDeviation(x2w);
    }
  } else if(p_mode == RESIDUAL_EVAL) {

    if (p_hasExciter) {
      p_exciter = getExciter();
      Efd = p_exciter->getFieldVoltage(); // Efd are called from Exciter
      //printf("\nEfd from exciter: %f\n", Efd);
    }
    
    if (p_hasGovernor) {
      p_governor = getGovernor();
      Pmech = p_governor->getMechanicalPower();
    }

    // Generator equations
    double Psiqpp = - x6Edp * (Xqpp - Xl) / (Xqp - Xl) - x5Psiqp * (Xqp - Xqpp) / (Xqp - Xl);
    double Psidpp = + x3Eqp * (Xdpp - Xl) / (Xdp - Xl) + x4Psidp * (Xdp - Xdpp) / (Xdp - Xl);
  
    double Telec = Psidpp * Iq - Psiqpp * Id;
    double TempD = (Xdp - Xdpp) / ((Xdp - Xl) * (Xdp - Xl)) * (-x4Psidp - (Xdp - Xl) * Id + x3Eqp);
    LadIfd = x3Eqp * (1 + Sat(x3Eqp)) + (Xd - Xdp) * (Id + TempD); // update Ifd later
    //LadIfd = x3Eqp + (Xd - Xdp) * (Id + TempD); // update Ifd later
   //printf("Psiqpp=%f,Psidpp=%f,Telec=%f,TempD=%f,LadIfd=%f\n",Psiqpp,Psidpp,Telec,TempD,LadIfd); 
   //printf("Id=%f, Iq=%f\n", Id, Iq);

    // RESIDUAL_EVAL for state 1 to 6
    values[x1d_idx] = x2w * OMEGA_S - dx1d;
    values[x2w_idx] = 1 / (2 * H) * ((Pmech - D * x2w) / (1 + x2w) - Telec) - dx2w; // Pmech can be called from Governor
    //printf("Efd = %f, LadIfd = %f\n", Efd, LadIfd);
    values[x3Eqp_idx] = (Efd - LadIfd) / Tdop - dx3Eqp; 
    values[x4Psidp_idx] = (-x4Psidp - (Xdp - Xl) * Id + x3Eqp) / Tdopp - dx4Psidp;
    values[x5Psiqp_idx] = (-x5Psiqp + (Xqp - Xl) * Iq + x6Edp) / Tqopp - dx5Psiqp;
    double TempQ = (Xqp - Xqpp) / ((Xqp - Xl) * (Xqp - Xl)) * (-x5Psiqp + (Xqp - Xl) * Iq + x6Edp);
    values[x6Edp_idx] = (-x6Edp + (Xq - Xqp) * (Iq - TempQ)) / Tqop - dx6Edp; 

    if (p_hasExciter) {
      p_exciter->setFieldCurrent(LadIfd);
    }
      
    if (p_hasGovernor) {
      p_governor->setRotorSpeedDeviation(x2w);
    }
  }
  
  return true;
}

/**
 * Return the generator current injection (in rectangular form) 
 * @param [output] IGD - real part of the generator current // SJin: match to Ir
 * @param [output] IGQ - imaginary part of the generator current // SJin: match to Ii 
*/
void GenrouGen::getCurrent(double *IGD, double *IGQ)
{
  double Vm = VD*VD + VQ*VQ;
  // Admittance
  B = -Xdpp / (Ra * Ra + Xdpp * Xdpp);
  G = Ra / (Ra * Ra + Xdpp * Xdpp);

  // Setup
  double Psiqpp = - x6Edp * (Xqpp - Xl) / (Xqp - Xl) - x5Psiqp * (Xqp - Xqpp) / (Xqp - Xl); 
  double Psidpp = + x3Eqp * (Xdpp - Xl) / (Xdp - Xl) + x4Psidp * (Xdp - Xdpp) / (Xdp - Xl); 
  Vd = -Psiqpp * (1 + x2w);
  Vq = +Psidpp * (1 + x2w);

  double Vdterm = VD * sin(x1d) - VQ * cos(x1d);
  double Vqterm = VD * cos(x1d) + VQ * sin(x1d);

  // dq Axis
  Id = (Vd - Vdterm) * G - (Vq - Vqterm) * B;
  Iq = (Vd - Vdterm) * B + (Vq - Vqterm) * G;

  // Generator current injections in the network
  Ir = + Id * sin(x1d) + Iq * cos(x1d);
  Ii = - Id * cos(x1d) + Iq * sin(x1d);

  *IGD = Ir; // SJin: To be confirmed
  *IGQ = Ii; // SJin: To be confirmed

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
    //printf("B = %f, G = %f, VD = %f, VQ = %f, x1d = %f, sin(x1d) = %f, cos(x1d) = %f\n", B, G, VD, VQ, x1d, sin(x1d), cos(x1d));
    // These are the partial derivatives of the generator currents (see getCurrent function) w.r.t to the voltage variables VD and VQ    
    /*row[idx] = 0; col[idx] = 0;
    values[idx] = sin(x1d)*(B*((Vd*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) + 1) - G*((Vd*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq))) - cos(x1d)*(B*((Vd*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vd*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) + 1));
    //values[idx] = -cos(x1d);//B*cos(x1d) + G*sin(x1d);
    //values[idx] = B*sin(x1d) - G*cos(x1d);
    //printf("values[%d] = %f\n", idx, values[idx]);
    idx++;
    row[idx] = 0; col[idx] = 1; 
    values[idx] = - cos(x1d)*(B*((Vq*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - 1) - G*((Vq*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq))) - sin(x1d)*(B*((Vq*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vq*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - 1));
    //values[idx] = sin(x1d); //G*cos(x1d) - B*sin(x1d);
    //values[idx] = B*cos(x1d) + G*sin(x1d);
    //printf("values[%d] = %f\n", idx, values[idx]);
    idx++;
    row[idx] = 1; col[idx] = 0;
    values[idx] = cos(x1d)*(B*((Vd*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) + 1) - G*((Vd*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq))) + sin(x1d)*(B*((Vd*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vd*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) + 1));
    //values[idx] = sin(x1d); //B*sin(x1d) - G*cos(x1d);
    //values[idx] = B*cos(x1d) + G*sin(x1d);
    //printf("values[%d] = %f\n", idx, values[idx]);
    idx++;
    row[idx] = 1; col[idx] = 1;
    values[idx] = sin(x1d)*(B*((Vq*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - 1) - G*((Vq*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq))) - cos(x1d)*(B*((Vq*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vq*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - 1));
    //values[idx] = cos(x1d); //B*cos(x1d) + G*sin(x1d);
    //values[idx] = G*cos(x1d) - B*sin(x1d);
    //printf("values[%d] = %f\n", idx, values[idx]);
    idx++;*/

    *nval = idx;
  } else if(p_mode == DFG_DV) {  // SJin: Jacobian matrix block Jfyi 
    // These are the partial derivatives of the generator equations w.r.t variables VD and VQ
    /*row[idx] = 1; col[idx] = 0;
    values[idx] = -((B*((Vd*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) + 1) - G*((Vd*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq)))*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl)) + (B*((Vd*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vd*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) + 1))*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp)))/(2*H);
    idx++;
    row[idx] = 1; col[idx] = 1;
    values[idx] = ((B*((Vq*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vq*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - 1))*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl)) - (B*((Vq*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - 1) - G*((Vq*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq)))*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp)))/(2*H);
    idx++;
    row[idx] = 2; col[idx] = 0;
    values[idx] = -((Xd - Xdp)*(B*((Vd*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vd*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) + 1) - ((Xdp - Xdpp)*(B*((Vd*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vd*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) + 1)))/(Xdp - Xl)))/Tdop;
    idx++;
    row[idx] = 2; col[idx] = 1;
    values[idx] = ((Xd - Xdp)*(G*((Vq*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq)) - B*((Vq*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - 1) + ((Xdp - Xdpp)*(B*((Vq*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - 1) - G*((Vq*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq))))/(Xdp - Xl)))/Tdop;
    idx++;
    row[idx] = 3; col[idx] = 0;
    values[idx] = -((Xdp - Xl)*(B*((Vd*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vd*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) + 1)))/Tdopp;
    idx++;
    row[idx] = 3; col[idx] = 1;
    values[idx] = -((Xdp - Xl)*(B*((Vq*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - 1) - G*((Vq*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq))))/Tdopp;
    idx++;
    row[idx] = 4; col[idx] = 0;
    values[idx] = -((Xl - Xqp)*(B*((Vd*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) + 1) - G*((Vd*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq))))/Tqopp;
    idx++;
    row[idx] = 4; col[idx] = 1;
    values[idx] = ((Xl - Xqp)*(B*((Vq*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vq*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - 1)))/Tqopp;
    idx++;
    row[idx] = 5; col[idx] = 0;
    values[idx] = ((Xq - Xqp)*(B*((Vd*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) + 1) - G*((Vd*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq)) + ((Xqp - Xqpp)*(B*((Vd*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vd*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) + 1) - G*((Vd*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vd*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq))))/(Xl - Xqp)))/Tqop;
    idx++;
    row[idx] = 5; col[idx] = 1;
    values[idx] = -((Xq - Xqp)*(B*((Vq*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vq*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - 1) + ((Xqp - Xqpp)*(B*((Vq*cos(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - (Vq*sin(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq)) + G*((Vq*cos(Theta)*cos(x1d))/sqrt(Vd*Vd+Vq*Vq) + (Vq*sin(Theta)*sin(x1d))/sqrt(Vd*Vd+Vq*Vq) - 1)))/(Xl - Xqp)))/Tqop;
    idx++;*/

    *nval = idx;
  } else if(p_mode == DIG_DX) { // SJin: Jacobian matrix block Jgx (CORRECT! Output is consistent with Finite difference Jacobian)
    // These are the partial derivatives of the generator currents (see getCurrent) w.r.t generator variables
    row[idx] = 0; col[idx] = 0;
    values[idx] = cos(x1d)*(B*(VQ*cos(x1d) - VD*sin(x1d) + ((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(x2w + 1)) - G*(VD*cos(x1d) + VQ*sin(x1d) - ((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(x2w + 1))) - sin(x1d)*(B*(VD*cos(x1d) + VQ*sin(x1d)) + G*(VQ*cos(x1d) - VD*sin(x1d))) - cos(x1d)*(B*(VQ*cos(x1d) - VD*sin(x1d)) - G*(VD*cos(x1d) + VQ*sin(x1d))) + sin(x1d)*(B*(VD*cos(x1d) + VQ*sin(x1d) - ((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(x2w + 1)) + G*(VQ*cos(x1d) - VD*sin(x1d) + ((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(x2w + 1))); 
    idx++;
    row[idx] = 0; col[idx] = 1;
    values[idx] = sin(x1d)*(B*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) + G*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))) + cos(x1d)*(B*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl)) - G*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))); 
    idx++;
    row[idx] = 0; col[idx] = 2;
    values[idx] = (B*cos(x1d)*(Xdpp - Xl)*(x2w + 1))/(Xdp - Xl) + (G*sin(x1d)*(Xdpp - Xl)*(x2w + 1))/(Xdp - Xl); 
    idx++;
    row[idx] = 0; col[idx] = 3;
    values[idx] = (B*cos(x1d)*(Xdp - Xdpp)*(x2w + 1))/(Xdp - Xl) + (G*sin(x1d)*(Xdp - Xdpp)*(x2w + 1))/(Xdp - Xl); 
    idx++;
    row[idx] = 0; col[idx] = 4;
    values[idx] = (G*cos(x1d)*(Xqp - Xqpp)*(x2w + 1))/(Xl - Xqp) - (B*sin(x1d)*(Xqp - Xqpp)*(x2w + 1))/(Xl - Xqp); 
    idx++;
    row[idx] = 0; col[idx] = 5;
    values[idx] = (B*sin(x1d)*(Xl - Xqpp)*(x2w + 1))/(Xl - Xqp) - (G*cos(x1d)*(Xl - Xqpp)*(x2w + 1))/(Xl - Xqp); 
    idx++;

    row[idx] = 1; col[idx] = 0;
    values[idx] = sin(x1d)*(B*(VQ*cos(x1d) - VD*sin(x1d)) - G*(VD*cos(x1d) + VQ*sin(x1d))) - cos(x1d)*(B*(VD*cos(x1d) + VQ*sin(x1d)) + G*(VQ*cos(x1d) - VD*sin(x1d))) + cos(x1d)*(B*(VD*cos(x1d) + VQ*sin(x1d) - ((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(x2w + 1)) + G*(VQ*cos(x1d) - VD*sin(x1d) + ((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(x2w + 1))) - sin(x1d)*(B*(VQ*cos(x1d) - VD*sin(x1d) + ((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(x2w + 1)) - G*(VD*cos(x1d) + VQ*sin(x1d) - ((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(x2w + 1))); 
    idx++;
    row[idx] = 1; col[idx] = 1;
    values[idx] = cos(x1d)*(B*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) + G*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))) - sin(x1d)*(B*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl)) - G*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))); 
    idx++;
    row[idx] = 1; col[idx] = 2;
    values[idx] = (G*cos(x1d)*(Xdpp - Xl)*(x2w + 1))/(Xdp - Xl) - (B*sin(x1d)*(Xdpp - Xl)*(x2w + 1))/(Xdp - Xl); 
    idx++;
    row[idx] = 1; col[idx] = 3;
    values[idx] = (G*cos(x1d)*(Xdp - Xdpp)*(x2w + 1))/(Xdp - Xl) - (B*sin(x1d)*(Xdp - Xdpp)*(x2w + 1))/(Xdp - Xl); 
    idx++;
    row[idx] = 1; col[idx] = 4;
    values[idx] = - (B*cos(x1d)*(Xqp - Xqpp)*(x2w + 1))/(Xl - Xqp) - (G*sin(x1d)*(Xqp - Xqpp)*(x2w + 1))/(Xl - Xqp); 
    idx++;
    row[idx] = 1; col[idx] = 5;
    values[idx] = (B*cos(x1d)*(Xl - Xqpp)*(x2w + 1))/(Xl - Xqp) + (G*sin(x1d)*(Xl - Xqpp)*(x2w + 1))/(Xl - Xqp); 
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
    values[idx] = (((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(B*(VD*cos(x1d) + VQ*sin(x1d)) + G*(VQ*cos(x1d) - VD*sin(x1d))) - ((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(B*(VQ*cos(x1d) - VD*sin(x1d)) - G*(VD*cos(x1d) + VQ*sin(x1d))))/(2*H);
    idx++;
    row[idx] = 1; col[idx] = 1; // row3 col3 for gen1
    values[idx] = -shift -((Pmech - D*x2w)/pow(x2w + 1, 2) + ((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(B*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) + G*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))) - ((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(B*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl)) - G*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))) + D/(x2w + 1))/(2*H);
    idx++;
    row[idx] = 1; col[idx] = 2;
    values[idx] = -(((B*(VQ*cos(x1d) - VD*sin(x1d) + ((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(x2w + 1)) - G*(VD*cos(x1d) + VQ*sin(x1d) - ((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(x2w + 1)))*(Xdpp - Xl))/(Xdp - Xl) - (B*(Xdpp - Xl)*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(x2w + 1))/(Xdp - Xl) + (G*(Xdpp - Xl)*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(x2w + 1))/(Xdp - Xl))/(2*H);
    idx++;
    row[idx] = 1; col[idx] = 3;
    values[idx] = -(((B*(VQ*cos(x1d) - VD*sin(x1d) + ((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(x2w + 1)) - G*(VD*cos(x1d) + VQ*sin(x1d) - ((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(x2w + 1)))*(Xdp - Xdpp))/(Xdp - Xl) - (B*(Xdp - Xdpp)*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(x2w + 1))/(Xdp - Xl) + (G*(Xdp - Xdpp)*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(x2w + 1))/(Xdp - Xl))/(2*H);
    idx++;
    row[idx] = 1; col[idx] = 4;
    values[idx] = (((B*(VD*cos(x1d) + VQ*sin(x1d) - ((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(x2w + 1)) + G*(VQ*cos(x1d) - VD*sin(x1d) + ((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(x2w + 1)))*(Xqp - Xqpp))/(Xl - Xqp) + (B*(Xqp - Xqpp)*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(x2w + 1))/(Xl - Xqp) + (G*(Xqp - Xqpp)*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(x2w + 1))/(Xl - Xqp))/(2*H);
    idx++;
    row[idx] = 1; col[idx] = 5;
    values[idx] = -(((B*(VD*cos(x1d) + VQ*sin(x1d) - ((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(x2w + 1)) + G*(VQ*cos(x1d) - VD*sin(x1d) + ((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(x2w + 1)))*(Xl - Xqpp))/(Xl - Xqp) + (B*(Xl - Xqpp)*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))*(x2w + 1))/(Xl - Xqp) + (G*(Xl - Xqpp)*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))*(x2w + 1))/(Xl - Xqp))/(2*H);
    idx++;

    row[idx] = 2; col[idx] = 0;
    values[idx] = ((Xd - Xdp)*(G*(VD*cos(x1d) + VQ*sin(x1d)) - B*(VQ*cos(x1d) - VD*sin(x1d)) + ((Xdp - Xdpp)*(B*(VQ*cos(x1d) - VD*sin(x1d)) - G*(VD*cos(x1d) + VQ*sin(x1d))))/(Xdp - Xl)))/Tdop;
    idx++;
    row[idx] = 2; col[idx] = 1;
    values[idx] = -((Xd - Xdp)*(G*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) - B*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl)) + ((Xdp - Xdpp)*(B*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl)) - G*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))))/(Xdp - Xl)))/Tdop;
    idx++;
    row[idx] = 2; col[idx] = 2;
    values[idx] = -shift -((Xd - Xdp)*(((B*(Xdpp - Xl)*(x2w + 1) + 1)*(Xdp - Xdpp))/pow(Xdp - Xl, 2) - (B*(Xdpp - Xl)*(x2w + 1))/(Xdp - Xl)) + 1)/Tdop;
    idx++;
    row[idx] = 2; col[idx] = 3;
    values[idx] = -((Xd - Xdp)*(((B*(Xdp - Xdpp)*(x2w + 1) - 1)*(Xdp - Xdpp))/pow(Xdp - Xl, 2) - (B*(Xdp - Xdpp)*(x2w + 1))/(Xdp - Xl)))/Tdop;
    idx++;
    row[idx] = 2; col[idx] = 4;
    values[idx] = (((G*(Xqp - Xqpp)*(x2w + 1))/(Xl - Xqp) - (G*(Xdp - Xdpp)*(Xqp - Xqpp)*(x2w + 1))/((Xdp - Xl)*(Xl - Xqp)))*(Xd - Xdp))/Tdop;
    idx++;
    row[idx] = 2; col[idx] = 5;
    values[idx] = -(((G*(Xl - Xqpp)*(x2w + 1))/(Xl - Xqp) - (G*(Xdp - Xdpp)*(Xl - Xqpp)*(x2w + 1))/((Xdp - Xl)*(Xl - Xqp)))*(Xd - Xdp))/Tdop;
    idx++;

    row[idx] = 3; col[idx] = 0;
    values[idx] = -((Xdp - Xl)*(B*(VQ*cos(x1d) - VD*sin(x1d)) - G*(VD*cos(x1d) + VQ*sin(x1d))))/Tdopp;
    idx++;
    row[idx] = 3; col[idx] = 1;
    values[idx] = ((Xdp - Xl)*(B*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl)) - G*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp))))/Tdopp;
    idx++;
    row[idx] = 3; col[idx] = 2;
    values[idx] = (B*(Xdpp - Xl)*(x2w + 1) + 1)/Tdopp;
    idx++;
    row[idx] = 3; col[idx] = 3;
    values[idx] = -shift + (B*(Xdp - Xdpp)*(x2w + 1) - 1)/Tdopp;
    idx++;
    row[idx] = 3; col[idx] = 4;
    values[idx] = (G*(Xdp - Xl)*(Xqp - Xqpp)*(x2w + 1))/(Tdopp*(Xl - Xqp));
    idx++;
    row[idx] = 3; col[idx] = 5;
    values[idx] = -(G*(Xdp - Xl)*(Xl - Xqpp)*(x2w + 1))/(Tdopp*(Xl - Xqp));
    idx++;

    row[idx] = 4; col[idx] = 0;
    values[idx] = ((Xl - Xqp)*(B*(VD*cos(x1d) + VQ*sin(x1d)) + G*(VQ*cos(x1d) - VD*sin(x1d))))/Tqopp;
    idx++;
    row[idx] = 4; col[idx] = 1;
    values[idx] = -((Xl - Xqp)*(B*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) + G*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))))/Tqopp;
    idx++;
    row[idx] = 4; col[idx] = 2;
    values[idx] = -(G*(Xdpp - Xl)*(Xl - Xqp)*(x2w + 1))/(Tqopp*(Xdp - Xl));
    idx++;
    row[idx] = 4; col[idx] = 3;
    values[idx] = -(G*(Xdp - Xdpp)*(Xl - Xqp)*(x2w + 1))/(Tqopp*(Xdp - Xl));
    idx++;
    row[idx] = 4; col[idx] = 4;
    values[idx] = -shift + (B*(Xqp - Xqpp)*(x2w + 1) - 1)/Tqopp;
    idx++;
    row[idx] = 4; col[idx] = 5;
    values[idx] = -(B*(Xl - Xqpp)*(x2w + 1) - 1)/Tqopp;
    idx++;

    row[idx] = 5; col[idx] = 0;
    values[idx] = -((Xq - Xqp)*(B*(VD*cos(x1d) + VQ*sin(x1d)) + G*(VQ*cos(x1d) - VD*sin(x1d)) + ((Xqp - Xqpp)*(B*(VD*cos(x1d) + VQ*sin(x1d)) + G*(VQ*cos(x1d) - VD*sin(x1d))))/(Xl - Xqp)))/Tqop;
    idx++;
    row[idx] = 5; col[idx] = 1;
    values[idx] = ((Xq - Xqp)*(B*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) + G*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl)) + ((Xqp - Xqpp)*(B*((x6Edp*(Xl - Xqpp))/(Xl - Xqp) - (x5Psiqp*(Xqp - Xqpp))/(Xl - Xqp)) + G*((x4Psidp*(Xdp - Xdpp))/(Xdp - Xl) + (x3Eqp*(Xdpp - Xl))/(Xdp - Xl))))/(Xl - Xqp)))/Tqop;
    idx++;
    row[idx] = 5; col[idx] = 2;
    values[idx] = (((G*(Xdpp - Xl)*(x2w + 1))/(Xdp - Xl) + (G*(Xdpp - Xl)*(Xqp - Xqpp)*(x2w + 1))/((Xdp - Xl)*(Xl - Xqp)))*(Xq - Xqp))/Tqop;
    idx++;
    row[idx] = 5; col[idx] = 3;
    values[idx] = (((G*(Xdp - Xdpp)*(x2w + 1))/(Xdp - Xl) + (G*(Xdp - Xdpp)*(Xqp - Xqpp)*(x2w + 1))/((Xdp - Xl)*(Xl - Xqp)))*(Xq - Xqp))/Tqop;
    idx++;
    row[idx] = 5; col[idx] = 4;
    values[idx] = -((Xq - Xqp)*(((B*(Xqp - Xqpp)*(x2w + 1) - 1)*(Xqp - Xqpp))/pow(Xl - Xqp, 2) + (B*(Xqp - Xqpp)*(x2w + 1))/(Xl - Xqp)))/Tqop;
    idx++;
    row[idx] = 5; col[idx] = 5;
    values[idx] = -shift + ((Xq - Xqp)*(((B*(Xl - Xqpp)*(x2w + 1) - 1)*(Xqp - Xqpp))/pow(Xl - Xqp, 2) + (B*(Xl - Xqpp)*(x2w + 1))/(Xl - Xqp)) - 1)/Tqop;
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

