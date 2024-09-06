/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   esst4b.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   Sep 6, 2024
 * 
 * @brief  
 * 
 * 
 */

#include <esst4b.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

#define TS_THRESHOLD 1
#define KI_THRESHOLD 1 // Check?

/**
 *  Basic constructor
 */
Esst4bExc::Esst4bExc(void)
{
  dx1Vm = 0;
  dx2Vcomp = 0;
  dx3Va = 0; 
  dx4Vr = 0; 
  dx1Vm_1 = 0;
  dx2Vcomp_1 = 0;
  dx3Va_1 = 0;
  dx4Vr_1 = 0;
}

/**
 *  Basic destructor
 */
Esst4bExc::~Esst4bExc(void)
{
}

void Esst4bExc::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxexc = 0;
  *nvar = nxexc;
}

void Esst4bExc::preStep(double time, double timestep)
{
  if(integrationtype != EXPLICIT) return;

  // block-based explicit implementation (predictor)

  double vabc[3],vdq0[3];

  vabc[0] = p_va; vabc[1] = p_vb; vabc[2] = p_vc;

  double delta = getGenerator()->getAngle();
  
  abc2dq0(vabc,p_time,delta,vdq0);
  double Vd, Vq;
  Vd = vdq0[0]; Vq = vdq0[1];
  
  Ec = sqrt(Vd*Vd + Vq*Vq);

  double u1, u2, u3, u4;
  
  if(!zero_TR) {
    Vmeas = Filter_blkR.getoutput(Ec, timestep, true);
  } else {
    Vmeas = Ec;
  }
  u1 = Vref - Vmeas + Vuel + Vs;
  
  if(!zero_KIR) {
    u2 = PIControl_blkR.getoutput(u1, timestep, true);
  } else {
    u2 = Kpr * u1;
  }
  
  if(!zero_TA) {
    u3 = Filter_blkA.getoutput(u2, timestep, true);
  } else {
    u3 = u2;
  }
  
  if(!zero_KIM) {
    u4 = PIControl_blmM.getoutput(u3 - Efd * Kg, timestep, true);
  } else {
    u4 = Kpm * (u3 - Efd * Kg);
  }
  u4 = LVGate_blk.getoutput(u4);
  double Vb = CalculateVb(Vterm, Theta, Ir, Ii, LadIfd); 
  Efd = u4 * Vb;

}

void Esst4bExc::postStep(double time)
{

}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to
 * Esst4bExc
 */
void Esst4bExc::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  if (!data->getValue(EXCITER_TR, &Tr, idx)) Tr = 0.0; // Tr
  //if (!data->getValue(EXCITER_KPR, &Kpr, idx)) 
  Kpr = 0.0; // TBD: Kpr
  //if (!data->getValue(EXCITER_KIR, &Kir, idx)) 
  Kir = 0.0; // TBD: Kir
  if (!data->getValue(EXCITER_VRMAX, &Vrmax, idx)) Vrmax = 0.0; // Vrmax
  if (!data->getValue(EXCITER_VRMIN, &Vrmin, idx)) Vrmin = 0.0; // Vrmin
  if (!data->getValue(EXCITER_TA, &Ta, idx)) Ta = 0.0; // Ta
  //if (!data->getValue(EXCITER_KPM, &Kpm, idx)) 
  Kpm = 0.0; // TBD: Kpm
  //if (!data->getValue(EXCITER_KIM, &Kim, idx)) 
  Kim = 0.0; // TBD: Kim
  //if (!data->getValue(EXCITER_VMMAX, &Vmmax, idx)) 
  Vmmax = 0.0; // TBD: Vmmax
  //if (!data->getValue(EXCITER_VMMIN, &Vmmin, idx)) 
  Vmmin = 0.0; // TBD: Vmmin
  //if (!data->getValue(EXCITER_KG, &Kg, idx)) 
  Kg = 0.0; // TBD: Kg
  //if (!data->getValue(EXCITER_KP, &Kp, idx)) 
  Kp = 0.0; // TBD: Kp
  //if (!data->getValue(EXCITER_KI, &KI, idx)) 
  KI = 0.0; // TBD: KI
  //if (!data->getValue(EXCITER_VBMAX, &Vbmax, idx)) 
  Vbmax = 0.0; // TBD: Vbmax
  //if (!data->getValue(EXCITER_KC, &Kc, idx)) 
  Kc = 0.0; // TBD: Kc
  //if (!data->getValue(EXCITER_XL, &Xl, idx)) 
  Xl = 0.0; // TBD: Xl
  //if (!data->getValue(EXCITER_KPANG, &Kpang, idx)) 
  Kpang = 0.0; // TBD: Kpang
  //if (!data->getValue(EXCITER_VGMAX, &Vgmax, idx)) 
  Vgmax = 0.0; // TBD: Vgmax

  //printf("Tr = %f, Vrmax = %f, Vrmin = %f, Ta = %f\n", Tr, Vrmax, Vrmin, Ta); exit(1); 

  if (integrationtype != IMPLICIT) {
    // block-based explicit implementation
    // right now we just hard code Vuel, Voel and Vs
    Vs = 0.0;
    Vuel = 0.0;
    Voel = 1000.0;
    OptionToModifyLimitsForInitialStateLimitViolation = true;
    zero_TR = false;
    zero_TA = false;
    zero_KIM = false;
    zero_KIR = false;
  }
}

/**
 * Saturation function
 * @ param x
 */
double Esst4bExc::Sat(double x)
{
    /*double a_ = S12 / S10 - 1.0 / 1.2;
    double b_ = (-2 * S12 / S10 + 2);
    double c_ = S12 / S10 - 1.2;
    double A = (-b_ - sqrt(b_ * b_ - 4 * a_ * c_)) / (2 * a_);
    double B = S10 / ((1.0 - A) * (1.0 - A));
    return B * ( x - A) * (x - A);*/
    //double B = log(SE2 / SE1)/(E2 - E1);
    //double A = SE1 / exp(B * E1);
    //return A * exp(B * x);
}

double Esst4bExc::sqr(double x)
{
  return x * x;
}

/**
 * FEX function
 * @ param IN
 */
double Esst4bExc::FEX(double IN)
{
  double result;
  if (IN <= 0.0) result = 1;
  else if (IN <= 0.433) result = 1 - 0.577 * IN;
  else if (IN <= 0.75) result = sqrt(0.75 - sqr(IN));
  else if (IN < 1.0) result = 1.732 * (1 - IN);
  else result = 0;
  return result;
}

/**
 * CalculateVb function
 * @ param Vterm, Theta, Ir, Ii, LadIfd
 */
double Esst4bExc::CalculateVb(double Vterm,
  double Theta, double Ir, double Ii, double LadIfd)
{
  double pi = 4.0 * atan(1.0);
  // Calculate complex values for VE Calculation
  Kpvr = Kp * cos(Kpang * 180 / pi); 
  Kpvi = Kp * sin(Kpang * 180 / pi); 
  Kpir = - Kpvi * Xl;
  Kpii = + Kpvr * Xl + KI;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  double Ve = sqr(Kpvr * Vrterm - Kpvi * Viterm + Kpir * Ir - Kpii * Ii) +
              sqr(Kpvr * Viterm + Kpvi * Vrterm + Kpir * Ii + Kpii * Ir);
  Ve = sqrt(Ve);
  double Vb = Ve * FEX(Kc * LadIfd / Ve); // FEX function from 1.8.1
  if (Vb > Vbmax) Vb = Vbmax;
  return Vb; 
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void Esst4bExc::init(gridpack::RealType* xin)
{
  gridpack::RealType *x = xin+offsetb; // exciter array starts from this location
  double Ec = sqrt(VD*VD + VQ*VQ);

  // Get initial field voltage
  Efd = getInitialFieldVoltage();

  if(integrationtype != IMPLICIT) {
    // Set up blocks
    // Set parameters for the first block
    /*if (Tr < TS_THRESHOLD * ts) zero_TR = true; // how to get ts?
    if (Ta < TS_THRESHOLD * ts) zero_TA = true;
    if (Kim < KI_THRESHOLD) zero_KIM = true;
    if (Kir < KI_THRESHOLD) zero_KIR = true;
    
    if (zero_TR) printf("Tr=%f is better at least %d times larger than timestep=%f.\n", Tr, TS_THRESHOLD, ts);
    if (zero_TA) printf("Ta=%f is better at least %d times larger than timestep=%f.\n", Ta, TS_THRESHOLD, ts);
    if (zero_KIM) printf("Kim=%f is less than %d times of timestep=%f, treated as zero.\n", Kim, TS_THRESHOLD, ts);
    if (zero_KIR) printf("Kir=%f is less than %d times of timestep=%f, treated as zero.\n", Kir, TS_THRESHOLD, ts);*/
    
    // printf("print: inside esst4b model, Ir=%f, Ii=%f\n", Ir, Ii);
    
    if (Kpr == 0.0 && Kir < KI_THRESHOLD) {
      Kpr = 1.0;
      printf("Kpr and Kir cannot be both zeros; reset Kpr=1.0.\n");
    }
    
    if (Kpm == 0.0 && Kim < KI_THRESHOLD) {
      Kpm = 1.0;
      printf("Kpm and Kim cannot be both zeros; reset Kpm=1.0.\n");
    }

    if (!zero_TR) {
      Filter_blkR.setparams(1.0, Tr);
    }
    /*---else {a gain=1 block}---*/

    if (!zero_KIR) {
      PIControl_blkR.setparams(Kpr, Kir, Vrmin, Vrmax, -10000.0, 10000.0);
    }
    /*---else {a gain=Kpr block}---*/
    
    if (!zero_TA) {
      Filter_blkA.setparams(1.0, Ta);
    }
    
    if (!zero_KIM) {
      PIControl_blmM.setparams(Kpm, Kim, Vmmin, Vmmax, -10000.0, 10000.0);
    }
    
    LVGate_blk.setparams(Voel);
    
    double u1, u2, u3, u4;

    //Vterm = mag; // Ec
    //Theta = ang; // what's ang?
    Vterm = Ec;

    double Vb = CalculateVb(Vterm, Theta, Ir, Ii, LadIfd); 

    u1 = Efd / Vb;

    // LV Gate?
    // assume LV gate always take feedforward path during initialization, may need to adjust Voel
    if (OptionToModifyLimitsForInitialStateLimitViolation) {
      if (u1 > Voel) Voel = u1 + 0.1;
      LVGate_blk.setparams(Voel);
    }
    
    if (!zero_KIM) {
      u2 = PIControl_blmM.init_given_y(u1);
      // Check limits here, but these would be 
      // initial state limit violations that are not possible!
      if (OptionToModifyLimitsForInitialStateLimitViolation) {
        if (u1 > Vmmax) Vmmax = u1+0.1;
        if (u1 < Vmmin) Vmmin = u1-0.1;
      }
    } else {
      u2 = u1/Kpm;
    }
    
    if (!zero_TA) {
      u3 = Filter_blkA.init_given_y(u2 + Efd * Kg);
    } else {
      u3 = u2 + Efd * Kg;
    }
    
    if (!zero_KIR) {
      u4 = PIControl_blkR.init_given_y(u3);
      // Check limits here, but these would be 
      // initial state limit violations that are not possible!
      if (OptionToModifyLimitsForInitialStateLimitViolation) {
        if (u3 > Vrmax) Vrmax = u3+0.1;
        if (u3 < Vrmin) Vrmin = u3-0.1;
      }
    } else {
      u4 = u3/Kpr;
    }

    if(!zero_TR) {
      Vmeas = Filter_blkR.init_given_u(Vcomp);
    } else 
      Vmeas = Vcomp;
    Vref = u4 + Vmeas - Vuel - Vs; 

  } else { // Implicit implementation
    /*Vterm = mag;
    presentMag = mag;
    Theta = ang;
    presentAng = ang;
    // State 1
    double Vb = CalculateVb(Vterm, Theta, Ir, Ii, LadIfd); // TBD: What's the init value of Ir and Ii?
    printf("esst4b: Efd = %f\n", Efd);
    double Vm = Efd/ Vb;
    // Check limits here, but these would be 
    // initial state limit violations that are not possible!
    if (OptionToModifyLimitsForInitialStateLimitViolation) { // TBD: inital value of this bool? 
      if (Vm > Vmmax) Vmmax = Vm;
      if (Vm < Vmmin) Vmmin = Vm;
    }
    double TempIn;
    if (Kim != 0) {
      x1Vm = Vm;
      TempIn = 0;
    } else {
      x1Vm = 0;
      TempIn = Vm / Kpm;
    }
    // State 3
    x3Va = Efd * Kg;
    if (x3Va > Vgmax) x3Va = Vgmax;
    x3Va = x3Va + TempIn;
    // State 4
    double Vr = x3Va;
    // Check limits here, but these would be
    // initial state limit violations that are not possible!
    if (OptionToModifyLimitsForInitialStateLimitViolation) {
      if (Vr > Vrmax) Vrmax = Vr;
      if (Vr < Vrmin) Vrmin = Vr;
    }
    if (Kir != 0) {
      x4Vr = Vr;
      TempIn = 0;
    } else {
      x4Vr = 0;
      TempIn = Vr / Kpr;
    }
    // State 2
    x2Vcomp = Vcomp;  // TBD: init value of Vcomp?
    // Vref
    Vref = Vcomp + TempIn;

    //printf("esst4b init:  %f\t%f\t%f\t%f\n", x1Vm, x2Vcomp, x3Va, x4Vr); 
    x[0] = x1Vm;
    x[1] = x2Vcomp;
    x[2] = x3Va;
    x[3] = x4Vr;*/
  }
}

/**
 * Write output from exciters to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Esst4bExc::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out exciter state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Esst4bExc::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto exciters
 * @param values array containing exciter state variables
*/
void Esst4bExc::setValues(gridpack::RealType *val)
{
  gridpack::RealType *values = val+offsetb; // exciter array starts from this location

  if(integrationtype == EXPLICIT) return;

  if(p_mode == XVECTOBUS) {
    x1Vm = values[0];
    x2Vcomp = values[1];
    x3Va = values[2];
    x4Vr = values[3];
  } else if(p_mode == XDOTVECTOBUS) {
    dx1Vm = values[0];
    dx2Vcomp = values[1];
    dx3Va = values[2];
    dx4Vr = values[3];
  }
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Esst4bExc::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // exciter array starts from this location

  if(integrationtype == EXPLICIT) return;

  // TBD
}

/**
   Non-zero pattern of the Jacobian (x denotes non-zero entry)
         Vmeas   xLL     Efd     delta    va    vb    vc
 eq.0 |    x                       x       x     x     x  
 eq.1 |    x      x              
 eq.2 |    x      x       x      

 Number of non-zeros = 5 + 2 + 3 = 10 
 * Get number of matrix values contributed by exciter
 * @return number of matrix values
 */
int Esst4bExc::matrixNumValues()
{
  int nmat = 0;
  if(integrationtype == IMPLICIT) nmat = 15; /* extra 5 spaces for padding */ // To be checked!
  return nmat;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Esst4bExc::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;

  if(integrationtype != IMPLICIT) {
    *nvals = ctr;
    return;
  }

  //TBD: implicit implementation
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
/*void gridpack::dynamic_simulation::Esst4bExc::predictor(double t_inc, bool flag)
{
  if (!flag) {
    x1Vm = x1Vm_1;
    x2Vcomp = x2Vcomp_1;
    x3Va = x3Va_1;
    x4Vr = x4Vr_1;
  }
  // State 2
  if (Tr < TS_THRESHOLD * t_inc) {
    x2Vcomp = Vcomp; // Must propogate input value instantaneously
    dx2Vcomp = 0;
  } else {
    dx2Vcomp = 1 / Tr * (Vcomp - x2Vcomp);
  }
  // State 4
  double TempIn = -x2Vcomp + Vref + Vstab;// + Vuel; // TBD: what is Vuel?
  double TempMax = Vrmax - TempIn * Kpr;
  double TempMin = Vrmin - TempIn * Kpr;
  if (x4Vr > TempMax) x4Vr = TempMax;
  else if (x4Vr < TempMin) x4Vr = TempMin;
  dx4Vr = Kir * TempIn;
  if (dx4Vr > 0 && x4Vr >= TempMax) dx4Vr = 0;
  else if (dx4Vr < 0 && x4Vr <= TempMin) dx4Vr = 0;
  //State 3
  TempIn = x4Vr + TempIn * Kpr;
  if (Ta < TS_THRESHOLD * t_inc) {
    x3Va = TempIn; // Must propogate input value instantaneously
    dx3Va = 0;
  } else {
    dx3Va = 1 / Tr * (TempIn - x3Va);
  }
  //State 1
  TempIn = Efd * Kg;
  if (TempIn > Vgmax) TempIn = Vgmax;
  TempIn = x3Va - TempIn;
  TempMax = Vmmax - TempIn * Kpm;
  TempMin = Vmmin - TempIn * Kpm;
  if (x1Vm > TempMax) x1Vm = TempMax;
  else if (x1Vm < TempMin) x1Vm = TempMin;
  dx1Vm = Kim * TempIn;
  if (dx1Vm > 0 && x1Vm >= TempMax) dx1Vm = 0;
  else if (dx1Vm < 0 && x1Vm <= TempMin) dx1Vm = 0;

  x1Vm_1 = x1Vm + dx1Vm * t_inc;
  x2Vcomp_1 = x2Vcomp + dx2Vcomp * t_inc;
  x3Va_1 = x3Va + dx3Va * t_inc;
  x4Vr_1 = x4Vr + dx4Vr * t_inc;

  printf("esst4b dx: %f\t%f\t%f\t%f\t\n", dx1Vm, dx2Vcomp, dx3Va, dx4Vr);
  printf("esst4b x: %f\t%f\t%f\t%f\n", x1Vm_1, x2Vcomp_1, x3Va_1, x4Vr_1);

  double Vb = CalculateVb(Vterm, Theta, Ir, Ii, LadIfd);
  //if (x1Vm > Voel) TempIn = Voel * Vb; // TBD: what is Voel?
  //else Efd = x1Vm_1 * Vb;
  Efd = x1Vm_1 * Vb; // TBD: temporailly

  printf("esst4b Efd: %f\n", Efd);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
/*void gridpack::dynamic_simulation::Esst4bExc::corrector(double t_inc, bool flag)
{
  // State 2
  if (Tr < TS_THRESHOLD * t_inc) {
    x2Vcomp_1 = Vcomp; // Must propogate input value instantaneously
    dx2Vcomp_1 = 0;
  } else {
    dx2Vcomp_1 = 1 / Tr * (Vcomp - x2Vcomp_1);
  }
  // State 4
  double TempIn = -x2Vcomp_1 + Vref + Vstab;// + Vuel; // TBD: what is Vuel?
  double TempMax = Vrmax - TempIn * Kpr;
  double TempMin = Vrmin - TempIn * Kpr;
  if (x4Vr_1 > TempMax) x4Vr_1 = TempMax;
  else if (x4Vr_1 < TempMin) x4Vr_1 = TempMin;
  dx4Vr_1 = Kir * TempIn;
  if (dx4Vr_1 > 0 && x4Vr_1 >= TempMax) dx4Vr_1 = 0;
  else if (dx4Vr_1 < 0 && x4Vr_1 <= TempMin) dx4Vr_1 = 0;
  //State 3
  TempIn = x4Vr_1 + TempIn * Kpr;
  if (Ta < TS_THRESHOLD * t_inc) {
    x3Va_1 = TempIn; // Must propogate input value instantaneously
    dx3Va_1 = 0;
  } else {
    dx3Va_1 = 1 / Tr * (TempIn - x3Va_1);
  }
  //State 1
  TempIn = Efd * Kg;
  if (TempIn > Vgmax) TempIn = Vgmax;
  TempIn = x3Va_1 - TempIn;
  TempMax = Vmmax - TempIn * Kpm;
  TempMin = Vmmin - TempIn * Kpm;
  if (x1Vm_1 > TempMax) x1Vm_1 = TempMax;
  else if (x1Vm_1 < TempMin) x1Vm_1 = TempMin;
  dx1Vm_1 = Kim * TempIn;
  if (dx1Vm_1 > 0 && x1Vm_1 >= TempMax) dx1Vm_1 = 0;
  else if (dx1Vm_1 < 0 && x1Vm_1 <= TempMin) dx1Vm_1 = 0;

  x1Vm_1 = x1Vm + (dx1Vm + dx1Vm_1) / 2.0 * t_inc;
  x2Vcomp_1 = x2Vcomp + (dx2Vcomp + dx2Vcomp_1) / 2.0 * t_inc;
  x3Va_1 = x3Va + (dx3Va + dx3Va_1) / 2.0 * t_inc;
  x4Vr_1 = x4Vr + (dx4Vr + dx4Vr_1) / 2.0 * t_inc;

  printf("esst4b dx: %f\t%f\t%f\t%f\t\n", dx1Vm_1, dx2Vcomp_1, dx3Va_1, dx4Vr_1);
  printf("esst4b x: %f\t%f\t%f\t%f\n", x1Vm_1, x2Vcomp_1, x3Va_1, x4Vr_1);
  
  double Vb = CalculateVb(Vterm, Theta, Ir, Ii, LadIfd);
  //if (x1Vm > Voel) TempIn = Voel * Vb; // TBD: what is Voel?
  //else Efd = x1Vm_1 * Vb;
  Efd = x1Vm_1 * Vb; // TBD: temporially

  printf("esst4b Efd: %f\n", Efd);
}*/

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
double Esst4bExc::getFieldVoltage()
{
  return Efd;
}

/** 
 * Get the value of the field voltage parameter
 * and its global location
 * @return value of field voltage
 */
double Esst4bExc::getFieldVoltage(int *Efd_gloc)
{
  if(integrationtype == IMPLICIT) {
    *Efd_gloc = p_gloc + 2;
  } else {
    *Efd_gloc = -1;
  }
  return Efd;
}

bool Esst4bExc::getFieldVoltagePartialDerivatives(int *xexc_loc,double *dEfd_dxexc,double *dEfd_dxgen)
{
  return false;
}

/**
 * Update the event function values
 */
void Esst4bExc::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType>& evalues)
{

}

/**
 * Event handler
 */
void Esst4bExc::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{

}

/**
 * Set event
 */
void Esst4bExc::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{

}

/**
 * Set the field voltage parameter inside the exciter
 * @param fldv value of the field voltage
 */
void Esst4bExc::setFieldVoltage(double fldv)
{
  Efd = fldv;
}

/**
 * Set the field current parameter inside the exciter
 * @param fldc value of the field current
 */
void Esst4bExc::setFieldCurrent(double fldc)
{
  LadIfd = fldc;
}

/** 
 * Get the value of the field current parameter
 * @return value of field current
 */
double Esst4bExc::getFieldCurrent()
{
  return 0.0;
}

/** 
 * Set the value of the Vterminal
 * @return value of field current
 */
void Esst4bExc::setVterminal(double mag)
{
  //Vterminal = mag;
}

/** 
 * Set the value of the omega
 * @return value of field current
 */
void Esst4bExc::setOmega(double omega)
{
  //w = omega;
}

void Esst4bExc::setVuel(double vtmp)
{
  Vuel = vtmp;
}

void Esst4bExc::setVs(double vtmp)
{
  Vs = vtmp;
}

void Esst4bExc::setIri(double vIr, double vIi)
{
  Ir = vIr;
  Ii = vIi;
}

void Esst4bExc::setVoel(double vtmp)
{
  Voel = vtmp;
}

/** 
 * Set the value of the Vcomp
 * @return value of the Vcomp
 */
void Esst4bExc::setVcomp(double vtmp)
{
  Vcomp = vtmp;
}
