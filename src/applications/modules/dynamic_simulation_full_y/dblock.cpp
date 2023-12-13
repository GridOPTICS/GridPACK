#include "dblock.hpp"

// ------------------------------------
// GainLimiter
// ------------------------------------

GainLimiter::GainLimiter(void)
{
}

void GainLimiter::setparams(double K,double ymin,double ymax)
{
  p_K = K;
  p_ymin = ymin;
  p_ymax = ymax;
}

void GainLimiter::setparams(double K,double ymin,double ymax,double dymin, double dymax)
{
  p_K = K;
  p_ymin = ymin;
  p_ymax = ymax;
  p_dymin = fabs(dymin);
  p_dymax = dymax;

  p_yprev = 100000; // Some large value to initialize
}

double GainLimiter::getoutput(double u)
{
  double yout = 0.0;

  yout = std::max(p_ymin,std::min(p_K*u,p_ymax));

  return yout;
}

double GainLimiter::getoutput(double u, double ymin, double ymax)
{
  double yout = 0.0;

  p_ymin = ymin;
  p_ymax = ymax;

  yout = getoutput(u);

  return yout;
}

double GainLimiter::getoutput(double u, double dt, bool isupdate)
{
  double yout;

  yout = std::max(p_ymin,std::min(p_K*u,p_ymax));

  if(p_yprev == 100000) {
    // First time this function is being called, so initialize yprev
    p_yprev = yout;
  }

  double ymax_ratelimited = p_yprev + p_dymax*dt;
  double ymin_ratelimited = p_yprev - p_dymin*dt;

  yout  = std::max(ymin_ratelimited,std::min(yout,ymax_ratelimited));

  if(isupdate) p_yprev = yout;

  return yout;
}



// ------------------------------------
// Slope
// ------------------------------------

Slope::Slope(void)
{
  // Initialize to zero slope
  p_u[0] = 0.0;
  p_u[1] = 1.0;
  p_y[0] = p_y[1] = 0.0;
  p_ya = p_yb = 0.0;
}

void Slope::setparams(double *u,double *y,double ya, double yb)
{
  p_u[0] = u[0];
  p_u[1] = u[1];
  p_y[0] = y[0];
  p_y[1] = y[1];
  p_ya   = ya;
  p_yb   = yb;

  p_dydu = (p_y[1] - p_y[0])/(p_u[1] - p_u[0]);
}

double Slope::getoutput(double u)
{
  double yout;

  if(u - p_u[0] < 1e-5) yout = p_ya;
  else if(p_u[1] - u < 1e-5) yout = p_yb;
  else yout = p_y[0] + p_dydu*(u - p_u[0]);

  return yout;
}

// ------------------------------------
// Deadband
// ------------------------------------

Deadband::Deadband(void)
{
  p_u1 = p_u2 = 0.0;
}

void Deadband::setparams(double u1,double u2)
{
  p_u1 = u1;
  p_u2 = u2;
}

double Deadband::getoutput(double u)
{
  double yout = 0.0;

  if(u - p_u1 > 1e-6 && p_u2 - u > 1e-6) yout = 0.0; // In the deadband
  else yout = u;
  
  return yout;
}

//-------------------
//   PieceWiseSlope
//-------------------

PiecewiseSlope::PiecewiseSlope(void)
{
  p_n = 0;
  p_increasing = true;
}

void PiecewiseSlope::setparams(int n, double *uin, double *yin)
{
  double u[2],y[2],ya,yb;
  int    i;

  // Check the slope direction - increasing or decreasing
  if(uin[1] > uin[0]) p_increasing = true;
  else p_increasing = false;
  
  for(i=0; i < n-1; i++) {
    // Check if the segment is valid, i.e., monotonically increasing or decreasing
    if((p_increasing && uin[i+1] > uin[i]) || (!p_increasing && uin[i+1] < uin[i])) {
      p_u[p_n]   = uin[i];
      p_u[p_n+1] = uin[i+1];
      p_y[p_n]   = yin[i];
      p_y[p_n+1] = yin[i+1];

      p_slopes[p_n].setparams(p_u+p_n,p_y+p_n,p_y[p_n],p_y[p_n+1]);
      p_n++;
    }
  }
}

double PiecewiseSlope::getoutput(double u)
{
  int i;
  double yout;

  if(p_increasing  && u <= p_u[0])   yout = p_y[0];
  if(p_increasing  && u >= p_u[p_n]) yout = p_y[p_n];
  if(!p_increasing && u >= p_u[0])   yout = p_y[0];
  if(!p_increasing && u <= p_u[p_n]) yout = p_y[p_n];

  for(i=0; i < p_n; i++) {
    if(p_u[i] <= u && u <= p_u[i+1]) {
      yout = p_slopes[i].getoutput(u);
    }
  }

  return yout;
}
    
  
double PiecewiseSlope::init_given_y(double yout)
{
  int i;
  double u;
  
  if(p_increasing  && yout <= p_y[0])   u = p_u[0];
  if(p_increasing  && yout >= p_y[p_n]) u = p_u[p_n];
  if(!p_increasing && yout >= p_y[0])   u = p_u[0];
  if(!p_increasing && yout <= p_y[p_n]) u = p_u[p_n];

  for (i=0; i<p_n; i++) {
    if (p_y[i] <= yout && yout <= p_y[i+1]) {
      u = p_u[i] + (p_u[i+1] - p_u[i]) / (p_y[i+1] - p_y[i]) * (yout - p_y[i]);
    }
  }

  return u;
}

//---yuan add below---//
/**
 * Basic constructor
 */
Deadband_backlash::Deadband_backlash(void)
{
  p_Db2=0.0;
  p_type=1;  // default: type 1 deadband backlash
  p_LastOutput=0.0;
}

/**
 * @param theDb2 
 * @param InitOutput 
 */
void Deadband_backlash::setparams(double theDb2, int DbType)
{
  p_Db2 = theDb2;
  p_type = DbType;
}

/**
 * @param InitOutput 
 */
void Deadband_backlash::init_given_y(double InitOutput)
{
  p_LastOutput = InitOutput;
}

/**
 * @param theInput
 */
double Deadband_backlash::getoutput(double theInput)
{
  double Result;
  if (p_type ==  1) {
	  if (theInput >= p_LastOutput - p_Db2 && theInput <= p_LastOutput + p_Db2) 
		Result = p_LastOutput;
	  else {
		Result = theInput;
		p_LastOutput = Result;
	  }
  } else {  // type 2: identical to BashLashClass in old implementation
	  if (theInput >= p_LastOutput - p_Db2 && theInput <= p_LastOutput + p_Db2) 
		Result = p_LastOutput;
	  else {
		if (theInput > p_LastOutput + p_Db2) Result = theInput - p_Db2;
		else Result = theInput + p_Db2;  // however, in old implementation, this was "Result = theInput -Db2; // Going down -Db2"
		p_LastOutput = Result;
	  }
  }
  return Result;
}


/**
 * Basic constructor
 */
Deadband_dbint::Deadband_dbint(void)
{
    p_Db1 = 0.0;
    p_Eps = 0.0;
    p_State = 1;
    p_InitialValue = 0.0;
}

/**
 * @param theDb1 
 * @param theEps
 */
void Deadband_dbint::setparams(double theDb1, double theEps)
{
  p_State = 0;
  p_Db1 = abs(theDb1); // Negatives not allowed
  p_Eps = abs(theEps);
  if (p_Eps > p_Db1) p_Eps = p_Db1; // Don't allow it to be larger
  //printf("theDb1 = %f, theEps = %f, theInit = %f\n", theDb1, theEps, theInit);
  //printf("Db1 = %f, Eps = %f, InitialValue = %f\n", Db1, Eps, InitialValue);
}

/**
 * @param theInitInput
 */
void Deadband_dbint::init_given_u(double theInitInput)
{
  p_InitialValue = theInitInput;
}

/**
 * @param theInitOutput
 */
void Deadband_dbint::init_given_y(double theInitOutput)
{
  p_InitialValue = theInitOutput;
}

/**
 * @param anInput
 */
double Deadband_dbint::getoutput(double anInput)
{
  double Result;
  anInput = anInput - p_InitialValue; // InitialValue almost always zero
  // Most common situation when deadband resulting in no output
  if (abs(anInput) <= p_Db1 && ((p_Eps == 0) || (p_State == 0))) 
    Result = 0;
  // No Hysteresis
  else if (p_Eps == 0) { // Must be eitehr > Db1 or < -Db1
    if (anInput > p_Db1) Result = anInput - p_Db1;
    else Result = anInput + p_Db1;
  // Check for bowtie hysteresis
  } else if (p_Db1 == p_Eps) {
    if (anInput >= p_Db1) {
      Result = anInput;
      p_State = 1; // On positive side
    } else if (anInput <= -p_Db1) {
      Result = anInput;
      p_State = 2;
    } else if (p_State == 1) {
      if (anInput > 0)
        Result = anInput; // Between 0 and Db1
      else { // <= 0
        Result = 0;
        p_State = 0;
      }
    } else if (p_State == 2) {
      if (anInput < 0)
        Result = anInput; // Between -Db1 and 0
      else { // >= 0
        Result = 0;
        p_State = 0; // Reset
      }
    }
  // Intentional Deadband with Partial Hysteresis
  // Do remaining case with Eps > 0 and Eps < Db1.
  // Note we covered the case with input < Db1 with State=0 above
  } else {
    if (Result >= p_Db1) {
      Result = anInput - p_Db1 + p_Eps;
      p_State = 1;
    } else if (anInput <= -p_Db1) {
      Result = anInput + p_Db1 - p_Eps;
      p_State = 2;
    } else { // Cover case of close to zero
      if (abs(anInput) <= p_Db1 - p_Eps) {
        Result = 0;
        p_State = 0;
      } else if (anInput > 0) { // Between Db1-Eps and Db1
        if (p_State == 1)
          Result = anInput - p_Db1 + p_Eps;
        else { // Going from negative to positive
          p_State = 0;
          Result = 0;
        }
      } else if (anInput < 0) { // Between -Db1 and -(Db1-Eps)
        if (p_State == 2)
          Result = anInput + p_Db1 - p_Eps;
        else { // Going from positive to negative
          p_State = 0;
          Result = 0;
        }
      }
    }
  }
  Result = Result + p_InitialValue; 
  return Result;
}
//---yuan add above---//
