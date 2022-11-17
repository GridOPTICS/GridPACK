#include "dblock.hpp"

// ------------------------------------
// Gain
// ------------------------------------

Gain::Gain(void)
{
}

void Gain::setparams(double K,double ymin,double ymax)
{
  p_K = K;
  p_ymin = ymin;
  p_ymax = ymax;
}

double Gain::getoutput(double u)
{
  double yout;

  yout = std::max(p_ymin,std::min(p_K*u,p_ymax));

  return yout;
}

// ------------------------------------
// Slope
// ------------------------------------

Slope::Slope(void)
{
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

