/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   utilfunctions.cpp
 * 
 * @brief  Utility functions source
 * 
 * 
 */
// -------------------------------------------------------------

#include <emtutilfunctions.hpp>

/**
 * abc2dq0 - Does a vector transform from abc to dq0 reference frame using the given angle theta
 * @param[input] xabc - the vector in abc reference frame represented as an array with 3 elements
 * @param[input] t - current time
 * @param[input] theta - the angle of rotation
 * @param[output] xdq0 - the transformed vector in dq0 reference frame
 */
void abc2dq0(double *xabc,double t,double theta,double *xdq0)
{
  double omegat = OMEGA_S*t;
  double omegat_minus = omegat - 2.0*PI/3.0;
  double omegat_plus = omegat + 2.0*PI/3.0;

  xdq0[0] = TWO_OVER_THREE*(sin(omegat)*xabc[0] + sin(omegat_minus)*xabc[1] + sin(omegat_plus)*xabc[2]);
  xdq0[1] = TWO_OVER_THREE*(cos(omegat)*xabc[0] + cos(omegat_minus)*xabc[1] + cos(omegat_plus)*xabc[2]);
  xdq0[2] = TWO_OVER_THREE*0.5*(xabc[0] + xabc[1] + xabc[2]);
}

/**
 * dq02abc - Does a vector transform from dq0 to abc reference frame using the given angle theta
 * @param[input] xdq0 - the vector in dq0 reference frame represented as an array with 3 elements
 * @param[input] t - current time
 * @param[input] theta - the angle of rotation
 * @param[output] xabc - the transformed vector in abc reference frame
 */
void dq02abc(double *xdq0,double t,double theta,double *xabc)
{
  double omegat = OMEGA_S*t;
  double omegat_minus = omegat - 2.0*PI/3.0;
  double omegat_plus = omegat + 2.0*PI/3.0;

  xabc[0] = sin(omegat)*xabc[0]       + cos(omegat)*xdq0[1]       + xdq0[2];
  xabc[1] = sin(omegat_minus)*xabc[0] + cos(omegat_minus)*xdq0[1] + xdq0[2];
  xabc[0] = sin(omegat_plus)*xabc[0]  + cos(omegat_plus)*xdq0[1]  + xdq0[2];  
}
