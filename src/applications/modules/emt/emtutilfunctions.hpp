/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   emtutilfunctions.hpp
 * 
 * @brief  Utility functions header
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _emtutilfunctions_h
#define _emtutilfunctions_h

#include <gridpack/applications/modules/emt/constants.hpp>

/**
 * abc2dq0 - Does a vector transform from abc to dq0 reference frame using the given angle theta
 * @param[input] xabc - the vector in abc reference frame represented as an array with 3 elements
 * @param[input] t - current time
 * @param[input] theta - the angle of rotation
 * @param[output] xdq0 - the transformed vector in dq0 reference frame
 */
void abc2dq0(double *xabc,double t,double theta,double *xdq0);

/**
 * dq02abc - Does a vector transform from dq0 to abc reference frame using the given angle theta
 * @param[input] xdq0 - the vector in dq0 reference frame represented as an array with 3 elements
 * @param[input] t - current time
 * @param[input] theta - the angle of rotation
 * @param[output] xabc - the transformed vector in abc reference frame
 */
void dq02abc(double *xdq0,double t,double theta,double *xabc);

/**
  getTdq0 - returns the transformation matrix used for abc to dq0 transform
  @param[input] t     - current time
  @param[input] theta - the transform angle
  @param[output] Tdq0 - abc2dq0 transformation matrix
*/
void getTdq0(double t, double theta, double Tdq0[][3]);

/**
  getdTdq0dtheta - returns the partial derivative of transformation matrix Tdq0 w.r.t angle theta
  @param[input] t             - current time
  @param[input] theta         - the transform angle
  @param[output] dTdq0_ddelta - \partial(Tdq0){delta}
*/
void getdTdq0dtheta(double t, double theta, double dTdq0dtheta[][3]);

/**
  getTdq0inv - returns the transformation matrix used for dq0 to abc transform
  @param[input] t     - current time
  @param[input] theta - the transform angle
  @param[output] Tdq0inv - dq02abc transformation matrix
*/
void getTdq0inv(double t, double theta, double Tdq0inv[][3]);

/**
  getdTdq0invdtheta - returns the partial derivative of transformation matrix Tdq0inv w.r.t angle theta
  @param[input] t             - current time
  @param[input] theta         - the transform angle
  @param[output] dTdq0inv_ddelta - \partial(Tdq0inv){delta}
*/
void getdTdq0invdtheta(double t, double theta, double dTdq0invddelta[][3]);


//inverse = matrix^-1
void inverse3x3(double matrix[3][3],double inverse[3][3]);

// C = alpha*A*B
void scaledmatmatmult3x3(double A[3][3],double B[3][3],double C[3][3],double alpha);

// C = A*B
void matmatmult3x3(double A[3][3],double B[3][3], double C[3][3]);

// y = alpha*A*x
void scaledmatvecmult3x3(double A[3][3], double *x,double *y,double alpha);

// y = A*x
void matvecmult3x3(double A[3][3],double *x, double *y);

// y = alpha*A^T*x
void scaledmattransposevecmult3x3(double A[3][3], double *x,double *y,double alpha);

// y = A^T*x
void mattransposevecmult3x3(double A[3][3],double *x,double *y);

// Dot product: result = u^T*v
void vecdot3(double u[3], double v[3], double *result);

// Multiplies a 1 x 3 vector with a 3 x 3. Result is a 1 X 3 vector
void vec3multmat3x3(double vec[3], double A[3][3],double *result);

// Multiplies a 1 x 3 vector with a 3 x 3. Result is a 1 X 3 vector
// y = alpha*A*x
void scaledvec3multmat3x3(double vec[3], double A[3][3],double *result,double alpha);


#endif
