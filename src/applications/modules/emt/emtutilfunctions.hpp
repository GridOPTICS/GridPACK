/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   utilfunctions.hpp
 * 
 * @brief  Utility functions header
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _utilfunctions_h
#define _utilfunctions_h

#include <constants.hpp>

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

#endif
