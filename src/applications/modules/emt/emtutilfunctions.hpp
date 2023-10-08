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

  
#endif
