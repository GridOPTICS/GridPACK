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


// 2x2 determinant
double determinant2x2(double a, double b, double c, double d)
{
  return a * d - b * c;
}

// 3x3 determinant
double determinant3x3(double matrix[3][3])
{
  double det = 0.0;
  
  // Calculate the determinant using the formula for 3x3 matrices
  det = matrix[0][0] * determinant2x2(matrix[1][1], matrix[1][2], matrix[2][1], matrix[2][2])
    - matrix[0][1] * determinant2x2(matrix[1][0], matrix[1][2], matrix[2][0], matrix[2][2])
    + matrix[0][2] * determinant2x2(matrix[1][0], matrix[1][1], matrix[2][0], matrix[2][1]);
  
  return det;
}

// 3x3 inverse
void inverse3x3(double matrix[3][3], double result[3][3])
{
  double det = determinant3x3(matrix);
  
  if (det == 0.0) {
    printf("Matrix is not invertible (determinant is zero).\n");
    return;
  }
  
  // Calculate the inverse using the formula for 3x3 matrices
  result[0][0] = determinant2x2(matrix[1][1], matrix[1][2], matrix[2][1], matrix[2][2]) / det;
  result[0][1] = -determinant2x2(matrix[0][1], matrix[0][2], matrix[2][1], matrix[2][2]) / det;
  result[0][2] = determinant2x2(matrix[0][1], matrix[0][2], matrix[1][1], matrix[1][2]) / det;
  result[1][0] = -determinant2x2(matrix[1][0], matrix[1][2], matrix[2][0], matrix[2][2]) / det;
  result[1][1] = determinant2x2(matrix[0][0], matrix[0][2], matrix[2][0], matrix[2][2]) / det;
  result[1][2] = -determinant2x2(matrix[0][0], matrix[0][2], matrix[1][0], matrix[1][2]) / det;
  result[2][0] = determinant2x2(matrix[1][0], matrix[1][1], matrix[2][0], matrix[2][1]) / det;
  result[2][1] = -determinant2x2(matrix[0][0], matrix[0][1], matrix[2][0], matrix[2][1]) / det;
  result[2][2] = determinant2x2(matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1]) / det;
}

// 3x3 scaled matrix multiply - multiplies two 3x3 matrices A, B; scales by given value alpha, and returns the resultant matrix C
// C = alpha*A*B
void scaledmatmatmult3x3(double A[3][3],double B[3][3],double C[3][3],double alpha)
{
  for(int i=0; i < 3; i++) {
    for(int j=0; j < 3; j++) {
      C[i][j] = 0.0;
      for(int k=0; k < 3; k++) {
	C[i][j] += A[i][k]*B[k][j];
      }
      C[i][j] *= alpha;
    }
  }
}

void matmatmult3x3(double A[3][3],double B[3][3],double C[3][3])
{
  scaledmatmatmult3x3(A,B,C,1.0);
}

// y = alpha*A*x
void scaledmatvecmult3x3(double A[3][3],double *x,double *y,double alpha)
{
  y[0] = y[1] = y[2] = 0.0;

  y[0] = alpha*A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
  y[1] = alpha*A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
  y[2] = alpha*A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];
}

// y = Ax
void matvecmult3x3(double A[3][3],double *x,double *y)
{
  scaledmatvecmult3x3(A,x,y,1.0);
}


