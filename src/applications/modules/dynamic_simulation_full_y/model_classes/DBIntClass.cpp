/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   DBIntClass.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 11, 2015
 * 
 * @brief  
 * 
 * 
 */

#include <cmath>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "DBIntClass.hpp"

/**
 * Basic constructor
 */
gridpack::dynamic_simulation::DBIntClass::DBIntClass(void)
{
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::DBIntClass::~DBIntClass(void)
{
}

/**
 * @param theDb1 
 * @param theEps
 * @param theInit
 */
void gridpack::dynamic_simulation::DBIntClass::Initialize(double theDb1, double theEps, double theInit)
{
  InitialValue = theInit;
  State = 0;
  Db1 = abs(theDb1); // Negatives not allowed
  Eps = abs(theEps);
  if (Eps > Db1) Eps = Db1; // Don't allow it to be larger
  //printf("theDb1 = %f, theEps = %f, theInit = %f\n", theDb1, theEps, theInit);
  //printf("Db1 = %f, Eps = %f, InitialValue = %f\n", Db1, Eps, InitialValue);
}

/**
 * @param anInput
 */
double gridpack::dynamic_simulation::DBIntClass::Output(double anInput)
{
  double Result;
  anInput = anInput - InitialValue; // InitialValue almost always zero
  // Most common situation when deadband resulting in no output
  if (abs(anInput) <= Db1 && ((Eps == 0) || (State == 0))) 
    Result = 0;
  // No Hysteresis
  else if (Eps == 0) { // Must be eitehr > Db1 or < -Db1
    if (anInput > Db1) Result = anInput - Db1;
    else Result = anInput + Db1;
  // Check for bowtie hysteresis
  } else if (Db1 == Eps) {
    if (anInput >= Db1) {
      Result = anInput;
      State = 1; // On positive side
    } else if (anInput <= -Db1) {
      Result = anInput;
      State = 2;
    } else if (State == 1) {
      if (anInput > 0)
        Result = anInput; // Between 0 and Db1
      else { // <= 0
        Result = 0;
        State = 0;
      }
    } else if (State == 2) {
      if (anInput < 0)
        Result = anInput; // Between -Db1 and 0
      else { // >= 0
        Result = 0;
        State = 0; // Reset
      }
    }
  // Intentional Deadband with Partial Hysteresis
  // Do remaining case with Eps > 0 and Eps < Db1.
  // Note we covered the case with input < Db1 with State=0 above
  } else {
    if (Result >= Db1) {
      Result = anInput - Db1 + Eps;
      State = 1;
    } else if (anInput <= -Db1) {
      Result = anInput + Db1 - Eps;
      State = 2;
    } else { // Cover case of close to zero
      if (abs(anInput) <= Db1 - Eps) {
        Result = 0;
        State = 0;
      } else if (anInput > 0) { // Between Db1-Eps and Db1
        if (State == 1)
          Result = anInput - Db1 + Eps;
        else { // Going from negative to positive
          State = 0;
          Result = 0;
        }
      } else if (anInput < 0) { // Between -Db1 and -(Db1-Eps)
        if (State == 2)
          Result = anInput + Db1 - Eps;
        else { // Going from positive to negative
          State = 0;
          Result = 0;
        }
      }
    }
  }
  Result = Result + InitialValue; 
  return Result;
}
