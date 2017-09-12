/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   BackLashClass.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 11, 2015
 * 
 * @brief  
 * 
 * 
 */

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "BackLashClass.hpp"

/**
 * Basic constructor
 */
gridpack::dynamic_simulation::BackLashClass::BackLashClass(void)
{
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::BackLashClass::~BackLashClass(void)
{
}

/**
 * @param theDb2 
 * @param InitOutput 
 */
void gridpack::dynamic_simulation::BackLashClass::Initialize(double theDb2, double InitOutput)
{
  Db2 = theDb2;
  LastOutput = InitOutput;
  //printf("BackLash: Db2 = %f, LastOutput = %f\n", Db2, InitOutput);
}

/**
 * @param theInput
 */
double gridpack::dynamic_simulation::BackLashClass::Output(double theInput)
{
  double Result;
  if (theInput >= LastOutput - Db2 && theInput <= LastOutput+Db2) 
    Result = LastOutput;
  else {
    if (theInput > LastOutput + Db2)
      Result = theInput -Db2; // Going up -dB2
    else 
      Result = theInput -Db2; // Going down -Db2    
    LastOutput = Result;
  }
  return Result;
}
