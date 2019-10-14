/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   BackLashClass.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   10/02/2019
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
BackLashClass::BackLashClass(void)
{
}

/**
 * Basic destructor
 */
BackLashClass::~BackLashClass(void)
{
}

/**
 * @param theDb2 
 * @param InitOutput 
 */
void BackLashClass::Initialize(double theDb2, double InitOutput)
{
  Db2 = theDb2;
  LastOutput = InitOutput;
  //printf("BackLash: Db2 = %f, LastOutput = %f\n", Db2, InitOutput);
}

/**
 * @param theInput
 */
double BackLashClass::Output(double theInput)
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
