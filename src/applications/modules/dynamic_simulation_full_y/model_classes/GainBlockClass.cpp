/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   GainBlockClass.cpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 24, 2015
 * 
 * @brief  
 * 
 * 
 */

#include "boost/smart_ptr/shared_ptr.hpp"
#include "GainBlockClass.hpp"

/**
 * Basic constructor
 */
gridpack::dynamic_simulation::GainBlockClass::GainBlockClass(void)
{
  Count = 0; // Default to 0 for now, need to be specified based on dyr
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::GainBlockClass::~GainBlockClass(void)
{
}

/**
 * @param theX 
 */
double gridpack::dynamic_simulation::GainBlockClass::XtoY(double theX)
{
  double Result;
  int I;
  if (Count == 0) 
    Result = theX;
  // TBD: need to work on the index
  /*else {
    if (theX < X[1]) Result = Y[1];
    else if (theX > X[Count]) {
      I = Count;
      if ((! ExtrapolatePastEnd) || (I < 2)) 
        Result = Y[I];
      else 
        Result = Y[I] + (Y[I] - Y[I-1]) / (X[I] - X[I-1]) * (theX - X[I]);
    } else {
      I = 2;
      while ((I < Count) && (X[I] < theX)) 
        I = I + 1;
      Result = Y[I-1] + (Y[I] - Y[I-1]) / (X[I] - X[I-1]) * (theX - X[I-1]);
    }
  }*/
  return Result;
}

/**
 * @param theY 
 */
double gridpack::dynamic_simulation::GainBlockClass::YtoX(double theY)
{
  double Result;
  int I;
  if (Count == 0) Result = theY; // return input if empty
  // TBD: need to work on the index
  /*else {
    if (Y[1] == Y[Count]) {
      Result = X[1]; // Flat curve
    } else if (Y[1] < Y[Count]) { // Assume curve going up monotonically
      if (theY <= Y[1]) Result = X[1]; // this shouldn't really happen
      else if (theY >= Y[Count]) {
        I = Count;
        if ((ExtrapolatePastEnd) && (I > 1)) {
          if (Y[I] > Y[I-1])
            Result = X[I] + (X[I] - X[I-1]) / (Y[I] - Y[I-1]) * (theY - Y[I]);
          else Result = X[I];
        } else
          Result = X[Count];
      } else {
        I = 2;
        while ((I < Count) && (Y[I] <= theY)) 
          I = I + 1;
          if (Y[I] > Y[I-1])
            Result = X[I-1] + (X[I] - X[I-1]) / (Y[I] - Y[I-1]) * (theY - Y[I-1]);
          else 
            Result = X[I]; // Could occur if curve is flat
      }
    } else if (Y[1] > Y[Count]) {// Assume curve going down monotonically
      if (theY >= Y[1]) Result = X[1];
      else if (theY <= Y[Count]) Result = X[Count];
      else {
        I = 2;
        while ((I < Count) && (Y[I] >= theY)) 
          I = I + 1; // Get first Y < theY
        if (Y[I] < Y[I-1]) 
          Result = X[I-1] + (X[I] - X[I-1]) / (Y[I] - Y[I-1]) * (theY - Y[I-1]);
        else 
          Result = X[I]; // Could occur if curve is flat
      }
    } 
  }*/
  return Result;
}
