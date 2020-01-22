/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   random_test.cpp
 * @author Bruce Palmer
 * @date   March 3, 2015
 * 
 * @brief  A simple test of the GridPACK random number generator wrappers
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------

#include <iostream>
#include <ga.h>
#include "gridpack/parallel/parallel.hpp"
#include "gridpack/parallel/random.hpp"
#include "gridpack/environment/environment.hpp"

#define MAX_VALS  1000000

#define MAX_BINS  100

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::Environment env(argc, argv);
  GA_Initialize();
  // Create an artificial scope so that all objects call their destructors
  // before GA_Terminate is called
  if (1) {

    int i;
    int inormal[MAX_BINS];
    double normal[MAX_BINS];
    int igaussian[MAX_BINS];
    double gaussian[MAX_BINS];
    for (i=0; i<MAX_BINS; i++) {
      inormal[i] = 0;
      igaussian[i] = 0;
    }
    double ninc = 1.0/static_cast<double>(MAX_BINS);
    double ginc = 6.0/static_cast<double>(MAX_BINS);
    double offset = 3.0;

    gridpack::random::Random random;

    int iseed = 993820;
    random.seed(iseed);
    double x;
    int idx;
    for (i=0; i<MAX_VALS; i++) {
      x = random.drand();
      idx = static_cast<int>(x/ninc);
      if (idx < MAX_BINS) inormal[idx]++;
      x = random.grand()+offset;
      idx = static_cast<int>(x/ginc);
      if (idx >= 0 && idx < MAX_BINS) igaussian[idx]++;
    }
    i = MAX_BINS;
    GA_Igop(inormal,i,"+");
    GA_Igop(igaussian,i,"+");

    //normalize distributions
    double rnorm = 0.0;
    for (i=0; i<MAX_BINS; i++) {
      rnorm += static_cast<double>(inormal[i]);
    }
    for (i=0; i<MAX_BINS; i++) {
      normal[i] = static_cast<double>(inormal[i])/(rnorm*ninc);
    }

    rnorm = 0.0;
    for (i=0; i<MAX_BINS; i++) {
      rnorm += static_cast<double>(igaussian[i]);
    }
    for (i=0; i<MAX_BINS; i++) {
      gaussian[i] = static_cast<double>(igaussian[i])/(rnorm*ginc);
    }

    //print out distributions
    
    if (GA_Nodeid() == 0) {
      printf("Print out results for uniform distribution\n\n");
      for (i=0; i<MAX_BINS; i++) {
        printf("%12.4f   %12.6f\n",ninc*static_cast<double>(i),normal[i]);
      }
      printf("\nPrint out results for normal distribution\n\n");
      for (i=0; i<MAX_BINS; i++) {
        printf("%12.4f,%12.6f\n",ginc*static_cast<double>(i)-offset,
            gaussian[i]);
      }
    }
  }
  return 0;
}

