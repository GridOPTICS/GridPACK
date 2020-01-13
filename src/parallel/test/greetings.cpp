/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   greetings.cpp
 * @author William A. Perkins
 * @date   2014-12-09 09:39:25 d3g096
 * 
 * @brief  A simple test of the GridPACK parallel environment
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created May  6, 2013 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------

#include <iostream>
#include <ga.h>
#include "gridpack/parallel/parallel.hpp"
#include "gridpack/environment/environment.hpp"

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::Environment env(argc, argv);

  // Limit scope so that program exits cleanly
  if (1) {
    gridpack::parallel::Communicator world;

    int nprocs = world.size();
    int me = world.rank();
    int i;

    for (i=0; i<nprocs; i++) {
      if (i == me) {
        printf("I am process %d of %d.\n",me,nprocs);
      }
    }
    world.barrier();

    // Create two communicators, each containing half the processors
    if (me == 0) printf("\nCreating communicators using split:\n\n");
    if (nprocs > 1) {
      int color;
      if (static_cast<double>(me)/static_cast<double>(nprocs) < 0.5) {
        color = 0;
      } else {
        color = 1;
      }
      gridpack::parallel::Communicator lcomm = world.split(color);
      for (i=0; i<nprocs; i++) {
        if (i == me) 
          printf("I am process %d (original process is %d) of %d on communicator %d.\n",
              lcomm.rank(),me,lcomm.size(),color);
      }
    }

    // Create a communictor to test copy constructor
    gridpack::parallel::Communicator gcomm;
    if (me == 0) printf("\nCreating communicators using divide:\n\n");
    if (nprocs > 1) {
      int nsize = nprocs/2;
      gridpack::parallel::Communicator lcomm = world.divide(nsize);
      gcomm = lcomm;
      for (i=0; i<nprocs; i++) {
        if (i == me) 
          printf("I am process %d (original process is %d) of %d.\n",
              lcomm.rank(),me,lcomm.size());
      }
    }

    if (me == 0) printf("\nTesting assignment operator for communicators:\n\n");
    if (nprocs > 1) {
      for (i=0; i<nprocs; i++) {
        if (i == me) 
          printf("(assignment) I am process %d (original process is %d) of %d.\n",
              gcomm.rank(),me,gcomm.size());
      }
    }

    gridpack::parallel::Communicator cpcomm(gcomm);
    if (me == 0) printf("\nTesting copy constructor for communicators:\n\n");
    if (nprocs > 1) {
      for (i=0; i<nprocs; i++) {
        if (i == me) 
          printf("(copy) I am process %d (original process is %d) of %d.\n",
              gcomm.rank(),me,gcomm.size());
      }
    }
    // Test sum operations
#define NVALS 10
    if (me == 0) printf("\nTesting summation operator for communicators\n\n");
    int nvals = NVALS;
    float xf[NVALS];
    int ok = 0;
    for (i=0; i<nvals; i++) {
      xf[i] = static_cast<float>(i);
    }
    world.sum(xf,nvals);
    for (i=0; i<nvals; i++) {
      if (xf[i] != static_cast<float>(i*nprocs)) {
        printf("(sum) failure for float value %d\n",i);
        ok = 1;
      }
    }
    double xd[NVALS];
    for (i=0; i<nvals; i++) {
      xd[i] = static_cast<double>(i);
    }
    world.sum(xd,nvals);
    for (i=0; i<nvals; i++) {
      if (xd[i] != static_cast<double>(i*nprocs)) {
        printf("(sum) failure for double value %d\n",i);
        ok = 1;
      }
    }
    int xi[NVALS];
    for (i=0; i<nvals; i++) {
      xi[i] = i;
    }
    world.sum(xi,nvals);
    for (i=0; i<nvals; i++) {
      if (xi[i] != i*nprocs) {
        printf("(sum) failure for int value %d\n",i);
        ok = 1;
      }
    }
    long xl[NVALS];
    for (i=0; i<nvals; i++) {
      xl[i] = static_cast<long>(i);
    }
    world.sum(xl,nvals);
    for (i=0; i<nvals; i++) {
      if (xl[i] != static_cast<long>(i*nprocs)) {
        printf("(sum) failure for long value %d\n",i);
        ok = 1;
      }
    }
    gridpack::ComplexType xc[NVALS];
    for (i=0; i<nvals; i++) {
      xc[i] = gridpack::ComplexType(static_cast<double>(i),
          static_cast<double>(i));
    }
    world.sum(xc,nvals);
    for (i=0; i<nvals; i++) {
      if (xc[i] != gridpack::ComplexType(static_cast<double>(i*nprocs),
            static_cast<double>(i*nprocs))) {
        printf("(sum) failure for complex value %d\n",i);
        ok = 1;
      }
    }
    if (me == 0 && ok == 0) {
      printf("\nSummation test for communicators passed\n\n");
    } else if (me == 0) {
      printf("\nSummation test for communicators failed\n\n");
    }
  }

  return 0;
}

