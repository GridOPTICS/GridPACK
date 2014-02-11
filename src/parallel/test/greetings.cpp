/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   greetings.cpp
 * @author William A. Perkins
 * @date   2014-02-10 14:39:36 d3g096
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

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);

  // Limit scope so that program exits cleanly
  if (1) {
    gridpack::parallel::Communicator world;

    int nprocs = world.size();
    int me = world.rank();
    int i;
    char buf[128];

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
      int nsize = gcomm.size();
      for (i=0; i<nprocs; i++) {
        if (i == me) 
          printf("(assignment) I am process %d (original process is %d) of %d.\n",
              gcomm.rank(),me,gcomm.size());
      }
    }

    gridpack::parallel::Communicator cpcomm(gcomm);
    if (me == 0) printf("\nTesting copy constructor for communicators:\n\n");
    if (nprocs > 1) {
      int nsize = gcomm.size();
      for (i=0; i<nprocs; i++) {
        if (i == me) 
          printf("(copy) I am process %d (original process is %d) of %d.\n",
              gcomm.rank(),me,gcomm.size());
      }
    }
  }
  return 0;
}

