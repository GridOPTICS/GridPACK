/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   task_test.cpp
 * @author Bruce Palmer
 * @date   February 10, 2014
 * 
 * @brief  A simple test of the GridPACK task manager
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created February 10, 2014
// Last Change: February 10, 2014
// -------------------------------------------------------------

#include <iostream>
#include <ga.h>
#include "gridpack/parallel/parallel.hpp"
#include "gridpack/parallel/task_manager.hpp"

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  GA_Initialize();
  // Create an artificial scope so that all objects call their destructors
  // before GA_Terminate is called
  if (1) {
    gridpack::parallel::Communicator world;

    int nprocs = world.size();
    int me = world.rank();
    int i;
    char buf[128];

    int ntasks = 4*nprocs; 

    gridpack::parallel::TaskManager tskmgr;

    tskmgr.set(ntasks);
    int itask;
    while(tskmgr.nextTask(&itask)) {
      printf("Evaluating task %d on processor %d out of %d\n",itask,me,nprocs);
    }

    // Create two communicators, each containing half the processors
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

      tskmgr.set(ntasks);
      while(tskmgr.nextTask(lcomm,&itask)) {
        printf("Evaluating task %d on processor %d (global id %d) in sub-communicator of size %d\n",
            itask,lcomm.rank(),me,lcomm.size());
      }
    }
  }

  GA_Terminate();
  return 0;
}

