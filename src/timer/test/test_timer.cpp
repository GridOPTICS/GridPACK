/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

#include <math.h>

#include "mpi.h"
#include "gridpack/timer/coarse_timer.hpp"

#define LOOPSIZE 10000

/**
 * This program is designed to perform some simple tests of the CoarseTimer
 * class. It times a loop and tests some failure modes
 */
void run()
{
  int i;
  double t;

  // Grab an instance of timer
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();

  // Create three timer categories
  int t_success = timer->createCategory("Succeed");
  int t_fail1 = timer->createCategory("Failure 1");
  int t_fail2 = timer->createCategory("Failuer 2");
  timer->start(t_success);
  timer->start(t_fail1);
  for (i=0; i<LOOPSIZE; i++) {
    if (i>0) {
      t = exp(1.0/static_cast<double>(i));
    } else {
      t = exp(1.0);
    }
  }
  timer->stop(t_success);
  timer->stop(t_fail2);
  timer->dump();
}

main (int argc, char **argv) {

  // Initialize MPI libraries
  int ierr = MPI_Init(&argc, &argv);
  int me;
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &me);
  if (me == 0) {
    printf("Testing Mapper Module\n");
  }

  run();

  ierr = MPI_Finalize();
}

