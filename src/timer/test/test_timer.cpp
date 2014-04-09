/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

#include <math.h>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include "mpi.h"
#include "gridpack/parallel/distributed.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/timer/local_timer.hpp"

#define LOOPSIZE 1000000

BOOST_AUTO_TEST_SUITE ( CoarseTimerTest )
/**
 * This program is designed to perform some simple tests of the CoarseTimer
 * class. It times a loop and tests some failure modes
 */

BOOST_AUTO_TEST_CASE( Timings )
{
  int i;
  double t;
  int me, nprocs;
  int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &me);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // Grab an instance of timer
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  BOOST_REQUIRE(timer != NULL);

  // Create three timer categories
  int t_success = timer->createCategory("CoarseTimer: Succeed");
  BOOST_CHECK_EQUAL(t_success, 0);
  int t_fail1 = timer->createCategory("CoarseTimer: Failure 1");
  BOOST_CHECK_EQUAL(t_fail1, 1);
  int t_fail2 = timer->createCategory("CoarseTimer: Failure 2");
  BOOST_CHECK_EQUAL(t_fail2, 2);
  timer->start(t_success);
  timer->start(t_fail1);
  double t_chk = MPI_Wtime();
  int nloop = LOOPSIZE;
  nloop += static_cast<int>(static_cast<double>(nloop*me)
        / static_cast<double>(nprocs));
  for (i=0; i<nloop; i++) {
    if (i>0) {
      t = exp(1.0/static_cast<double>(i));
    } else {
      t = exp(1.0);
    }
  }
  t_chk = MPI_Wtime()-t_chk;
  timer->stop(t_success);
  timer->stop(t_fail2);
  double tavg,tmax,tmin;
  ierr=MPI_Allreduce(&t_chk, &tavg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  tavg = tavg/static_cast<double>(nprocs);
  ierr=MPI_Allreduce(&t_chk, &tmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  ierr=MPI_Allreduce(&t_chk, &tmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  if (me == 0) {
    printf("Average time from MPI_Wtime(): %12.6f\n",tavg);
    printf("Maximum time from MPI_Wtime(): %12.6f\n",tmax);
    printf("Minimum time from MPI_Wtime(): %12.6f\n\n",tmin);
  }
  timer->dump();

  gridpack::parallel::Communicator world;
  gridpack::utility::LocalTimer ltime(world);

  // Create three local timer categories
  t_success = ltime.createCategory("LocalTimer: Succeed");
  BOOST_CHECK_EQUAL(t_success, 0);
  t_fail1 = ltime.createCategory("LocalTimer: Failure 1");
  BOOST_CHECK_EQUAL(t_fail1, 1);
  t_fail2 = ltime.createCategory("LocalTimer: Failure 2");
  BOOST_CHECK_EQUAL(t_fail2, 2);
  ltime.start(t_success);
  ltime.start(t_fail1);
  for (i=0; i<nloop; i++) {
    if (i>0) {
      t = exp(1.0/static_cast<double>(i));
    } else {
      t = exp(1.0);
    }
  }
  ltime.stop(t_success);
  ltime.stop(t_fail2);
  ltime.dump();

}

BOOST_AUTO_TEST_SUITE_END( )

bool init_function(void)
{
  return true;
}

// Main program
int main (int argc, char **argv) {

  // Initialize parallel environment
  gridpack::parallel::Environment env(argc, argv);
  gridpack::parallel::Communicator world;
  int me = world.rank();
  if (me == 0) {
    printf("Testing Mapper Module\n");
  }

  return ::boost::unit_test::unit_test_main( &init_function, argc, argv );
}
