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

#define LOOPSIZE 10000

BOOST_AUTO_TEST_SUITE ( CoarseTimerTest )
/**
 * This program is designed to perform some simple tests of the CoarseTimer
 * class. It times a loop and tests some failure modes
 */

BOOST_AUTO_TEST_CASE( Timings )
{
  int i;
  double t;

  // Grab an instance of timer
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  BOOST_REQUIRE(timer != NULL);

  // Create three timer categories
  int t_success = timer->createCategory("Succeed");
  BOOST_CHECK_EQUAL(t_success, 0);
  int t_fail1 = timer->createCategory("Failure 1");
  BOOST_CHECK_EQUAL(t_fail1, 1);
  int t_fail2 = timer->createCategory("Failuer 2");
  BOOST_CHECK_EQUAL(t_fail2, 2);
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

BOOST_AUTO_TEST_SUITE_END( )

bool init_function(void)
{
  return true;
}

// Main program
int main (int argc, char **argv) {

  // Initialize parallel environment
  gridpack::parallel::Environment env(argc, argv);
  boost::mpi::communicator world;
  int me = world.rank();
  if (me == 0) {
    printf("Testing Mapper Module\n");
  }

  return ::boost::unit_test::unit_test_main( &init_function, argc, argv );
}
