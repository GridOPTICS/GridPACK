/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   PTI23_test.cpp
 * @author Kevin Glass
 * @date   2014-02-11 10:28:46 d3g096
 * 
 * @brief  Test PTI23_parser capability. Currently not implemented.
 * 
 * 
 */

#include <iostream>
#include <string>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <gridpack/parser/PTI23_parser.hpp>
#include <gridpack/configuration/configuration.hpp>
#include <gridpack/timer/coarse_timer.hpp>

#include "mpi.h"
#include <macdecls.h>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#define EPSILON     0.0000001
#define TOLERANCE(x, y, eps) (x - EPSILON > y && x + EPSILON < y)

BOOST_AUTO_TEST_SUITE(Parser)

BOOST_AUTO_TEST_CASE( SerialInputTest )
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  std::string                line;

  // read configuration file
  gridpack::utility::Configuration *config =
    gridpack::utility::Configuration::configuration();
  gridpack::parallel::Communicator world;
  config->open("test.xml",world);
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.TestParser");
  std::string filename = cursor->get("networkConfiguration",
      "No network configuration specified");

  int me;
  int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &me);
  int t_read = timer->createCategory("Serial Read");
  timer->start(t_read);
  if (me == 0) {
    std::ifstream            input;
    input.open(filename.c_str());
    if (!input.is_open()) {
      throw gridpack::Exception("failed to open case data file");
    }
    while(std::getline(input,line)) {
    }
    input.close();
  }
  timer->stop(t_read);
  timer->dump();
}

BOOST_AUTO_TEST_SUITE_END()

// -------------------------------------------------------------
// init_function
// -------------------------------------------------------------
bool init_function()
{
  return true;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  gridpack::parallel::Communicator world;
  int me = world.rank();
  if (me == 0) {
    printf("Testing Serial Input\n");
  }
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  return result;
}


