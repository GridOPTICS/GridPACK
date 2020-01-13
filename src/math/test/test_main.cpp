// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   test_main.cpp
 * @author William A. Perkins
 * @date   2016-12-16 09:10:59 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December 16, 2016 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------

#include "gridpack/utilities/exception.hpp"
#include "gridpack/parallel/parallel.hpp"
#include "gridpack/environment/environment.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "math.hpp"


#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

/// The configuration used for these tests
gridpack::utility::Configuration::CursorPtr test_config;

bool init_function()
{
  gridpack::parallel::Communicator world;
  boost::scoped_ptr<gridpack::utility::Configuration> 
    config(gridpack::utility::Configuration::configuration());
  
  config->enableLogging();
  config->open("gridpack.xml", world);

  test_config = config->getCursor("GridPACK.MathTests");

  return true;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::Environment env(argc, argv);
  gridpack::parallel::Communicator world;

  // In the GridPACK setup, CTest determines unit test success or
  // failure by a phrase produced by Boost::test. When run in
  // parallel, these phrases can get corrupted. The following is to
  // make sure at least one uncorrupted occurance of the phrase
  // occurs.

  int lresult = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  lresult = (lresult == boost::exit_success ? 0 : 1);
  int gresult;

  boost::mpi::all_reduce(world, lresult, gresult, std::plus<int>());
  if (world.rank() == 0) {
    if (gresult == 0) {
      std::cout << "No errors detected" << std::endl;
    } else {
      std::cout << "failure detected" << std::endl;
    }
  }
  return gresult;
}


