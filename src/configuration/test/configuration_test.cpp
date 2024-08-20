/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   configuration_test.cpp
 * @author William A. Perkins
 * @date   2024-08-20 06:36:13 d3g096
 * 
 * @brief  A test of Configurable and Configuration
 * 
 * 
 */

#include <iostream>
#include <iterator>
#include <memory>
#include <boost/lexical_cast.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>
#include <boost/serialization/string.hpp>

#include "gridpack/parallel/parallel.hpp"
#include "gridpack/environment/environment.hpp"
#include "gridpack/utilities/uncopyable.hpp"
#include "configurable.hpp"


BOOST_AUTO_TEST_SUITE ( ConfigurationTest ) 

// -------------------------------------------------------------
//  class ConfigurableThing
// -------------------------------------------------------------
class ConfigurableThing 
  : public gridpack::utility::Configurable,
    private gridpack::utility::Uncopyable
{
public:

  /// Default constructor.
  ConfigurableThing(void)
    : gridpack::utility::Configurable("Thing"), 
      gridpack::utility::Uncopyable(),
      string1(), string2(), real(0.0), integer(0), flag(true)
  {}

  /// Destructor
  ~ConfigurableThing(void)
  {}

  std::string string1;
  std::string string2;
  double real;
  int integer;
  bool flag;

protected:
  
  void 
  p_configure(gridpack::utility::Configuration::CursorPtr props)
  {
    if (props != NULL) {
      string1 = props->get("String1", string1);
      string2 = props->get("String2", string2);
      real = props->get("Float", real);
      integer = props->get("Integer", integer);
      flag = props->get("Flag", flag);
    }
  }
};


BOOST_AUTO_TEST_CASE( Configurable )
{
  std::unique_ptr<gridpack::utility::Configuration> 
    config(gridpack::utility::Configuration::configuration());

  BOOST_REQUIRE(config.get() != NULL);

  gridpack::utility::Configuration::CursorPtr cursor =
       config->getCursor("GridPACK");

  BOOST_REQUIRE(cursor != NULL);

  std::unique_ptr<ConfigurableThing> thing(new ConfigurableThing);

  thing->configure(cursor);

  std::cout << "string1: \"" << thing->string1 << "\"" << std::endl;
  std::cout << "string2: \"" << thing->string2 << "\"" << std::endl;
  std::cout << "   real: " << thing->real << std::endl;
  std::cout << "integer: " << thing->integer << std::endl;
  std::cout << "   flag: " << thing->flag << std::endl;

  BOOST_CHECK_EQUAL(thing->string1, "A simple string.");
  BOOST_CHECK_EQUAL(thing->integer, 123);
  BOOST_CHECK(!thing->flag);
}

BOOST_AUTO_TEST_SUITE_END( )

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
  gridpack::Environment env(argc, argv);
  gridpack::parallel::Communicator world;

  gridpack::utility::Configuration *config =
    gridpack::utility::Configuration::configuration();
  
  config->enableLogging();
  config->open("configuration_test.xml", world);

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
