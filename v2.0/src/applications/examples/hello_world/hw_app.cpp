/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hw_app.cpp
 * @author Bruce Palmer
 * @date   2014-01-28 10:31:00 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "hw_app.hpp"
#include "hw_factory.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/serial_io/serial_io.hpp"


// Calling program for hello world application

/**
 * Basic constructor
 */
gridpack::hello_world::HWApp::HWApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::hello_world::HWApp::~HWApp(void)
{
}

/**
 * Execute application
 * @param argc number of arguments
 * @param argv list of character strings
 */
void gridpack::hello_world::HWApp::execute(int argc, char** argv)
{
  // load input file
  gridpack::parallel::Communicator world;
  boost::shared_ptr<HWNetwork> network(new HWNetwork(world));

  // read configuration file
  std::string filename = "10x10.raw";

  // Read in external PTI file with network configuration
  gridpack::parser::PTI23_parser<HWNetwork> parser(network);
  parser.parse(filename.c_str());

  // partition network
  network->partition();

  // create factory
  gridpack::hello_world::HWFactory factory(network);
  factory.load();

  // Create serial IO object to export data from buses
  gridpack::serial_io::SerialBusIO<HWNetwork> busIO(128,network);
  busIO.header("\nMessage from buses\n");
  busIO.write();


  // Create serial IO object to export data from branches
  gridpack::serial_io::SerialBranchIO<HWNetwork> branchIO(128,network);
  branchIO.header("\nMessage from branches\n");
  branchIO.write();
}
