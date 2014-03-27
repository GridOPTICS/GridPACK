/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   rg_app.cpp
 * @author Bruce Palmer
 * @date   2014-02-05 08:25:44 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include "gridpack/math/matrix.hpp"
#include "gridpack/math/vector.hpp"
#include "gridpack/math/linear_solver.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "rg_app.hpp"
#include "rg_factory.hpp"


// Calling program for resistor grid application

/**
 * Basic constructor
 */
gridpack::resistor_grid::RGApp::RGApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::resistor_grid::RGApp::~RGApp(void)
{
}

/**
 * Execute application
 * @param argc number of arguments
 * @param argv list of character strings
 */
void gridpack::resistor_grid::RGApp::execute(int argc, char** argv)
{
  // read configuration file
  gridpack::parallel::Communicator world;
  gridpack::utility::Configuration *config =
    gridpack::utility::Configuration::configuration();
  config->open("input.xml",world);
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.ResistorGrid");

  // create network and read in external PTI file with network configuration
  boost::shared_ptr<RGNetwork> network(new RGNetwork(world));
  gridpack::parser::PTI23_parser<RGNetwork> parser(network);
  std::string filename;
  if (!cursor->get("networkConfiguration",&filename)) {
    filename = "small.raw";
  }
  parser.parse(filename.c_str());

  // partition network
  network->partition();

  // create factory and load parameters from input file to network components
  gridpack::resistor_grid::RGFactory factory(network);
  factory.load();

  // set network components using factory
  factory.setComponents();

  // create mapper to generate voltage matrix
  gridpack::mapper::FullMatrixMap<RGNetwork> vMap(network);
  boost::shared_ptr<gridpack::math::Matrix> V = vMap.mapToMatrix();

  // create mapper to generate RHS vector
  gridpack::mapper::BusVectorMap<RGNetwork> rMap(network);
  boost::shared_ptr<gridpack::math::Vector> R = rMap.mapToVector();

  // create solution vector by cloning R
  boost::shared_ptr<gridpack::math::Vector> X(R->clone());

  // create linear solver and solve equations
  gridpack::math::LinearSolver solver(*V);
  solver.configure(cursor);
  solver.solve(*R, *X);

  // push solution back on to buses
  rMap.mapToBus(X);

  // Create serial IO object to export data from buses
  gridpack::serial_io::SerialBusIO<RGNetwork> busIO(128,network);
  char ioBuf[128];

  busIO.header("\nVoltages on buses\n\n");
  busIO.write();


  gridpack::serial_io::SerialBranchIO<RGNetwork> branchIO(128,network);
  branchIO.header("\nCurrent on branches\n\n");
  branchIO.write();
}
