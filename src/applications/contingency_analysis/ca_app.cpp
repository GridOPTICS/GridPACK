/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ca_app.cpp
 * @author Yousu Chen 
 * @date   January 20, 2014
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#include "gridpack/math/matrix.hpp"
#include "gridpack/math/vector.hpp"
#include "gridpack/math/linear_solver.hpp"
#include "gridpack/math/linear_matrix_solver.hpp"
#include "gridpack/applications/contingency_analysis/ca_app.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/applications/contingency_analysis/ca_factory.hpp"

// Calling program for contingency analysis application

/**
 * Basic constructor
 */
gridpack::contingency_analysis::CAApp::CAApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::contingency_analysis::CAApp::~CAApp(void)
{
}

/**
 * Execute application
 */
void gridpack::contingency_analysis::CAApp::execute(int argc, char** argv)
{
  gridpack::parallel::Communicator world;
  boost::shared_ptr<CANetwork> network(new CANetwork(world));

  // read configuration file
  gridpack::utility::Configuration *config = gridpack::utility::Configuration::configuration();
  if (argc >= 2 && argv[1] != NULL) {
    char inputfile[256];
    sprintf(inputfile,"%s",argv[1]);
    config->open(inputfile,world);
  } else {
    config->open("input.xml",world);
  }
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Contingency_analysis");
  std::string filename = cursor->get("networkConfiguration",
      "No network configuration specified");

  // load input file
  gridpack::parser::PTI23_parser<CANetwork> parser(network);
  parser.parse(filename.c_str());

  // partition network
  network->partition();

  // Create serial IO object to export data from buses or branches
  gridpack::serial_io::SerialBusIO<CANetwork> busIO(128, network);
  gridpack::serial_io::SerialBranchIO<CANetwork> branchIO(128, network);
  char ioBuf[128];

  // create factory
  gridpack::contingency_analysis::CAFactory factory(network);
  factory.load();

  // set network components using factory
  factory.setComponents();

  // set YBus components so that you can create Y matrix  
  factory.setYBus();

  factory.setMode(YBUS);
  gridpack::mapper::FullMatrixMap<CANetwork> ybusMap(network);
  boost::shared_ptr<gridpack::math::Matrix> orgYbus = ybusMap.mapToMatrix();
  branchIO.header("\n=== orginal ybus: ============\n");
  orgYbus->print();
}

