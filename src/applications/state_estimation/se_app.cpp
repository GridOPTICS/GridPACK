/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_app.cpp
 * @author Yousu Chen 
 * @date   2/24/2014 
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
#include "gridpack/applications/state_estimation/se_app.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/applications/state_estimation/se_factory.hpp"

// Calling program for state estimation application

/**
 * Basic constructor
 */
gridpack::state_estimation::SEApp::SEApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::state_estimation::SEApp::~SEApp(void)
{
}

/**
 * Execute application
 */
void gridpack::state_estimation::SEApp::execute(int argc, char** argv)
{
  gridpack::parallel::Communicator world;
  boost::shared_ptr<SENetwork> network(new SENetwork(world));

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
  cursor = config->getCursor("Configuration.State_estimation");
  std::string filename;
  if (!cursor->get("networkConfiguration",&filename)) {
     printf("No network configuration specified\n");
     return;
  }

  // load input file
  gridpack::parser::PTI23_parser<SENetwork> parser(network);
  parser.parse(filename.c_str());

  // partition network
  network->partition();

  // Create serial IO object to export data from buses or branches
  gridpack::serial_io::SerialBusIO<SENetwork> busIO(128, network);
  gridpack::serial_io::SerialBranchIO<SENetwork> branchIO(128, network);
  char ioBuf[128];

  // create factory
  gridpack::state_estimation::SEFactory factory(network);
  factory.load();

  // set network components using factory
  factory.setComponents();

  // Set up bus data exchange buffers. Need to decide what data needs to be exchanged
  factory.setExchange();

  // Create bus data exchange
  network->initBusUpdate();

  // set YBus components so that you can create Y matrix  
  factory.setYBus();

  factory.setMode(YBUS);
  gridpack::mapper::FullMatrixMap<SENetwork> ybusMap(network);
  boost::shared_ptr<gridpack::math::Matrix> orgYbus = ybusMap.mapToMatrix();
  ///branchIO.header("\n=== orginal ybus: ============\n");
  ///orgYbus->print();

}

