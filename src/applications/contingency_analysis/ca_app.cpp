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

#include <iostream>
#include "gridpack/math/matrix.hpp"
#include "gridpack/math/vector.hpp"
#include "gridpack/math/linear_solver.hpp"
#include "gridpack/math/newton_raphson_solver.hpp"
#include "gridpack/math/nonlinear_solver.hpp"
#include "ca_app.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "ca_factory.hpp"
#include "gridpack/timer/coarse_timer.hpp"

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
void gridpack::contingency_analysis::CAApp::execute(
    gridpack::parallel::Communicator comm,
    gridpack::contingency_analysis::Contingency contingency,
    int argc, char** argv)
{
  gridpack::utility::CoarseTimer *timer = 
    gridpack::utility::CoarseTimer::instance();
  boost::shared_ptr<CANetwork> network(new CANetwork(comm));

  // read configuration file
  gridpack::utility::Configuration *config = gridpack::utility::Configuration::configuration();
  if (argc >= 2 && argv[1] != NULL) {
    char inputfile[256];
    sprintf(inputfile,"%s",argv[1]);
    config->open(inputfile,comm);
  } else {
    config->open("input.xml",comm);
  }
  printf("Got to 1\n");
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Contingency_analysis");
  std::string filename = cursor->get("networkConfiguration",
      "No network configuration specified");

  // load input file
  printf("Got to 2\n");
  gridpack::parser::PTI23_parser<CANetwork> parser(network);
  parser.parse(filename.c_str());
  printf("Got to 3\n");

  // partition network
  network->partition();
  printf("Got to 4\n");

  // Create serial IO object to export data from buses or branches
  gridpack::serial_io::SerialBusIO<CANetwork> busIO(128, network);
  gridpack::serial_io::SerialBranchIO<CANetwork> branchIO(128, network);
  char ioBuf[128];
  printf("Got to 5\n");

  // create factory
  gridpack::contingency_analysis::CAFactory factory(network);
  factory.load();
  printf("Got to 6\n");

  // set network components using factory
  factory.setComponents();
  printf("Got to 7\n");
 
  factory.setExchange();
  printf("Got to 8\n");

  // initialize bus data exchange
  network->initBusUpdate();
  printf("Got to 9\n");

  // set YBus components so that you can create Y matrix  
  factory.setYBus();
  printf("Got to 10\n");

  factory.setMode(YBus);
  gridpack::mapper::FullMatrixMap<CANetwork> ybusMap(network);
  printf("Got to 11\n");
  boost::shared_ptr<gridpack::math::Matrix> orgYbus = ybusMap.mapToMatrix();
  printf("Got to 12\n");
  branchIO.header("\n=== orginal ybus: ============\n");
  orgYbus->print();
  printf("Got to 13\n");

  //////////////////////////////////////////////////////////////
  factory.setMode(S_Cal);

  // make Sbus components to create S vector
  factory.setSBus();

  // Set PQ
  factory.setMode(RHS);
  gridpack::mapper::BusVectorMap<CANetwork> vMap(network);
  boost::shared_ptr<gridpack::math::Vector> PQ = vMap.mapToVector();
  //PQ->print();
  factory.setMode(Jacobian);
  gridpack::mapper::FullMatrixMap<CANetwork> jMap(network);
  boost::shared_ptr<gridpack::math::Matrix> J = jMap.mapToMatrix();
  //J->print(); 

  // Create X vector by cloning PQ
  boost::shared_ptr<gridpack::math::Vector> X(PQ->clone());

  // Convergence and iteration parameters
  double tolerance;
  int max_iteration;
  ComplexType tol;

  // These need to eventually be set using configuration file
  tolerance = 1.0e-5;
  max_iteration = 50;

  // Create linear solver
  gridpack::math::LinearSolver solver(*J);
  solver.configure(cursor);

  tol = 2.0*tolerance;
  int iter = 0;

  // First iteration
  X->zero(); //might not need to do this
  solver.solve(*PQ, *X);
  tol = PQ->normInfinity();
  //J->print();
  //PQ->print();
  //X->print();

  while (real(tol) > tolerance && iter < max_iteration) {
    // Push current values in X vector back into network components
    // Need to implement setValues method in PFBus class in order for this to
    // work
    factory.setMode(RHS);
    vMap.mapToBus(X);

    // Exchange data between ghost buses (I don't think we need to exchange data
    // between branches)
    network->updateBuses();

    // Create new versions of Jacobian and PQ vector
    vMap.mapToVector(PQ);
    factory.setMode(Jacobian);
    jMap.mapToMatrix(J);

    // Create linear solver
    X->zero(); //might not need to do this
#if 1
//    sprintf(dbgfile,"j%d.bin",iter+1);
//    J->saveBinary(dbgfile);
//    sprintf(dbgfile,"pq%d.bin",iter+1);
//    PQ->saveBinary(dbgfile);
    solver.solve(*PQ, *X);
#else
//    sprintf(dbgfile,"j%d.bin",iter+1);
//    J->saveBinary(dbgfile);
//    sprintf(dbgfile,"pq%d.bin",iter+1);
//    PQ->saveBinary(dbgfile);
    gridpack::math::LinearSolver isolver(*J);
    isolver.configure(cursor);
    isolver.solve(*PQ, *X);
#endif

    tol = PQ->normInfinity();
    sprintf(ioBuf,"\nIteration %d Tol: %12.6e\n",iter+1,real(tol));
    busIO.header(ioBuf);
    iter++;
  }
  // Push final result back onto buses
  factory.setMode(RHS);
  vMap.mapToBus(X);

  branchIO.header("\n   Branch Power Flow\n");
  branchIO.header("\n        Bus 1       Bus 2            P"
                  "                    Q\n");
  branchIO.write();


  busIO.header("\n   Bus Voltages and Phase Angles\n");
  busIO.header("\n   Bus Number      Phase Angle      Voltage Magnitude\n");
  busIO.write();
}

