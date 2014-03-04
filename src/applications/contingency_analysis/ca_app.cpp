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
 * @param comm communicator that hosts contingency application
 */
gridpack::contingency_analysis::CAApp::CAApp(gridpack::parallel::Communicator comm)
{
  p_network.reset(new CANetwork(comm));
  p_factory.reset(new CAFactory(p_network));
}

/**
 * Basic destructor
 */
gridpack::contingency_analysis::CAApp::~CAApp(void)
{
}

/**
 * Initialize application by reading in grid network, partioning it and
 * setting up buffers and indices
 * @param argc number of arguments
 * @param argv list of character strings
 */ 
void gridpack::contingency_analysis::CAApp::init(int argc, char** argv)
{
  //boost::shared_ptr<CANetwork> p_network(new CANetwork(comm));
  gridpack::parallel::Communicator comm = p_network->communicator();

  // read configuration file
  gridpack::utility::Configuration *config = gridpack::utility::Configuration::configuration();
  if (argc >= 2 && argv[1] != NULL) {
    char inputfile[256];
    sprintf(inputfile,"%s",argv[1]);
    config->open(inputfile,comm);
  } else {
    config->open("input.xml",comm);
  }
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Contingency_analysis");
  std::string filename = cursor->get("networkConfiguration",
      "No network configuration specified");

  // load input file
  gridpack::parser::PTI23_parser<CANetwork> parser(p_network);
  parser.parse(filename.c_str());

  // partition network
  p_network->partition();

  // create factory
  //boost::shared_ptr<CAFactory> p_factory(new CAFactory(p_network));
  p_factory->load();

  // set network components using factory
  p_factory->setComponents();
 
  p_factory->setExchange();

  // initialize bus data exchange
  p_network->initBusUpdate();

}
/**
 * Execute application
 */
void gridpack::contingency_analysis::CAApp::execute(
    gridpack::contingency_analysis::Contingency contingency)
{
  gridpack::parallel::Communicator comm = p_network->communicator();
  gridpack::utility::CoarseTimer *timer = 
    gridpack::utility::CoarseTimer::instance();

  // get configuration file
  gridpack::utility::Configuration *config = gridpack::utility::Configuration::configuration();
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Contingency_analysis");
  std::string filename = cursor->get("networkConfiguration",
      "No network configuration specified");

  //Set voltage and phase angle to initial values
  p_factory->resetVoltage();

  // Create serial IO object to export data from buses or branches
  gridpack::serial_io::SerialBusIO<CANetwork> busIO(128, p_network);
  gridpack::serial_io::SerialBranchIO<CANetwork> branchIO(128, p_network);
  char ioBuf[128];
  sprintf(ioBuf,"%s.out",contingency.p_name.c_str());
  busIO.open(ioBuf);
  branchIO.setStream(busIO.getStream());

  // set contingency
  p_factory->setContingency(contingency);

  // set YBus components so that you can create Y matrix  
  p_factory->setYBus();

  p_factory->setMode(gridpack::powerflow::YBus);
  gridpack::mapper::FullMatrixMap<CANetwork> ybusMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> orgYbus = ybusMap.mapToMatrix();
  //branchIO.header("\n=== orginal ybus: ============\n");
  //orgYbus->print();

  //////////////////////////////////////////////////////////////
  p_factory->setMode(gridpack::powerflow::S_Cal);

  // make Sbus components to create S vector
  p_factory->setSBus();

  // Set PQ
  p_factory->setMode(gridpack::powerflow::RHS);
  gridpack::mapper::BusVectorMap<CANetwork> vMap(p_network);
  boost::shared_ptr<gridpack::math::Vector> PQ = vMap.mapToVector();
  branchIO.header("\n=== PQ: ============\n");
  PQ->print();
  p_factory->setMode(gridpack::powerflow::Jacobian);
  gridpack::mapper::FullMatrixMap<CANetwork> jMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> J = jMap.mapToMatrix();
  //J->print(); 

  // Create X vector by cloning PQ
  boost::shared_ptr<gridpack::math::Vector> X(PQ->clone());

  // Convergence and iteration parameters
  double tolerance;
  int max_iteration;
  ComplexType tol;

  // These need to eventually be set using configuration file
  tolerance = 1.0e-6;
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
    p_factory->setMode(gridpack::powerflow::RHS);
    vMap.mapToBus(X);

    // Exchange data between ghost buses (I don't think we need to exchange data
    // between branches)
    p_network->updateBuses();

    // Create new versions of Jacobian and PQ vector
    vMap.mapToVector(PQ);
    p_factory->setMode(gridpack::powerflow::Jacobian);
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
  p_factory->setMode(gridpack::powerflow::RHS);
  vMap.mapToBus(X);

  branchIO.header("\n   Branch Power Flow\n");
  branchIO.header("\n        Bus 1       Bus 2            P"
                  "                    Q\n");
  branchIO.write();


  busIO.header("\n   Bus Voltages and Phase Angles\n");
  busIO.header("\n   Bus Number      Phase Angle      Voltage Magnitude\n");
  busIO.write();

  busIO.close();
  // clear contingency
  p_factory->clearContingency(contingency);
}

