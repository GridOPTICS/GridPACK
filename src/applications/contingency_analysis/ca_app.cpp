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
#include "gridpack/timer/local_timer.hpp"

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
  gridpack::utility::CoarseTimer *timer = 
    gridpack::utility::CoarseTimer::instance();
  int t_parse = timer->createCategory("Parse Input");
  timer->start(t_parse);
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
  timer->stop(t_parse);

  // partition network
  int t_part = timer->createCategory("Partition Network");
  timer->start(t_part);
  p_network->partition();
  timer->stop(t_part);

  // create factory
  int t_init = timer->createCategory("Initialize Network");
  timer->start(t_init);
  p_factory->load();

  // set network components using factory
  p_factory->setComponents();
 
  p_factory->setExchange();

  // initialize bus data exchange
  p_network->initBusUpdate();
  timer->stop(t_init);
}
/**
 * Execute application
 */
void gridpack::contingency_analysis::CAApp::execute(
    gridpack::contingency_analysis::Contingency contingency)
{
  gridpack::parallel::Communicator comm = p_network->communicator();
  gridpack::utility::LocalTimer ltimer(comm);
  gridpack::utility::CoarseTimer *timer = 
    gridpack::utility::CoarseTimer::instance();
  int lt_task = ltimer.createCategory("Evaluate Contingency");
  int t_task = timer->createCategory("Evaluate Contingency");
  timer->start(t_task);
  ltimer.start(lt_task);

  // get configuration file
  gridpack::utility::Configuration *config = gridpack::utility::Configuration::configuration();
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Contingency_analysis");
  std::string filename = cursor->get("networkConfiguration",
      "No network configuration specified");
  double maxV = cursor->get("maxVoltage",1.1); 
  double minV = cursor->get("minVoltage",0.9); 
  // Convergence and iteration parameters
  double tolerance = cursor->get("tolerance",1.0e-6);
  int max_iteration = cursor->get("maxIteration",50);
  ComplexType tol;

  //Set voltage and phase angle to initial values
  p_factory->resetVoltage();

  // Create serial IO object to export data from buses or branches
  gridpack::serial_io::SerialBusIO<CANetwork> busIO(128, p_network);
  gridpack::serial_io::SerialBranchIO<CANetwork> branchIO(512, p_network);
  char ioBuf[512];
  sprintf(ioBuf,"%s.out",contingency.p_name.c_str());
  busIO.open(ioBuf);
  branchIO.setStream(busIO.getStream());

  sprintf(ioBuf,"\nMaximum voltage limit: %f\n",maxV);
  busIO.header(ioBuf);
  sprintf(ioBuf,"\nMinimum voltage limit: %f\n",minV);
  busIO.header(ioBuf);
  sprintf(ioBuf,"\nMaximum number of iterations: %d\n",max_iteration);
  busIO.header(ioBuf);
  sprintf(ioBuf,"\nConvergence tolerance: %f\n",tolerance);
  busIO.header(ioBuf);

  // set contingency
  p_factory->setContingency(contingency);

  // check contingency for isolated buses
  int t_lone = timer->createCategory("Check for Lone Bus");
  int lt_lone = ltimer.createCategory("Check for Lone Bus");
  timer->start(t_lone);
  ltimer.start(lt_lone);
  if (p_factory->checkLoneBus(busIO.getStream().get())) {
    sprintf(ioBuf,"\nIsolated bus found for contingency %s\n",
        contingency.p_name.c_str());
    busIO.header(ioBuf);
    // Exchange bus status between processors
    p_network->updateBuses();
  }
  timer->stop(t_lone);
  ltimer.stop(lt_lone);

  // set YBus components so that you can create Y matrix  
  int t_matv = timer->createCategory("Set Matrices and Vectors");
  int lt_matv = ltimer.createCategory("Set Matrices and Vectors");
  timer->start(t_matv);
  ltimer.start(lt_matv);
  p_factory->setYBus();

#if 0
  p_factory->setMode(gridpack::powerflow::YBus);
  gridpack::mapper::FullMatrixMap<CANetwork> ybusMap(p_network);
  boost::shared_ptr<gridpack::math::Matrix> orgYbus = ybusMap.mapToMatrix();
#endif

  //////////////////////////////////////////////////////////////
  p_factory->setMode(gridpack::powerflow::S_Cal);

  // make Sbus components to create S vector
  p_factory->setSBus();
  timer->stop(t_matv);
  ltimer.stop(lt_matv);

  // Set PQ
  int t_cmap = timer->createCategory("Create Mapper");
  int lt_cmap = ltimer.createCategory("Create Mapper");
  timer->start(t_cmap);
  ltimer.start(lt_cmap);
  p_factory->setMode(gridpack::powerflow::RHS);
  gridpack::mapper::BusVectorMap<CANetwork> vMap(p_network);
  timer->stop(t_cmap);
  ltimer.stop(lt_cmap);
  int t_mapm = timer->createCategory("Map Matrix");
  int lt_mapm = ltimer.createCategory("Map Matrix");
  int t_mapv = timer->createCategory("Map Vector");
  int lt_mapv = ltimer.createCategory("Map Vector");
  timer->start(t_mapv);
  ltimer.start(lt_mapv);
  boost::shared_ptr<gridpack::math::Vector> PQ = vMap.mapToVector();
  timer->stop(t_mapv);
  ltimer.stop(lt_mapv);
  timer->start(t_cmap);
  ltimer.start(lt_cmap);
  p_factory->setMode(gridpack::powerflow::Jacobian);
  gridpack::mapper::FullMatrixMap<CANetwork> jMap(p_network);
  timer->stop(t_cmap);
  ltimer.stop(lt_cmap);
  timer->start(t_mapm);
  ltimer.start(lt_mapm);
  boost::shared_ptr<gridpack::math::Matrix> J = jMap.mapToMatrix();
  timer->stop(t_mapm);
  ltimer.stop(lt_mapm);

  // Create X vector by cloning PQ
  boost::shared_ptr<gridpack::math::Vector> X(PQ->clone());

  // Create linear solver
  gridpack::math::LinearSolver solver(*J);
  int t_solv = timer->createCategory("Solve Powerflow Equations");
  int lt_solv = ltimer.createCategory("Solve Powerflow Equations");
  timer->start(t_solv);
  ltimer.start(lt_solv);
  solver.configure(cursor);

  tol = 2.0*tolerance;
  int iter = 0;

  // First iteration
  X->zero(); //might not need to do this
  try {
    solver.solve(*PQ, *X);
  } catch (const gridpack::Exception e) {
    busIO.header("Solver failure\n\n");
    p_factory->clearContingency(contingency);
    timer->stop(t_solv);
    ltimer.stop(lt_solv);
    timer->stop(t_task);
    ltimer.stop(lt_task);
    ltimer.dump(busIO.getStream());
    return;
  }
  tol = PQ->normInfinity();
  timer->stop(t_solv);
  ltimer.stop(lt_solv);

  int t_updt = timer->createCategory("Update Buses");
  int lt_updt = ltimer.createCategory("Update Buses");

  bool converged = false;
  if (real(tol) <= tolerance) converged = true;
  while (real(tol) > tolerance && iter < max_iteration) {
    // Push current values in X vector back into network components
    timer->start(t_mapv);
    ltimer.start(lt_mapv);
    p_factory->setMode(gridpack::powerflow::RHS);
    vMap.mapToBus(X);
    timer->stop(t_mapv);
    ltimer.stop(lt_mapv);

    // Exchange data between ghost buses (I don't think we need to exchange data
    // between branches)
    timer->start(t_updt);
    ltimer.start(lt_updt);
    p_network->updateBuses();
    timer->stop(t_updt);
    ltimer.stop(lt_updt);

    // Create new versions of Jacobian and PQ vector
    timer->start(t_mapv);
    ltimer.start(lt_mapv);
    vMap.mapToVector(PQ);
    timer->stop(t_mapv);
    ltimer.stop(lt_mapv);
    timer->start(t_mapm);
    ltimer.start(lt_mapm);
    p_factory->setMode(gridpack::powerflow::Jacobian);
    jMap.mapToMatrix(J);
    timer->stop(t_mapm);
    ltimer.stop(lt_mapm);

    // Create linear solver
    timer->start(t_solv);
    ltimer.start(lt_solv);
    X->zero(); //might not need to do this
    try {
      solver.solve(*PQ, *X);
    } catch (const gridpack::Exception e) {
      busIO.header("Solver failure\n\n");
      p_factory->clearContingency(contingency);
      timer->stop(t_solv);
      ltimer.stop(lt_solv);
      timer->stop(t_task);
      ltimer.stop(lt_task);
      ltimer.dump(busIO.getStream());
      return;
    }
    timer->stop(t_solv);
    ltimer.stop(lt_solv);

    tol = PQ->normInfinity();
    if (real(tol) <= tolerance) converged = true;
    sprintf(ioBuf,"\nIteration %d Tol: %12.6e\n",iter+1,real(tol));
    busIO.header(ioBuf);
    iter++;
  }
  // Push final result back onto buses
  timer->start(t_mapv);
  ltimer.start(lt_mapv);
  p_factory->setMode(gridpack::powerflow::RHS);
  vMap.mapToBus(X);
  timer->stop(t_mapv);
  ltimer.stop(lt_mapv);

  int t_out = timer->createCategory("Write Output to File");
  int lt_out = ltimer.createCategory("Write Output to File");
  timer->start(t_out);
  ltimer.start(lt_out);
  if (converged) {
    // Check for any violations
    bool bus_ok, branch_ok;
    p_factory->checkContingencies(minV,maxV,&bus_ok,&branch_ok);
    if (!bus_ok || !branch_ok) {
      sprintf(ioBuf,"\n   Violation found for contingency %s\n",
         contingency.p_name.c_str());
    } else {
      sprintf(ioBuf,"\n   No violations found for contingency %s\n",
         contingency.p_name.c_str());
    }
    busIO.header(ioBuf);

    if (!branch_ok) {
      branchIO.header("\n   Branch Violations\n");
      branchIO.header("\n        Bus 1       Bus 2     Tag          P"
	  "                    Q          Rating   Overloading\n");
      branchIO.write("flow");
    }

    if (!bus_ok) {
      busIO.header("\n   Bus Voltages Violations\n");
      busIO.header("\n   Bus Number      Voltage Magnitude\n");
      busIO.write("violations_only");
    }
  }

  if (converged) {
    sprintf(ioBuf,"\nContingency evaluation converged\n");
  } else {
    sprintf(ioBuf,"\nContingency evaluation did NOT converge\n");
  }
  busIO.header(ioBuf);

  timer->stop(t_out);
  ltimer.stop(lt_out);
  // clear contingency
  p_factory->clearContingency(contingency);
  p_factory->clearLoneBus();
  ltimer.stop(lt_task);
  ltimer.dump(busIO.getStream());
  busIO.close();
  timer->stop(t_task);
}

