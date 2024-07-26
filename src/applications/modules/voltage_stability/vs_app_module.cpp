/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   vs_app.cpp
 * @author Bruce Palmer
 * @date   2018-06-20 11:07:20 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "vs_app_module.hpp"
#include "vs_factory_module.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/parser/PTI33_parser.hpp"
#include "gridpack/parser/PTI34_parser.hpp"
#include "gridpack/parser/PTI35_parser.hpp"
#include "gridpack/parser/MAT_parser.hpp"
#include "gridpack/export/PSSE34Export.hpp"
#include "gridpack/export/PSSE33Export.hpp"
#include "gridpack/export/PSSE23Export.hpp"
#include "gridpack/parser/GOSS_parser.hpp"
#include "gridpack/math/math.hpp"
#include "vs_helper.hpp"
#include "gridpack/utilities/string_utils.hpp"
#include <iostream>
#include <fstream>

#define USE_REAL_VALUES

/**
 * Basic constructor
 */
gridpack::voltage_stability::VSAppModule::VSAppModule(void)
{
  p_no_print = false;
  current_increment = 0.0;
  p_bPVAnlyDone = true;
  zone = 0;
  sink_area = 0;
  src_area = 0;
  gt = 0.0;
  lt = 0.0;
  PV_header = false;
}

/**
 * Basic destructor
 */
gridpack::voltage_stability::VSAppModule::~VSAppModule(void)
{
}

enum Parser{PTI23, PTI33, PTI34, PTI35, MAT_POWER, GOSS};

/**
 * Read in and partition the powerflow network. The input file is read
 * directly from the Powerflow block in the configuration file so no
 * external file names or parameters need to be passed to this routine
 * @param network pointer to a VSNetwork object. This should not have any
 * buses or branches defined on it.
 * @param config point to open configuration file
 * @param idx index of configuration to use if set to a non-negative value
 */
void gridpack::voltage_stability::VSAppModule::readNetwork(
    boost::shared_ptr<VSNetwork> &network,
    gridpack::utility::Configuration *config, int idx)
{
  p_network = network;
  p_comm = network->communicator();
  p_config = config;

  // load input file
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);

  // read configuration file
  config->enableLogging(&std::cout);

  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Powerflow");
  if (cursor == NULL) {
    printf("No Powerflow block detected in input deck\n");
  }
  std::string filename;
  int filetype = PTI23;
  if (idx == -1) {
    if (!cursor->get("networkConfiguration",&filename)) {
      if (cursor->get("networkConfiguration_v33",&filename)) {
        filetype = PTI33;
      } else if (cursor->get("networkConfiguration_v34",&filename)) {
        filetype = PTI34;
      } else if (cursor->get("networkConfiguration_v35",&filename)) {
        filetype = PTI35;
      } else if (cursor->get("networkConfiguration_mat",&filename)) {
        filetype = MAT_POWER;
      } else if (cursor->get("networkConfiguration_GOSS",&filename)) {
        filetype = GOSS;
      } else {
        printf("No network configuration file specified\n");
        return;
      }
    }
  } else if (idx >= 0) {
    gridpack::utility::Configuration::CursorPtr network_cursor;
    network_cursor = config->getCursor(
        "Configuration.Powerflow.networkFiles");
    gridpack::utility::Configuration::ChildCursors files;
    if (network_cursor) network_cursor->children(files);
    if (idx < files.size()) {
      if (!files[idx]->get("networkConfiguration",&filename)) {
        if (files[idx]->get("networkConfiguration_v33",&filename)) {
          filetype = PTI33;
        } else if (files[idx]->get("networkConfiguration_v34",&filename)) {
          filetype = PTI34;
        } else if (files[idx]->get("networkConfiguration_v35",&filename)) {
          filetype = PTI35;
        } else if (cursor->get("networkConfiguration_mat",&filename)) {
          filetype = MAT_POWER;
        } else {
          printf("Unknown network configuration file specified\n");
          return;
        }
      }
    } else {
      printf("Unknown file index\n");
      return;
    }
  } else {
    printf("No network configuration file specified\n");
    return;
  }
  // Convergence and iteration parameters
  p_tolerance = cursor->get("tolerance",1.0e-6);
  p_qlim = cursor->get("qlim",0);
  p_max_iteration = cursor->get("maxIteration",50);
  ComplexType tol;
  // Phase shift sign
  max_increment = cursor->get("maxIncrement",0.0);
  increment = cursor->get("TransferIncrement",0.0);
  src_area = cursor->get("SourceArea",0);
  sink_area = cursor->get("SinkArea",0);
  double phaseShiftSign = cursor->get("phaseShiftSign",1.0);

  int t_pti = timer->createCategory("Powerflow: Network Parser");
  timer->start(t_pti);
  if (filetype == PTI23) {
    gridpack::parser::PTI23_parser<VSNetwork> parser(network);
#ifdef USE_GOSS
    char sbuf[256], sbuf2[256];
    sprintf(sbuf,"{ \"simulation_id\": \"%s\"}",simID.c_str());
    sprintf(sbuf2,"reply.%s.%s\n",filename.c_str(),p_simID.c_str());
    p_goss_client.publish(networkFile,sbuf,sbuf2);
    std::vector<std::string> fileVec = p_goss_client.subscribeFileAsVector(std::string(sbuf2));
    parser.parse(fileVec);
#else
    try {
      parser.parse(filename.c_str());
    } catch (const gridpack::Exception e) {
      std::string w(e.what());
      if (!p_no_print) {
        char ebuf[512];
        sprintf(ebuf,"p[%d] unable to open network file: %s with error: %s\n",
            filename.c_str(),w.c_str());
        if (p_comm.rank() == 0) {
          printf("%s",ebuf);
        }
      }
      timer->stop(t_total);
    }
#endif
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  } else if (filetype == PTI33) {
    gridpack::parser::PTI33_parser<VSNetwork> parser(network);
#ifdef USE_GOSS
    char sbuf[256], sbuf2[256];
    sprintf(sbuf,"{ \"simulation_id\": \"%s\"}",simID.c_str());
    sprintf(sbuf2,"reply.%s.%s\n",filename.c_str(),p_simID.c_str());
    p_goss_client.publish(networkFile,sbuf,sbuf2);
    std::vector<std::string> fileVec = p_goss_client.subscribeFileAsVector(std::string(sbuf2));
    parser.parse(fileVec);
#else
    parser.parse(filename.c_str());
#endif
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  } else if (filetype == PTI34) {
    gridpack::parser::PTI34_parser<VSNetwork> parser(network);
#ifdef USE_GOSS
    char sbuf[256], sbuf2[256];
    sprintf(sbuf,"{ \"simulation_id\": \"%s\"}",simID.c_str());
    sprintf(sbuf2,"reply.%s.%s\n",filename.c_str(),p_simID.c_str());
    p_goss_client.publish(networkFile,sbuf,sbuf2);
    std::vector<std::string> fileVec = p_goss_client.subscribeFileAsVector(std::string(sbuf2));
    parser.parse(fileVec);
#else
    parser.parse(filename.c_str());
#endif
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  } else if (filetype == PTI35) {
    gridpack::parser::PTI35_parser<VSNetwork> parser(network);
#ifdef USE_GOSS
    char sbuf[256], sbuf2[256];
    sprintf(sbuf,"{ \"simulation_id\": \"%s\"}",simID.c_str());
    sprintf(sbuf2,"reply.%s.%s\n",filename.c_str(),p_simID.c_str());
    p_goss_client.publish(networkFile,sbuf,sbuf2);
    std::vector<std::string> fileVec = p_goss_client.subscribeFileAsVector(std::string(sbuf2));
    parser.parse(fileVec);
#else
    parser.parse(filename.c_str());
#endif
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  } else if (filetype == MAT_POWER) {
    gridpack::parser::MAT_parser<VSNetwork> parser(network);
#ifdef USE_GOSS
    char sbuf[256], sbuf2[256];
    sprintf(sbuf,"{ \"simulation_id\": \"%s\"}",simID.c_str());
    sprintf(sbuf2,"reply.%s.%s\n",filename.c_str(),p_simID.c_str());
    p_goss_client.publish(networkFile,sbuf,sbuf2);
    std::vector<std::string> fileVec = p_goss_client.subscribeFileAsVector(std::string(sbuf2));
    parser.parse(fileVec);
#else
    parser.parse(filename.c_str());
#endif
  } else if (filetype == GOSS) {
    gridpack::parser::GOSS_parser<VSNetwork> parser(network);
    parser.parse(filename.c_str());
  }
  timer->stop(t_pti);

  // Create serial IO object to export data from buses
  p_busIO.reset(new gridpack::serial_io::SerialBusIO<VSNetwork>(512,network));

  // Create serial IO object to export data from branches
  p_branchIO.reset(new gridpack::serial_io::SerialBranchIO<VSNetwork>(512,network));
  char ioBuf[128];

  if (!p_no_print) {
    sprintf(ioBuf,"\nMaximum number of iterations: %d\n",p_max_iteration);
    p_busIO->header(ioBuf);
    sprintf(ioBuf,"\nConvergence tolerance: %f\n",p_tolerance);
    p_busIO->header(ioBuf);
  }

  // partition network
  int t_part = timer->createCategory("Powerflow: Partition");
  timer->start(t_part);
  network->partition();
  timer->stop(t_part);
  timer->stop(t_total);
}

/**
 * Set up exchange buffers and other internal parameters and initialize
 * network components using data from data collection
 */
void gridpack::voltage_stability::VSAppModule::initialize()
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);

  // create factory
  p_factory.reset(new gridpack::voltage_stability::VSFactoryModule(p_network));
  int t_load = timer->createCategory("Powerflow: Factory Load");
  timer->start(t_load);
  p_factory->load();
  // p_factory->dumpData();
  timer->stop(t_load);

  // set network components using factory
  int t_setc = timer->createCategory("Powerflow: Factory Set Components");
  timer->start(t_setc);
  p_factory->setComponents();
  timer->stop(t_setc);

  // Set up bus data exchange buffers. Need to decide what data needs to be
  // exchanged
  int t_setx = timer->createCategory("Powerflow: Factory Set Exchange");
  timer->start(t_setx);
  p_factory->setExchange();
  timer->stop(t_setx);

  // Create bus data exchange
  int t_updt = timer->createCategory("Powerflow: Bus Update");
  timer->start(t_updt);
  p_network->initBusUpdate();
  timer->stop(t_updt);
  timer->stop(t_total);
}

/**
 * Reinitialize calculation from data collections
 */
void gridpack::voltage_stability::VSAppModule::reload()
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_load = timer->createCategory("Powerflow: Factory Load");
  timer->start(t_load);
  p_factory->load();
  timer->stop(t_load);
}

/**
 * Execute the iterative solve portion of the application using a
 * hand-coded Newton-Raphson solver
 * @return false if an error was encountered in the solution
 */
bool gridpack::voltage_stability::VSAppModule::solve()
{
  bool ret = true;
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);
  p_factory->clearViolations();
  gridpack::ComplexType tol = 2.0*p_tolerance;
  int iter = 0;
  bool repeat = true;
  int int_repeat = 0;
  while (repeat) {
    iter = 0;
    tol = 2.0*p_tolerance;
    int_repeat ++;
    char ioBuf[128];
    if (!p_no_print) {
      sprintf (ioBuf," repeat time = %d \n", int_repeat);
      p_busIO->header(ioBuf);
    }

    // set YBus components so that you can create Y matrix
    int t_fact = timer->createCategory("Powerflow: Factory Operations");
    timer->start(t_fact);
    p_factory->setYBus();
    timer->stop(t_fact);

    int t_cmap = timer->createCategory("Powerflow: Create Mappers");
    timer->start(t_cmap);
    p_factory->setMode(VS_YBus); 

#if 0
    gridpack::mapper::FullMatrixMap<VSNetwork> mMap(p_network);
#endif
    timer->stop(t_cmap);
    int t_mmap = timer->createCategory("Powerflow: Map to Matrix");
    timer->start(t_mmap);
#if 0
    gridpack::mapper::FullMatrixMap<VSNetwork> mMap(p_network);
    boost::shared_ptr<gridpack::math::Matrix> Y = mMap.mapToMatrix();
    p_busIO->header("\nY-matrix values\n");
    //  Y->print();
    Y->save("Ybus.m");
#endif
    timer->stop(t_mmap);

    // make Sbus components to create S vector
    timer->start(t_fact);
    p_factory->setSBus();
    timer->stop(t_fact);
    //  p_busIO->header("\nIteration 0\n");

    // Set PQ
    timer->start(t_cmap);
    p_factory->setMode(VS_RHS); 
    gridpack::mapper::BusVectorMap<VSNetwork> vMap(p_network);
    timer->stop(t_cmap);
    int t_vmap = timer->createCategory("Powerflow: Map to Vector");
    timer->start(t_vmap);

#ifdef USE_REAL_VALUES
    boost::shared_ptr<gridpack::math::RealVector> PQ = vMap.mapToRealVector();
#else
    boost::shared_ptr<gridpack::math::Vector> PQ = vMap.mapToVector();
#endif
    timer->stop(t_vmap);
    gridpack::ComplexType  tol_org = PQ->normInfinity();
    if (!p_no_print) {
      sprintf (ioBuf,"\n----------test Iteration 0, before VS solve, Tol: %12.6e \n", real(tol_org));
      p_busIO->header(ioBuf);
    }
    //  PQ->print();
    timer->start(t_cmap);
    p_factory->setMode(VS_Jacobian);
    gridpack::mapper::FullMatrixMap<VSNetwork> jMap(p_network);
    timer->stop(t_cmap);
    timer->start(t_mmap);

#ifdef USE_REAL_VALUES
    boost::shared_ptr<gridpack::math::RealMatrix> J = jMap.mapToRealMatrix();
#else
    boost::shared_ptr<gridpack::math::Matrix> J = jMap.mapToMatrix();
#endif
    timer->stop(t_mmap);
    //  p_busIO->header("\nJacobian values\n");
    //  J->print();

    // Create X vector by cloning PQ
#ifdef USE_REAL_VALUES
    boost::shared_ptr<gridpack::math::RealVector> X(PQ->clone());
#else
    boost::shared_ptr<gridpack::math::Vector> X(PQ->clone());
#endif

    gridpack::utility::Configuration::CursorPtr cursor;
    cursor = p_config->getCursor("Configuration.Powerflow");
    // Create linear solver
    int t_csolv = timer->createCategory("Powerflow: Create Linear Solver");
    timer->start(t_csolv);
#ifdef USE_REAL_VALUES
    gridpack::math::RealLinearSolver solver(*J);
#else
    gridpack::math::LinearSolver solver(*J);
#endif
    solver.configure(cursor);
    timer->stop(t_csolv);

    // First iteration
    X->zero(); //might not need to do this
    //p_busIO->header("\nCalling solver\n");
    int t_lsolv = timer->createCategory("Powerflow: Solve Linear Equation");
    timer->start(t_lsolv);
    //    char dbgfile[32];
    //    sprintf(dbgfile,"j0.bin");
    //    J->saveBinary(dbgfile);
    //    sprintf(dbgfile,"pq0.bin");
    //    PQ->saveBinary(dbgfile);
    try {
      solver.solve(*PQ, *X);
    } catch (const gridpack::Exception e) {
      std::string w(e.what());
      if (!p_no_print) {
        sprintf(ioBuf,"p[%d] hit exception: %s\n",
            p_network->communicator().rank(),
            w.c_str());
        p_busIO->header(ioBuf);
        p_busIO->header("Solver failure\n\n");
      }
      timer->stop(t_lsolv);
      timer->stop(t_total);

      return false;
    }
    timer->stop(t_lsolv);
    tol = PQ->normInfinity();
    // Create timer for map to bus
    int t_bmap = timer->createCategory("Powerflow: Map to Bus");
    int t_updt = timer->createCategory("Powerflow: Bus Update");
    if (!p_no_print) {
      sprintf(ioBuf,"\nIteration %d Tol: %12.6e\n",iter+1,real(tol));
      p_busIO->header(ioBuf);
    }										

    while (real(tol) > p_tolerance && iter < p_max_iteration) {
      // Push current values in X vector back into network components
      // Need to implement setValues method in VSBus class in order for this to
      // work
      timer->start(t_bmap);
      p_factory->setMode(VS_RHS);
      vMap.mapToBus(X);
      timer->stop(t_bmap);

      // Exchange data between ghost buses (I don't think we need to exchange data
      // between branches)
      timer->start(t_updt);
      //  p_factory->checkQlimViolations();
      p_network->updateBuses();
      timer->stop(t_updt);

      // Create new versions of Jacobian and PQ vector
      timer->start(t_vmap);
#ifdef USE_REAL_VALUES
      vMap.mapToRealVector(PQ);
#else
      vMap.mapToVector(PQ);
#endif
      //   p_busIO->header("\nnew PQ vector at iter %d\n",iter);
      //   PQ->print();
      timer->stop(t_vmap);
      timer->start(t_mmap);
      p_factory->setMode(VS_Jacobian);
#ifdef USE_REAL_VALUES
      jMap.mapToRealMatrix(J);
#else
      jMap.mapToMatrix(J);
#endif
      timer->stop(t_mmap);

      // Create linear solver
      timer->start(t_lsolv);
      X->zero(); //might not need to do this
      //    sprintf(dbgfile,"j%d.bin",iter+1);
      //    J->saveBinary(dbgfile);
      //    sprintf(dbgfile,"pq%d.bin",iter+1);
      //    PQ->saveBinary(dbgfile);
      try {
        solver.solve(*PQ, *X);
      } catch (const gridpack::Exception e) {
        std::string w(e.what());
        if (!p_no_print) {
          sprintf(ioBuf,"p[%d] hit exception: %s\n",
              p_network->communicator().rank(),
              w.c_str());
          p_busIO->header(ioBuf);
          p_busIO->header("Solver failure\n\n");
        }
        timer->stop(t_lsolv);
        timer->stop(t_total);
        return false;
      }
      timer->stop(t_lsolv);

      tol = PQ->normInfinity();
      if (!p_no_print) {
        sprintf(ioBuf,"\nIteration %d Tol: %12.6e\n",iter+1,real(tol));
        p_busIO->header(ioBuf);
      }
      iter++;
      if (real(tol)> 100.0*real(tol_org)){
        ret = false;
        if (!p_no_print) {
          sprintf (ioBuf,"\n-------------current iteration tol bigger than 100 times of original tol, power flow diverge\n");
          p_busIO->header(ioBuf);
        }
        break;
      }
    }

    if (iter >= p_max_iteration) ret = false;
    if (p_qlim == 0) {
      repeat = false;
    } else {
      if (p_factory->checkQlimViolations()) {
        repeat =false;
      } else {
        if (!p_no_print) {
          sprintf (ioBuf,"There are Qlim violations at iter =%d\n", iter);
          p_busIO->header(ioBuf);
        }
      }
    }
    // Push final result back onto buses
    timer->start(t_bmap);
    p_factory->setMode(VS_RHS);
    vMap.mapToBus(X);
    timer->stop(t_bmap);

    // Make sure that ghost buses have up-to-date values before printing out
    // results
    timer->start(t_updt);
    p_network->updateBuses();
    timer->stop(t_updt);
  }
  timer->stop(t_total);
  return ret;

}
/**
 * Execute the iterative solve portion of the application using a library
 * non-linear solver
 * @return false if an error was caught in the solution algorithm
 */
bool gridpack::voltage_stability::VSAppModule::nl_solve()
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);
  p_factory->clearViolations();

  int t_fact = timer->createCategory("Powerflow: Factory Operations");
  timer->start(t_fact);
  p_factory->setYBus();
  p_factory->setSBus();
  timer->stop(t_fact);

  // Solve the problem
  bool useNewton(false);
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Powerflow");
  useNewton = cursor->get("UseNewton", useNewton);


  int t_lsolv = timer->createCategory("Powerflow: Non-Linear Solver");
  timer->start(t_lsolv);

  VSSolverHelper helper(p_factory, p_network);
  math::RealNonlinearSolver::JacobianBuilder jbuildf = boost::ref(helper);
  math::RealNonlinearSolver::FunctionBuilder fbuildf = boost::ref(helper);

  boost::scoped_ptr<math::RealNonlinearSolver> solver;
  if (useNewton) {
    math::RealNewtonRaphsonSolver *tmpsolver =
      new math::RealNewtonRaphsonSolver(*(helper.J), jbuildf, fbuildf);
    tmpsolver->tolerance(p_tolerance);
    tmpsolver->maximumIterations(p_max_iteration);
    solver.reset(tmpsolver);
  } else {
    solver.reset(new math::RealNonlinearSolver(*(helper.J),
          jbuildf, fbuildf));
  }

  bool ret = true;
  try {
    solver->configure(cursor);
    solver->solve(*helper.X);
    helper.update(*helper.X);
  } catch (const Exception& e) {
    std::cerr << e.what() << std::endl;
    timer->stop(t_lsolv);
    timer->stop(t_total);
    ret = false;
  }

  timer->stop(t_lsolv);
  timer->stop(t_total);
  return ret;
}


/**
 * Write out results of powerflow calculation to standard output or a file
 */
void gridpack::voltage_stability::VSAppModule::write()
{
  if (p_no_print) return;
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);
  int t_write = timer->createCategory("Powerflow: Write Results");
  timer->start(t_write);
  p_branchIO->header("\n   Branch Power Flow\n");
  p_branchIO->header("\n        Bus 1       Bus 2   CKT         P"
                  "                    Q                   FVSI\n");
  p_branchIO->write();
  //p_branchIO->write("record");


  p_busIO->header("\n   Generator Power\n");
  p_busIO->header("\n   Bus Number  GenID        Pgen              Qgen\n");
  p_busIO->write("power");
  p_busIO->header("\n   Bus Voltages and Phase Angles\n");
  p_busIO->header("\n   Bus Number      Phase Angle      Voltage Magnitude\n");
  p_busIO->write();
  //p_busIO->write("record");
  timer->stop(t_write);
  timer->stop(t_total);
}

void gridpack::voltage_stability::VSAppModule::writeBus(const char *signal)
{
  if (p_no_print) return;
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);
  int t_write = timer->createCategory("Powerflow: Write Results");
  timer->start(t_write);
  timer->start(t_total);
  p_busIO->write(signal);
  timer->stop(t_write);
  timer->stop(t_total);
}

void gridpack::voltage_stability::VSAppModule::writeCABus()
{
  if (p_no_print) return;
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Contingency: Total Application");
  timer->start(t_total);
  int t_write = timer->createCategory("Contingency: Write Results");
  timer->start(t_write);
  p_busIO->write("ca");
//  p_busIO->write(signal);
  timer->stop(t_write);
  timer->stop(t_total);
}

void gridpack::voltage_stability::VSAppModule::writeBranch(const char *signal)
{
  if (p_no_print) return;
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);
  int t_write = timer->createCategory("Powerflow: Write Results");
  timer->start(t_write);
  timer->start(t_total);
  p_branchIO->write(signal);
  timer->stop(t_write);
  timer->stop(t_total);
}

void gridpack::voltage_stability::VSAppModule::writeCABranch()
{
  if (p_no_print) return;
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Contingency: Total Application");
  timer->start(t_total);
  int t_write = timer->createCategory("Contingency: Write Results");
  timer->start(t_write);
  p_branchIO->write("flow");
  timer->stop(t_write);
  timer->stop(t_total);
}

std::vector<std::string> gridpack::voltage_stability::VSAppModule::writeBusString(
    const char *signal)
{
  std::vector<std::string> ret;
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Contingency: Total Application");
  timer->start(t_total);
  int t_write = timer->createCategory("Contingency: Write Results");
  timer->start(t_write);
  ret = p_busIO->writeStrings(signal);
  timer->stop(t_write);
  timer->stop(t_total);
  return ret;
}

std::vector<std::string> gridpack::voltage_stability::VSAppModule::writeBranchString(
    const char *signal)
{
  std::vector<std::string> ret;
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Contingency: Total Application");
  timer->start(t_total);
  int t_write = timer->createCategory("Contingency: Write Results");
  timer->start(t_write);
  ret = p_branchIO->writeStrings(signal);
  timer->stop(t_write);
  timer->stop(t_total);
  return ret;
}

void gridpack::voltage_stability::VSAppModule::writeHeader(const char *msg)
{
  if (p_no_print) return;
  p_busIO->header(msg);
}

/**
 * Redirect output from standard out
 * @param filename name of file to write results to
 */
void gridpack::voltage_stability::VSAppModule::open(const char *filename)
{
  if (p_no_print) return;
  p_busIO->open(filename);
  p_branchIO->setStream(p_busIO->getStream());
}

void gridpack::voltage_stability::VSAppModule::close()
{
  if (p_no_print) return;
  p_busIO->close();
  p_branchIO->setStream(p_busIO->getStream());
}

/**
 * Print string. This can be used to direct output to the file opened using
 * the open command
 * @param buf string to be printed
 */
void gridpack::voltage_stability::VSAppModule::print(const char *buf)
{
  if (p_no_print) return;
  p_busIO->header(buf);
}

/**
 * Export final configuration to PSS/E v34 formatted file
 * @param filename name of file to store network configuration
 */
void gridpack::voltage_stability::VSAppModule::exportPSSE34(std::string &filename)
{
  if (p_no_print) return;
  gridpack::expnet::PSSE34Export<VSNetwork> exprt(p_network);
  exprt.writeFile(filename);
}

/**
 * Export final configuration to PSS/E v33 formatted file
 * @param filename name of file to store network configuration
 */
void gridpack::voltage_stability::VSAppModule::exportPSSE33(std::string &filename)
{
  if (p_no_print) return;
  gridpack::expnet::PSSE33Export<VSNetwork> exprt(p_network);
  exprt.writeFile(filename);
}

/**
 * Export final configuration to PSS/E v23 formatted file
 * @param filename name of file to store network configuration
 */
void gridpack::voltage_stability::VSAppModule::exportPSSE23(std::string &filename)
{
  //if (p_no_print) return;
  gridpack::expnet::PSSE23Export<VSNetwork> exprt(p_network);
  exprt.writeFile(filename);
}

/**
 * Save results of powerflow calculation to data collection objects
 */
void gridpack::voltage_stability::VSAppModule::saveData()
{
  p_factory->saveData();
}

/**
 * Save results of voltage_stability calculation to data collection objects
 * added by Renke, also modify the original bus mag, ang, 
 * and the original generator PG QG in the datacollection
 */
void gridpack::voltage_stability::VSAppModule::saveDataAlsotoOrg()
{
  p_factory->saveDataAlsotoOrg();
}

/**
 * get the power flow solution for the specific bus, vmag and v angle
 * @param bus original number, bus solution vmag and v angle
 * @return false if location of bus is not found in
 * network
 */

bool gridpack::voltage_stability::VSAppModule::getVSSolutionSingleBus(
    int bus_number, double &bus_mag, double &bus_angle)
{
	bool ret = true;
	std::vector<int> vec_busintidx;
	int ibus, nbus;
	gridpack::voltage_stability::VSBus *bus;
	
	vec_busintidx = p_network->getLocalBusIndices(bus_number);
	nbus = vec_busintidx.size();
	if (nbus == 0) ret = false;
	for(ibus=0; ibus<nbus; ibus++){
		bus = dynamic_cast<gridpack::voltage_stability::VSBus*>
		(p_network->getBus(vec_busintidx[ibus]).get());  //->getOriginalIndex()
		//printf("----renke debug VSAppModule::getVSSolutionSingleBus, \n");
		bus_mag=bus->getVoltage();
		double anglerads = bus->getPhase();
		double pi = 4.0*atan(1.0);
		bus_angle = 180.0*anglerads/pi;	
	}
	
	return ret;
}

/**
 * Set a contingency
 * @param event data describing location and type of contingency
 * @return false if location of contingency is not found in
 * network
 */
bool gridpack::voltage_stability::VSAppModule::setContingency(
    gridpack::voltage_stability::VSContingency &event)
{
  bool ret = true;
  if (event.p_type == VS_Generator) {
    int ngen = event.p_busid.size();
    int i, j, idx, jdx;
    for (i=0; i<ngen; i++) {
      idx = event.p_busid[i];
      std::string tag = event.p_genid[i];
      std::vector<int> lids = p_network->getLocalBusIndices(idx);
      if (lids.size() == 0) ret = false;
      gridpack::voltage_stability::VSBus *bus;
      for (j=0; j<lids.size(); j++) {
        jdx = lids[j];
        bus = dynamic_cast<gridpack::voltage_stability::VSBus*>(
            p_network->getBus(jdx).get());
        event.p_saveGenStatus[i] = bus->getGenStatus(tag);
        bus->setGenStatus(tag, false);
      }
    }
  } else if (event.p_type == VS_Branch) {
    int to, from;
    int nline = event.p_to.size();
    int i, j, idx, jdx;
    for (i=0; i<nline; i++) {
      to = event.p_to[i];
      from = event.p_from[i];
      std::string tag = event.p_ckt[i];
      std::vector<int> lids = p_network->getLocalBranchIndices(from,to);
      if (lids.size() == 0) ret = false;
      gridpack::voltage_stability::VSBranch *branch;
      for (j=0; j<lids.size(); j++) {
        jdx = lids[j];
        branch = dynamic_cast<gridpack::voltage_stability::VSBranch*>(
            p_network->getBranch(jdx).get());
        event.p_saveLineStatus[i] = branch->getBranchStatus(tag);
        branch->setBranchStatus(tag, false);
      }
    }
  } else {
    ret = false;
  }
  if (ret) {
    p_contingency_name = event.p_name;
  } else {
    p_contingency_name.clear();
  }
  p_factory->checkLoneBus();
  return ret;
}

/**
 * Return system to the state before the contingency
 * @param event data describing location and type of contingency
 * @return false if location of contingency is not found in network
 */
bool gridpack::voltage_stability::VSAppModule::unSetContingency(
    gridpack::voltage_stability::VSContingency &event)
{
  p_factory->clearLoneBus();
  bool ret = true;
  if (event.p_type == VS_Generator) {
    int ngen = event.p_busid.size();
    int i, j, idx, jdx;
    for (i=0; i<ngen; i++) {
      idx = event.p_busid[i];
      std::string tag = event.p_genid[i];
      std::vector<int> lids = p_network->getLocalBusIndices(idx);
      if (lids.size() == 0) ret = false;
      gridpack::voltage_stability::VSBus *bus;
      for (j=0; j<lids.size(); j++) {
        jdx = lids[j];
        bus = dynamic_cast<gridpack::voltage_stability::VSBus*>(
            p_network->getBus(jdx).get());
        bus->setGenStatus(tag, event.p_saveGenStatus[i]);
      }
    }
  } else if (event.p_type == VS_Branch) {
    int to, from;
    int nline = event.p_to.size();
    int i, j, idx, jdx;
    for (i=0; i<nline; i++) {
      to = event.p_to[i];
      from = event.p_from[i];
      std::string tag = event.p_ckt[i];
      std::vector<int> lids = p_network->getLocalBranchIndices(from,to);
      if (lids.size() == 0) ret = false;
      gridpack::voltage_stability::VSBranch *branch;
      for (j=0; j<lids.size(); j++) {
        jdx = lids[j];
        branch = dynamic_cast<gridpack::voltage_stability::VSBranch*>(
            p_network->getBranch(jdx).get());
        branch->setBranchStatus(tag,event.p_saveLineStatus[i]);
      }
    }
  } else {
    ret = false;
  }
  return ret;
}

/**
 * Set voltage limits on all buses
 * @param Vmin lower bound on voltages
 * @param Vmax upper bound on voltages
 */
void gridpack::voltage_stability::VSAppModule::setVoltageLimits(double Vmin, double Vmax)
{
  p_factory->setVoltageLimits(Vmin, Vmax);
}

/**
 * Check to see if there are any voltage violations in the network
 * @param area area number. If this parameter is included, only check for
 * violations in this area
 * @return true if no violations found
 */
bool gridpack::voltage_stability::VSAppModule::checkVoltageViolations()
{
  return p_factory->checkVoltageViolations();
}
bool gridpack::voltage_stability::VSAppModule::checkVoltageViolations(
 int area)
{
  return p_factory->checkVoltageViolations(area);
}

/**
 * Set "ignore" parameter on all buses with violations so that subsequent
 * checks are not counted as violations
 */
void gridpack::voltage_stability::VSAppModule::ignoreVoltageViolations()
{
  p_factory->ignoreVoltageViolations();
}

/**
 * Clear "ignore" parameter on all buses
 */
void gridpack::voltage_stability::VSAppModule::clearVoltageViolations()
{
  p_factory->clearVoltageViolations();
}

/**
 * Check to see if there are any line overload violations in the
 * network The last call checks for overloads on specific lines.
 * @param area area number. If this parameter is included, only check for
 * violations in this area
 * @param bus1 original index of "from" bus for branch
 * @param bus2 original index of "to" bus for branch
 * @param tags line IDs for individual lines
 * @param violations true if violation detected on branch, false otherwise
 * @return true if no violations found
 */
bool gridpack::voltage_stability::VSAppModule::checkLineOverloadViolations()
{
  return p_factory->checkLineOverloadViolations();
}
bool gridpack::voltage_stability::VSAppModule::checkLineOverloadViolations(int area)
{
  return p_factory->checkLineOverloadViolations(area);
}
bool gridpack::voltage_stability::VSAppModule::checkLineOverloadViolations(
    std::vector<int> &bus1, std::vector<int> &bus2,
    std::vector<std::string> &tags, std::vector<bool> &violations)
{
  return p_factory->checkLineOverloadViolations(bus1,bus2,tags,violations);
}

/**
 * Check to see if there are any Q limit violations in the network
 * @param area only check for violations in specified area
 * @return true if no violations found
 */
bool gridpack::voltage_stability::VSAppModule::checkQlimViolations()
{
  return p_factory->checkQlimViolations();
}
bool gridpack::voltage_stability::VSAppModule::checkQlimViolations(int area)
{
  return p_factory->checkQlimViolations(area);
}

/**
 * Clear changes that were made for Q limit violations and reset
 * system to its original state
 */
void gridpack::voltage_stability::VSAppModule::clearQlimViolations()
{
  p_factory->clearQlimViolations();
}

/**
 * Reset voltages to values in network configuration file
 */
void gridpack::voltage_stability::VSAppModule::resetVoltages()
{
  p_factory->resetVoltages();
}

/**
 * Scale generator real power. If zone less than 1 then scale all
 * generators in the area.
 * @param scale factor to scale real power generation
 * @param area index of area for scaling generation
 * @param zone index of zone for scaling generation
 */
void gridpack::voltage_stability::VSAppModule::scaleGeneratorRealPower(
    double scale, int area, int zone)
{
  p_factory->scaleGeneratorRealPower(scale,area,zone);
}

/**
 * Scale load power. If zone less than 1 then scale all
 * loads in the area.
 * @param scale factor to scale load real power
 * @param area index of area for scaling load
 * @param zone index of zone for scaling load
 */
void gridpack::voltage_stability::VSAppModule::scaleLoadPower(
    double scale, int area, int zone)
{
  p_factory->scaleLoadPower(scale,area,zone);
}

/**
 * Increment generators real power based off specified value. 
 * Increment generators in specified area.
 * @param transfer value to increment generators real power
 * @param area index of area for incrementing generation
 * @param zone index of zone for incrementing generation
 * @param total power generation of an area
 */
void gridpack::voltage_stability::VSAppModule::IncrementGeneratorRealPower(
    double inc, int area, int zone, double gt)
{
  p_factory->IncrementGeneratorRealPower(inc,area,zone,gt);
}

/**
 * Increment load power based off specified value. 
 * Increment loads in specified area.
 * @param transfer value to increment load real power
 * @param area index of area for incrementing load
 * @param zone index of zone for incrementing load
 * @param total active power demand of the area
 */
void gridpack::voltage_stability::VSAppModule::IncrementLoadPower(
    double inc, int area, int zone, double lt)
{
  p_factory->IncrementLoadPower(inc,area,zone,lt);
}

/**
 * Return the total real power load for all loads in the zone. If zone
 * less than 1, then return the total load real power for the area
 * @param area index of area
 * @param zone index of zone
 * @return total load
 */
double gridpack::voltage_stability::VSAppModule::getTotalLoadRealPower(int area, int zone)
{
  return p_factory->getTotalLoadRealPower(area,zone);
}

/**
 * Return the current real power generation and the maximum and minimum total
 * power generation for all generators in the zone. If zone is less than 1
 * then return values for all generators in the area
 * @param area index of area
 * @param zone index of zone
 * @param total total real power generation
 * @param pmin minimum allowable real power generation
 * @param pmax maximum available real power generation
 */
void gridpack::voltage_stability::VSAppModule::getGeneratorMargins(int area,
    int zone, double *total, double *pmin, double *pmax)
{
  p_factory->getGeneratorMargins(area,zone,total,pmin,pmax);
}

/**
 * Reset power of loads and generators to original values
 */
void gridpack::voltage_stability::VSAppModule::resetPower()
{
  p_factory->resetPower();
}

/**
 * Write real time path rating diagnostics
 * @param src_area generation area
 * @param src_zone generation zone
 * @param load_area load area
 * @param load_zone load zone
 * @param gen_scale scale factor for generation
 * @param load_scale scale factor for loads
 * @param file name of file containing diagnostics
 */
void gridpack::voltage_stability::VSAppModule::writeRTPRDiagnostics(
    int src_area, int src_zone, int load_area,
    int load_zone, double gen_scale, double load_scale, const char *file)
{
  if (p_no_print) return;
  p_factory->setRTPRParams(src_area,src_zone,load_area,load_zone,
      gen_scale,load_scale);
  p_busIO->open(file);
  double gtotal, ltotal, pmin, pmax, scaled;
  p_factory->getGeneratorMargins(src_area, src_zone,&gtotal,&pmin,&pmax);
  ltotal = p_factory->getTotalLoadRealPower(load_area,load_zone);
  if (gen_scale > 0.0) {
    scaled = gtotal + gen_scale*(pmax-gtotal);
  } else {
    scaled = gtotal + gen_scale*(gtotal-pmin);
  }
  char sbuf[128];
  sprintf(sbuf,"Total Generation:         %16.4f\n",gtotal);
  p_busIO->header(sbuf);
  sprintf(sbuf,"  Minimum Generation:     %16.4f\n",pmin);
  p_busIO->header(sbuf);
  sprintf(sbuf,"  Maximum Generation:     %16.4f\n",pmax);
  p_busIO->header(sbuf);
  sprintf(sbuf,"  Generator Scale Factor: %16.4f\n",gen_scale);
  p_busIO->header(sbuf);
  sprintf(sbuf,"  Scaled Generation:      %16.4f\n",scaled);
  p_busIO->header(sbuf);
  p_busIO->header("\nIndividual Scaled Generators\n");
  sprintf(sbuf,"\n     Bus ID   Status Area Zone     Real Power   Scaled Power"
      "           Pmin           Pmax\n\n");
  p_busIO->header(sbuf);
  p_busIO->write("src_gen");
  sprintf(sbuf,"\nTotal Load:               %16.4f\n",ltotal);
  p_busIO->header(sbuf);
  sprintf(sbuf,"  Load Scale Factor:      %16.4f\n",load_scale);
  p_busIO->header(sbuf);
  sprintf(sbuf,"  Scaled Load:            %16.4f\n",load_scale*ltotal);
  p_busIO->header(sbuf);
  p_busIO->header("\nIndividual Scaled Loads\n");
  sprintf(sbuf,"\n     Bus ID   Status Area Zone     Real Power   Scaled Power"
      " Reactive Power   Scaled Power\n\n");
  p_busIO->header(sbuf);
  p_busIO->write("sink_load");
  p_busIO->close();
}

/**
 * Get strings documenting contingency failures. Strings are *not* terminated
 * with a carriage return
 */
std::vector<std::string> gridpack::voltage_stability::VSAppModule::getContingencyFailures()
{
  std::vector<std::string> ret;
  std::vector<gridpack::voltage_stability::VSFactoryModule::Violation> violations;
  violations = p_factory->getViolations();
  int nsize = violations.size();
  int i;
  char sbuf[128];
  for (i=0; i<nsize; i++) {
    std::string string;
    if (violations[i].bus_violation) {
      sprintf(sbuf,"     Bus voltage violation on bus %d",violations[i].bus1);
      string = sbuf;
    } else if (violations[i].line_violation) {
      sprintf(sbuf,"     Branch overload violation on branch [%d,%d] for line %s",
          violations[i].bus1,violations[i].bus2,violations[i].tag);
      string = sbuf;
    }
    if (!p_no_print) {
      printf("DEBUG VIOLATION: (%s)\n",string.c_str());
    }
    ret.push_back(string);
  }
  return ret;
}

/**
 * User rate B parameter for line overload violations
 * @param flag if true, use RATEB parameter
 */
void gridpack::voltage_stability::VSAppModule::useRateB(bool flag)
{
  p_factory->useRateB(flag);
}

/**
 * Suppress all output from power flow module
 * @param flag if true, suppress printing
 */
void gridpack::voltage_stability::VSAppModule::suppressOutput(bool flag)
{
  p_no_print = flag;
}

#ifdef USE_GOSS
/**
 * Set GOSS client if one already exists
 * @param client pointer to existing GOSS client
 */
void gridpack::voltage_stability::VSAppModule::setGOSSClient(
    gridpack::goss::GOSSClient &client, std:;string simID)
{
  p_goss_client = client;
  p_simID = simID;
}
#endif
/**
 * Modify generator parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param gen_id two character token specifying generator on bus
 * @param genParam string representing dictionary name of data element
 *                to be modified
 * @param value new value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::voltage_stability::VSAppModule::modifyDataCollectionGenParam(
    int bus_id, std::string gen_id, std::string genParam, double value)
{
  gridpack::utility::StringUtils util;
  std::string clean_gen_id = util.clean2Char(gen_id);
  return p_modifyDataCollectionGenParam<double>(bus_id,clean_gen_id,genParam,value);
}
bool gridpack::voltage_stability::VSAppModule::modifyDataCollectionGenParam(
    int bus_id, std::string gen_id, std::string genParam, int value)
{
  gridpack::utility::StringUtils util;
  std::string clean_gen_id = util.clean2Char(gen_id);
  return p_modifyDataCollectionGenParam<int>(bus_id,clean_gen_id,genParam,value);
}

/**
 * Modify load parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param load_id two character token specifying load on bus
 * @param loadParam string representing dictionary name of data element
 *                to be modified
 * @param value new value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::voltage_stability::VSAppModule::modifyDataCollectionLoadParam(
    int bus_id, std::string load_id, std::string loadParam, double value)
{
  gridpack::utility::StringUtils util;
  std::string clean_load_id = util.clean2Char(load_id);
  return p_modifyDataCollectionLoadParam<double>(bus_id,clean_load_id,loadParam,value);
}
bool gridpack::voltage_stability::VSAppModule::modifyDataCollectionLoadParam(
    int bus_id, std::string load_id, std::string loadParam, int value)
{
  gridpack::utility::StringUtils util;
  std::string clean_load_id = util.clean2Char(load_id);
  return p_modifyDataCollectionLoadParam<int>(bus_id,clean_load_id,loadParam,value);
}

/**
 * Modify parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param busParam string representing dictionary name of data element
 *                to be modified
 * @param value new value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::voltage_stability::VSAppModule::modifyDataCollectionBusParam(
    int bus_id, std::string busParam, double value)
{
  return p_modifyDataCollectionBusParam<double>(bus_id,busParam,value);
}
bool gridpack::voltage_stability::VSAppModule::modifyDataCollectionBusParam(
    int bus_id, std::string busParam, int value)
{
  return p_modifyDataCollectionBusParam<int>(bus_id,busParam,value);
}

/**
 * Modify parameters in data collection for specified branch
 * @param bus1, bus2 bus IDs for from and to bus
 * @param ckt two character token specifying branch
 * @param branchParam string representing dictionary name of data element
 *                to be modified
 * @param value new value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::voltage_stability::VSAppModule::modifyDataCollectionBranchParam(
    int bus1, int bus2, std::string ckt,
    std::string branchParam, double value)
{
  return p_modifyDataCollectionBranchParam<double>(bus1,bus2,ckt,branchParam,value);
}
bool gridpack::voltage_stability::VSAppModule::modifyDataCollectionBranchParam(
    int bus1, int bus2, std::string ckt,
    std::string branchParam, int value)
{
  return p_modifyDataCollectionBranchParam<int>(bus1,bus2,ckt,branchParam,value);
}

/**
 * Get generator parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param gen_id two character token specifying generator on bus
 * @param genParam string representing dictionary name of data element
 *                to be modified
 * @param value value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::voltage_stability::VSAppModule::getDataCollectionGenParam(
    int bus_id, std::string gen_id,
    std::string genParam, double *value)
{
  return p_getDataCollectionGenParam<double>(bus_id, gen_id, genParam, value);
}
bool gridpack::voltage_stability::VSAppModule::getDataCollectionGenParam(
    int bus_id, std::string gen_id,
    std::string genParam, int *value)
{
  return p_getDataCollectionGenParam<int>(bus_id, gen_id, genParam, value);
}

/**
 * Get load parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param load_id two character token specifying load on bus
 * @param loadParam string representing dictionary name of data element
 *                to be modified
 * @param value value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::voltage_stability::VSAppModule::getDataCollectionLoadParam(
    int bus_id, std::string load_id,
    std::string loadParam, double *value)
{
  return p_getDataCollectionLoadParam<double>(bus_id, load_id, loadParam, value);
}
bool gridpack::voltage_stability::VSAppModule::getDataCollectionLoadParam(
    int bus_id, std::string load_id,
    std::string loadParam, int *value)
{
  return p_getDataCollectionLoadParam<int>(bus_id, load_id, loadParam, value);
}

/**
 * Get parameters in data collection for specified bus
 * @param bus_id bus ID
 * @param busParam string representing dictionary name of data element
 *                to be modified
 * @param value value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::voltage_stability::VSAppModule::getDataCollectionBusParam(
    int bus_id, std::string busParam, double *value)
{
  return p_getDataCollectionBusParam<double>(bus_id, busParam, value);
}
bool gridpack::voltage_stability::VSAppModule::getDataCollectionBusParam(
    int bus_id, std::string busParam, int *value)
{
  return p_getDataCollectionBusParam<int>(bus_id, busParam, value);
}

/**
 * Get parameters in data collection for specified branch
 * @param bus1, bus2 bus IDs for from and to bus
 * @param ckt two character token specifying branch
 * @param branchParam string representing dictionary name of data element
 *                to be modified
 * @param value value of parameter
 * @return return false if parameter is not found
 */
bool gridpack::voltage_stability::VSAppModule::getDataCollectionBranchParam(
    int bus1, int bus2, std::string ckt,
    std::string branchParam, double *value)
{
  return p_getDataCollectionBranchParam<double>(bus1, bus2, ckt, branchParam, value);
}
bool gridpack::voltage_stability::VSAppModule::getDataCollectionBranchParam(
    int bus1, int bus2, std::string ckt,
    std::string branchParam, int *value)
{
  return p_getDataCollectionBranchParam<int>(bus1, bus2, ckt, branchParam, value);
}

/**
 * Check to see if PV Analysis is complete
 * @return return true if parameter is not found
 */
bool gridpack::voltage_stability::VSAppModule::isPVAnlyDone()
{
  if (current_increment >= (max_increment)){
	  p_bPVAnlyDone = true;
  }
	return p_bPVAnlyDone;
}

/**
 * Set up PV Curve internal parameters and initialize
 */
void gridpack::voltage_stability::VSAppModule::InitializePVCurve()
{
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Powerflow");
  std::string filename;
  double bus_mag = 0.0;
  gt = p_factory->getTotalGenRealPower(src_area,zone);
  lt = getTotalLoadRealPower(sink_area, zone);
  p_bPVAnlyDone = false;
  if (!cursor->get("PVAnalysisData",&filename)) {
     printf("No PV Analysis output data file specified\n");
  }
  else{
     static int numBus = p_network->numBuses();
     int i;
     double bus_varray[numBus] = {};
     std::ofstream file;
     file.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
     file.is_open();
     if(!PV_header){
       file << "Incremental Transfer (MW),";
       for (i=0; i<numBus; i++) {
         if (p_network->getActiveBus(i)) {
           gridpack::voltage_stability::VSBus *bus =
             dynamic_cast<gridpack::voltage_stability::VSBus*>
             (p_network->getBus(i).get());
           file << "Bus " << bus->getOriginalIndex() << ",";
         }
       }
     file << "\n";
     PV_header = true;
     file.close();
    }
  }
}

/**
 * Execute one transfer increment
 */
void gridpack::voltage_stability::VSAppModule::IncrementPVCurveStep()
{
  IncrementGeneratorRealPower(increment, src_area, zone, gt);
  IncrementLoadPower(increment, sink_area, zone, lt);
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Powerflow");
  std::string filename;
  double bus_mag = 0.0;
  if (!cursor->get("PVAnalysisData",&filename)) {
     printf("No PV Analysis output data file specified\n");
  }
  else{
     static int numBus = p_network->numBuses();
     int i;
     double bus_varray[numBus] = {};
     std::ofstream file;
     file.open(filename.c_str(), std::ios_base::app);
     file.is_open();
     file << current_increment << ",";
     for (i=0; i<numBus; i++) {
       if (p_network->getActiveBus(i)) {
         gridpack::voltage_stability::VSBus *bus =
           dynamic_cast<gridpack::voltage_stability::VSBus*>
           (p_network->getBus(i).get());
         file << bus->getVoltage() << ",";
       }
     }
    file << "\n";
    file.close();
    
  }
  gt = gt + increment;
  lt = lt + increment;
  current_increment = current_increment + increment;
}
