/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_app.cpp
 * @author Bruce Palmer
 * @date   2015-05-06 12:51:38 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "pf_app_module.hpp"
#include "pf_factory_module.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/parser/PTI33_parser.hpp"
#include "gridpack/parser/GOSS_parser.hpp"
#include "gridpack/math/math.hpp"
#include "pf_helper.hpp"

#define USE_REAL_VALUES

/**
 * Basic constructor
 */
gridpack::powerflow::PFAppModule::PFAppModule(void)
{
}

/**
 * Basic destructor
 */
gridpack::powerflow::PFAppModule::~PFAppModule(void)
{
}

enum Parser{PTI23, PTI33, GOSS};

/**
 * Read in and partition the powerflow network. The input file is read
 * directly from the Powerflow block in the configuration file so no
 * external file names or parameters need to be passed to this routine
 * @param network pointer to a PFNetwork object. This should not have any
 * buses or branches defined on it.
 * @param config point to open configuration file
 */
void gridpack::powerflow::PFAppModule::readNetwork(
    boost::shared_ptr<PFNetwork> &network,
    gridpack::utility::Configuration *config)
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
  if (!cursor->get("networkConfiguration",&filename)) {
    if (cursor->get("networkConfiguration_v33",&filename)) {
      filetype = PTI33;
    } else if (cursor->get("networkConfiguration_GOSS",&filename)) {
      filetype = GOSS;
    } else {
      printf("No network configuration file specified\n");
      return;
    }
  }
  // Convergence and iteration parameters
  p_tolerance = cursor->get("tolerance",1.0e-6);
  p_max_iteration = cursor->get("maxIteration",50);
  ComplexType tol;
  // Phase shift sign
  double phaseShiftSign = cursor->get("phaseShiftSign",1.0);

  int t_pti = timer->createCategory("Powerflow: Network Parser");
  timer->start(t_pti);
  if (filetype == PTI23) {
    gridpack::parser::PTI23_parser<PFNetwork> parser(network);
    parser.parse(filename.c_str());
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  } else if (filetype == PTI33) {
    gridpack::parser::PTI33_parser<PFNetwork> parser(network);
    parser.parse(filename.c_str());
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  } else if (filetype == GOSS) {
    gridpack::parser::GOSS_parser<PFNetwork> parser(network);
    parser.parse(filename.c_str());
  }
  timer->stop(t_pti);

  // Create serial IO object to export data from buses
  p_busIO.reset(new gridpack::serial_io::SerialBusIO<PFNetwork>(512,network));

  // Create serial IO object to export data from branches
  p_branchIO.reset(new gridpack::serial_io::SerialBranchIO<PFNetwork>(512,network));
  char ioBuf[128];

  sprintf(ioBuf,"\nMaximum number of iterations: %d\n",p_max_iteration);
  p_busIO->header(ioBuf);
  sprintf(ioBuf,"\nConvergence tolerance: %f\n",p_tolerance);
  p_busIO->header(ioBuf);

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
void gridpack::powerflow::PFAppModule::initialize()
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);

  // create factory
  p_factory.reset(new gridpack::powerflow::PFFactoryModule(p_network));
  int t_load = timer->createCategory("Powerflow: Factory Load");
  timer->start(t_load);
  p_factory->load();
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
void gridpack::powerflow::PFAppModule::reload()
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
bool gridpack::powerflow::PFAppModule::solve()
{
  bool ret = true;
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);

  // set YBus components so that you can create Y matrix
  int t_fact = timer->createCategory("Powerflow: Factory Operations");
  timer->start(t_fact);
  p_factory->setYBus();
  timer->stop(t_fact);

  int t_cmap = timer->createCategory("Powerflow: Create Mappers");
  timer->start(t_cmap);
  p_factory->setMode(YBus); 
#if 0
  gridpack::mapper::FullMatrixMap<PFNetwork> mMap(p_network);
#endif
  timer->stop(t_cmap);
  int t_mmap = timer->createCategory("Powerflow: Map to Matrix");
  timer->start(t_mmap);
#if 0
  boost::shared_ptr<gridpack::math::Matrix> Y = mMap.mapToMatrix();
  p_busIO->header("\nY-matrix values\n");
  Y->print();
//  Y->save("Ybus.m");
#endif
  timer->stop(t_mmap);

  timer->start(t_fact);
  p_factory->setMode(S_Cal);
  timer->stop(t_fact);
  timer->start(t_cmap);
  gridpack::mapper::BusVectorMap<PFNetwork> vvMap(p_network);
  timer->stop(t_cmap);
  int t_vmap = timer->createCategory("Powerflow: Map to Vector");

  // make Sbus components to create S vector
  timer->start(t_fact);
  p_factory->setSBus();
  timer->stop(t_fact);
//  p_busIO->header("\nIteration 0\n");

  // Set PQ
  timer->start(t_cmap);
  p_factory->setMode(RHS); 
  gridpack::mapper::BusVectorMap<PFNetwork> vMap(p_network);
  timer->stop(t_cmap);
  timer->start(t_vmap);
#ifdef USE_REAL_VALUES
  boost::shared_ptr<gridpack::math::RealVector> PQ = vMap.mapToRealVector();
#else
  boost::shared_ptr<gridpack::math::Vector> PQ = vMap.mapToVector();
#endif
  timer->stop(t_vmap);
//  PQ->print();
  timer->start(t_cmap);
  p_factory->setMode(Jacobian);
  gridpack::mapper::FullMatrixMap<PFNetwork> jMap(p_network);
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

  gridpack::ComplexType tol = 2.0*p_tolerance;
  int iter = 0;

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
    p_busIO->header("Solver failure\n\n");
    timer->stop(t_lsolv);
    return false;
  }
  timer->stop(t_lsolv);
  tol = PQ->normInfinity();

  // Create timer for map to bus
  int t_bmap = timer->createCategory("Powerflow: Map to Bus");
  int t_updt = timer->createCategory("Powerflow: Bus Update");
  char ioBuf[128];

  while (real(tol) > p_tolerance && iter < p_max_iteration) {
    // Push current values in X vector back into network components
    // Need to implement setValues method in PFBus class in order for this to
    // work
    timer->start(t_bmap);
    p_factory->setMode(RHS);
    vMap.mapToBus(X);
    timer->stop(t_bmap);

    // Exchange data between ghost buses (I don't think we need to exchange data
    // between branches)
    timer->start(t_updt);
    p_network->updateBuses();
    timer->stop(t_updt);

    // Create new versions of Jacobian and PQ vector
    timer->start(t_vmap);
#ifdef USE_REAL_VALUES
    vMap.mapToRealVector(PQ);
#else
    vMap.mapToVector(PQ);
#endif
//    p_busIO->header("\nnew PQ vector\n");
//    PQ->print();
    timer->stop(t_vmap);
    timer->start(t_mmap);
    p_factory->setMode(Jacobian);
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
      p_busIO->header("Solver failure\n\n");
      timer->stop(t_lsolv);
      return false;
    }
    timer->stop(t_lsolv);

    tol = PQ->normInfinity();
    sprintf(ioBuf,"\nIteration %d Tol: %12.6e\n",iter+1,real(tol));
    p_busIO->header(ioBuf);
    iter++;
  }

  if (iter >= p_max_iteration) ret = false;

  // Push final result back onto buses
  timer->start(t_bmap);
  p_factory->setMode(RHS);
  vMap.mapToBus(X);
  timer->stop(t_bmap);

  // Make sure that ghost buses have up-to-date values before printing out
  // results
  timer->start(t_updt);
  p_network->updateBuses();
  timer->stop(t_updt);
  timer->stop(t_total);
  return ret;
}
/**
 * Execute the iterative solve portion of the application using a library
 * non-linear solver
 * @return false if an error was caught in the solution algorithm
 */
bool gridpack::powerflow::PFAppModule::nl_solve()
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);

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

  PFSolverHelper helper(p_factory, p_network);
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
    ret = false;
  }

  timer->stop(t_lsolv);
  timer->stop(t_total);
  return ret;
}


/**
 * Write out results of powerflow calculation to standard output or a file
 */
void gridpack::powerflow::PFAppModule::write()
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);
  int t_write = timer->createCategory("Powerflow: Write Results");
  timer->start(t_write);
  p_branchIO->header("\n   Branch Power Flow\n");
  p_branchIO->header("\n        Bus 1       Bus 2   CKT         P"
                  "                    Q\n");
  p_branchIO->write();


  p_busIO->header("\n   Generator Power\n");
  p_busIO->header("\n   Bus Number  GenID        Pgen              Qgen\n");
  p_busIO->write("power");
  p_busIO->header("\n   Bus Voltages and Phase Angles\n");
  p_busIO->header("\n   Bus Number      Phase Angle      Voltage Magnitude\n");
  p_busIO->write();
  timer->stop(t_write);
  timer->stop(t_total);
}

void gridpack::powerflow::PFAppModule::writeBus(const char *signal)
{
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

void gridpack::powerflow::PFAppModule::writeBranch(const char *signal)
{
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

void gridpack::powerflow::PFAppModule::writeHeader(const char *msg)
{
  p_busIO->header(msg);
}

/**
 * Redirect output from standard out
 * @param filename name of file to write results to
 */
void gridpack::powerflow::PFAppModule::open(const char *filename)
{
  p_busIO->open(filename);
  p_branchIO->setStream(p_busIO->getStream());
}

void gridpack::powerflow::PFAppModule::close()
{
  p_busIO->close();
  p_branchIO->setStream(p_busIO->getStream());
}

/**
 * Print string. This can be used to direct output to the file opened using
 * the open command
 * @param buf string to be printed
 */
void gridpack::powerflow::PFAppModule::print(const char *buf)
{
  p_busIO->header(buf);
}

/**
 * Save results of powerflow calculation to data collection objects
 */
void gridpack::powerflow::PFAppModule::saveData()
{
  p_factory->saveData();
}

/**
 * Set a contingency
 * @param event data describing location and type of contingency
 * @return false if location of contingency is not found in
 * network
 */
bool gridpack::powerflow::PFAppModule::setContingency(
    gridpack::powerflow::Contingency &event)
{
  bool ret = true;
  if (event.p_type == Generator) {
    int ngen = event.p_busid.size();
    int i, j, idx, jdx;
    for (i=0; i<ngen; i++) {
      idx = event.p_busid[i];
      std::string tag = event.p_genid[i];
      std::vector<int> lids = p_network->getLocalBusIndices(idx);
      if (lids.size() == 0) ret = false;
      gridpack::powerflow::PFBus *bus;
      for (j=0; j<lids.size(); j++) {
        jdx = lids[j];
        bus = dynamic_cast<gridpack::powerflow::PFBus*>(
            p_network->getBus(jdx).get());
        event.p_saveGenStatus[i] = bus->getGenStatus(tag);
        bus->setGenStatus(tag, false);
      }
    }
  } else if (event.p_type == Branch) {
    int to, from;
    int nline = event.p_to.size();
    int i, j, idx, jdx;
    for (i=0; i<nline; i++) {
      to = event.p_to[i];
      from = event.p_from[i];
      std::string tag = event.p_ckt[i];
      std::vector<int> lids = p_network->getLocalBranchIndices(from,to);
      if (lids.size() == 0) ret = false;
      gridpack::powerflow::PFBranch *branch;
      for (j=0; j<lids.size(); j++) {
        jdx = lids[j];
        branch = dynamic_cast<gridpack::powerflow::PFBranch*>(
            p_network->getBranch(jdx).get());
        event.p_saveLineStatus[i] = branch->getBranchStatus(tag);
        branch->setBranchStatus(tag, false);
      }
    }
  } else {
    ret = false;
  }
  p_factory->checkLoneBus();
  return ret;
}

/**
 * Return system to the state before the contingency
 * @param event data describing location and type of contingency
 * @return false if location of contingency is not found in network
 */
bool gridpack::powerflow::PFAppModule::unSetContingency(
    gridpack::powerflow::Contingency &event)
{
  bool ret = true;
  if (event.p_type == Generator) {
    int ngen = event.p_busid.size();
    int i, j, idx, jdx;
    for (i=0; i<ngen; i++) {
      idx = event.p_busid[i];
      std::string tag = event.p_genid[i];
      std::vector<int> lids = p_network->getLocalBusIndices(idx);
      if (lids.size() == 0) ret = false;
      gridpack::powerflow::PFBus *bus;
      for (j=0; j<lids.size(); j++) {
        jdx = lids[j];
        bus = dynamic_cast<gridpack::powerflow::PFBus*>(
            p_network->getBus(jdx).get());
        bus->setGenStatus(tag, event.p_saveGenStatus[i]);
      }
    }
  } else if (event.p_type == Branch) {
    int to, from;
    int nline = event.p_to.size();
    int i, j, idx, jdx;
    for (i=0; i<nline; i++) {
      to = event.p_to[i];
      from = event.p_from[i];
      std::string tag = event.p_ckt[i];
      std::vector<int> lids = p_network->getLocalBranchIndices(from,to);
      if (lids.size() == 0) ret = false;
      gridpack::powerflow::PFBranch *branch;
      for (j=0; j<lids.size(); j++) {
        jdx = lids[j];
        branch = dynamic_cast<gridpack::powerflow::PFBranch*>(
            p_network->getBranch(jdx).get());
        branch->setBranchStatus(tag,event.p_saveLineStatus[i]);
      }
    }
  } else {
    ret = false;
  }
  p_factory->clearLoneBus();
  return ret;
}

/**
 * Check to see if there are any voltage violations in the network
 * @param area area number. If this parameter is included, only check for
 * violations in this area
 * @param minV maximum voltage limit
 * @param maxV maximum voltage limit
 * @return true if no violations found
 */
bool gridpack::powerflow::PFAppModule::checkVoltageViolations(
  double Vmin, double Vmax)
{
  return p_factory->checkVoltageViolations(Vmin,Vmax);
}
bool gridpack::powerflow::PFAppModule::checkVoltageViolations(
 int area,    double Vmin, double Vmax)
{
  return p_factory->checkVoltageViolations(area, Vmin,Vmax);
}

/**
 * Set "ignore" parameter on all buses with violations so that subsequent
 * checks are not counted as violations
 * @param minV maximum voltage limit
 * @param maxV maximum voltage limit
 */
void gridpack::powerflow::PFAppModule::ignoreVoltageViolations(double Vmin,
    double Vmax)
{
  p_factory->ignoreVoltageViolations(Vmin,Vmax);
}

/**
 * Clear "ignore" parameter on all buses
 */
void gridpack::powerflow::PFAppModule::clearVoltageViolations()
{
  p_factory->clearVoltageViolations();
}

/**
 * Check to see if there are any line overload violations in the
 * network
 * @param area area number. If this parameter is included, only check for
 * violations in this area
 * @return true if no violations found
 */
bool gridpack::powerflow::PFAppModule::checkLineOverloadViolations()
{
  return p_factory->checkLineOverloadViolations();
}
bool gridpack::powerflow::PFAppModule::checkLineOverloadViolations(int area)
{
  return p_factory->checkLineOverloadViolations(area);
}

/**
 * Reset voltages to values in network configuration file
 */
void gridpack::powerflow::PFAppModule::resetVoltages()
{
  p_factory->resetVoltages();
}
