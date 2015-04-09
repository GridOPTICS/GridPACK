/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_app.cpp
 * @author Bruce Palmer
 * @date   2014-01-28 11:39:16 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/include/gridpack.hpp"
#include "pf_app_module.hpp"
#include "pf_factory_module.hpp"

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

  int t_pti = timer->createCategory("Powerflow: Network Parser");
  timer->start(t_pti);
  if (filetype == PTI23) {
    gridpack::parser::PTI23_parser<PFNetwork> parser(network);
    parser.parse(filename.c_str());
  } else if (filetype == PTI33) {
    gridpack::parser::PTI33_parser<PFNetwork> parser(network);
    parser.parse(filename.c_str());
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
 * Execute the iterative solve portion of the application
 */
void gridpack::powerflow::PFAppModule::solve()
{
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
  p_busIO->header("\nIteration 0\n");

  // Set PQ
  timer->start(t_cmap);
  p_factory->setMode(RHS); 
  gridpack::mapper::BusVectorMap<PFNetwork> vMap(p_network);
  timer->stop(t_cmap);
  timer->start(t_vmap);
  boost::shared_ptr<gridpack::math::Vector> PQ = vMap.mapToVector();
  timer->stop(t_vmap);
//  PQ->print();
  timer->start(t_cmap);
  p_factory->setMode(Jacobian);
  gridpack::mapper::FullMatrixMap<PFNetwork> jMap(p_network);
  timer->stop(t_cmap);
  timer->start(t_mmap);
  boost::shared_ptr<gridpack::math::Matrix> J = jMap.mapToMatrix();
  timer->stop(t_mmap);
//  p_busIO->header("\nJacobian values\n");
//  J->print();

  // Create X vector by cloning PQ
  boost::shared_ptr<gridpack::math::Vector> X(PQ->clone());

  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Powerflow");
  // Create linear solver
  int t_csolv = timer->createCategory("Powerflow: Create Linear Solver");
  timer->start(t_csolv);
  gridpack::math::LinearSolver solver(*J);
  solver.configure(cursor);
  timer->stop(t_csolv);

  gridpack::ComplexType tol = 2.0*p_tolerance;
  int iter = 0;

  // First iteration
  X->zero(); //might not need to do this
  p_busIO->header("\nCalling solver\n");
  int t_lsolv = timer->createCategory("Powerflow: Solve Linear Equation");
  timer->start(t_lsolv);
//    char dbgfile[32];
//    sprintf(dbgfile,"j0.bin");
//    J->saveBinary(dbgfile);
//    sprintf(dbgfile,"pq0.bin");
//    PQ->saveBinary(dbgfile);
  solver.solve(*PQ, *X);
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
    vMap.mapToVector(PQ);
//    p_busIO->header("\nnew PQ vector\n");
//    PQ->print();
    timer->stop(t_vmap);
    timer->start(t_mmap);
    p_factory->setMode(Jacobian);
    jMap.mapToMatrix(J);
    timer->stop(t_mmap);

    // Create linear solver
    timer->start(t_lsolv);
    X->zero(); //might not need to do this
//    sprintf(dbgfile,"j%d.bin",iter+1);
//    J->saveBinary(dbgfile);
//    sprintf(dbgfile,"pq%d.bin",iter+1);
//    PQ->saveBinary(dbgfile);
    solver.solve(*PQ, *X);
    timer->stop(t_lsolv);

    tol = PQ->normInfinity();
    sprintf(ioBuf,"\nIteration %d Tol: %12.6e\n",iter+1,real(tol));
    p_busIO->header(ioBuf);
    iter++;
  }

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
}

/**
 * Execute the iterative solve portion of the application
 */
/* void gridpack::powerflow::PFAppModule::solve_step1()
{
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
  timer->stop(t_total);
}

//
// solve for RTPT
//
void gridpack::powerflow::PFAppModule::solve_updatePg(std::vector<pathStress> pstress, std::vector<double> p_slice3Values)
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);

  if (pstress.size() == 3) {
    for (int i = 0; i < p_slice3Values.size(); i++) {
      p_factory.reset(new gridpack::powerflow::PFFactoryModule(p_network));
      double deltaP3 = p_slice3Values[i];
      for (int j = 0; j < pstress[2].sourceArea.size(); j++) {
        int ps_busID = pstress[2].sourceArea[j].busID;
        std::string ps_genID = pstress[2].sourceArea[j].genID;
        double ps_pg = deltaP3 * pstress[2].sourceArea[j].participation;
  	p_factory->updatePg(ps_busID, ps_genID, ps_pg);
      }
      for (int j = 0; j < pstress[2].sinkArea.size(); j++) {
        int ps_busID = pstress[2].sinkArea[j].busID;
        std::string ps_genID = pstress[2].sinkArea[j].genID;
        double ps_pg = -deltaP3 * pstress[2].sinkArea[j].participation;
  	p_factory->updatePg(ps_busID, ps_genID, ps_pg);
      }
      solve_step2();
      write();
    }
  }
}

//
// solve for RTPT
//
void gridpack::powerflow::PFAppModule::solve_step2()
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Powerflow: Total Application");
  timer->start(t_total);

  // make Sbus components to create S vector
  int t_fact = timer->createCategory("Powerflow: Factory Operations");
  timer->start(t_fact);
  p_factory->setSBus();
  timer->stop(t_fact);
  p_busIO->header("\nIteration 0\n");

  // Set PQ
  int t_cmap = timer->createCategory("Powerflow: Create Mappers");
  timer->start(t_cmap);
  p_factory->setMode(RHS); 
  gridpack::mapper::BusVectorMap<PFNetwork> vMap(p_network);
  int t_vmap = timer->createCategory("Powerflow: Map to Vector");
  timer->stop(t_cmap);
  timer->start(t_vmap);
  boost::shared_ptr<gridpack::math::Vector> PQ = vMap.mapToVector();
  timer->stop(t_vmap);
//  PQ->print();
  timer->start(t_cmap);
  p_factory->setMode(Jacobian);
  gridpack::mapper::FullMatrixMap<PFNetwork> jMap(p_network);
  timer->stop(t_cmap);
  int t_mmap = timer->createCategory("Powerflow: Map to Matrix");
  timer->start(t_mmap);
  boost::shared_ptr<gridpack::math::Matrix> J = jMap.mapToMatrix();
  timer->stop(t_mmap);
//  p_busIO->header("\nJacobian values\n");
//  J->print();

  // Create X vector by cloning PQ
  boost::shared_ptr<gridpack::math::Vector> X(PQ->clone());

  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = p_config->getCursor("Configuration.Powerflow");
  // Create linear solver
  int t_csolv = timer->createCategory("Powerflow: Create Linear Solver");
  timer->start(t_csolv);
  gridpack::math::LinearSolver solver(*J);
  solver.configure(cursor);
  timer->stop(t_csolv);

  gridpack::ComplexType tol = 2.0*p_tolerance;
  int iter = 0;

  // First iteration
  X->zero(); //might not need to do this
  p_busIO->header("\nCalling solver\n");
  int t_lsolv = timer->createCategory("Powerflow: Solve Linear Equation");
  timer->start(t_lsolv);
//    char dbgfile[32];
//    sprintf(dbgfile,"j0.bin");
//    J->saveBinary(dbgfile);
//    sprintf(dbgfile,"pq0.bin");
//    PQ->saveBinary(dbgfile);
  solver.solve(*PQ, *X);
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
    vMap.mapToVector(PQ);
//    p_busIO->header("\nnew PQ vector\n");
//    PQ->print();
    timer->stop(t_vmap);
    timer->start(t_mmap);
    p_factory->setMode(Jacobian);
    jMap.mapToMatrix(J);
    timer->stop(t_mmap);

    // Create linear solver
    timer->start(t_lsolv);
    X->zero(); //might not need to do this
//    sprintf(dbgfile,"j%d.bin",iter+1);
//    J->saveBinary(dbgfile);
//    sprintf(dbgfile,"pq%d.bin",iter+1);
//    PQ->saveBinary(dbgfile);
    solver.solve(*PQ, *X);
    timer->stop(t_lsolv);

    tol = PQ->normInfinity();
    sprintf(ioBuf,"\nIteration %d Tol: %12.6e\n",iter+1,real(tol));
    p_busIO->header(ioBuf);
    iter++;
  }

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
}
*/

/**
 * Write out results of powerflow calculation to standard output
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


  p_busIO->header("\n   Bus Voltages and Phase Angles\n");
  p_busIO->header("\n   Bus Number      Phase Angle      Voltage Magnitude\n");
  p_busIO->write();
  timer->stop(t_write);
  timer->stop(t_total);
}

/**
 * Save results of powerflow calculation to data collection objects
 */
void gridpack::powerflow::PFAppModule::saveData()
{
  p_factory->saveData();
}
