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
#include "pf_app.hpp"
#include "pf_factory.hpp"


// Calling program for powerflow application

/**
 * Basic constructor
 */
gridpack::powerflow::PFApp::PFApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::powerflow::PFApp::~PFApp(void)
{
}

enum Parser{PTI23, PTI33};

/**
 * Execute application
 * @param argc number of arguments
 * @param argv list of character strings
 */
void gridpack::powerflow::PFApp::execute(int argc, char** argv)
{
  // load input file
  gridpack::parallel::Communicator world;
  boost::shared_ptr<PFNetwork> network(new PFNetwork(world));

  // read configuration file
  gridpack::utility::Configuration *config
    = gridpack::utility::Configuration::configuration();
  config->enableLogging(&std::cout);
  bool opened;
  if (argc >= 2 && argv[1] != NULL) {
    char inputfile[256];
    sprintf(inputfile,"%s",argv[1]);
    opened = config->open(inputfile,world);
  } else {
    opened = config->open("input.xml",world);
  }
  if (!opened) return;

  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Powerflow");
  std::string filename;
  int filetype = PTI23;
  if (!cursor->get("networkConfiguration",&filename)) {
    if (cursor->get("networkConfiguration_v33",&filename)) {
      filetype = PTI33;
    } else {
      printf("No network configuration file specified\n");
      return;
    }
  }
  // Convergence and iteration parameters
  double tolerance = cursor->get("tolerance",1.0e-6);
  int max_iteration = cursor->get("maxIteration",50);
  ComplexType tol;
  // Phase shift sign
  double phaseShiftSign = cursor->get("phaseShiftSign",1.0);

  printf("Filename: (%s)\n",filename.c_str());
  if (filetype == PTI23) {
    printf("Using V23 parser\n");
    gridpack::parser::PTI23_parser<PFNetwork> parser(network);
    parser.parse(filename.c_str());
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  } else if (filetype == PTI33) {
    printf("Using V33 parser\n");
    gridpack::parser::PTI33_parser<PFNetwork> parser(network);
    parser.parse(filename.c_str());
    if (phaseShiftSign == -1.0) {
      parser.changePhaseShiftSign();
    }
  }
  printf("NBUS: %d NBRANCH: %d\n",network->numBuses(),
      network->numBranches());

  // Create serial IO object to export data from buses
  gridpack::serial_io::SerialBusIO<PFNetwork> busIO(8192,network);
  char ioBuf[128];

  sprintf(ioBuf,"\nMaximum number of iterations: %d\n",max_iteration);
  busIO.header(ioBuf);
  sprintf(ioBuf,"\nConvergence tolerance: %f\n",tolerance);
  busIO.header(ioBuf);

  // partition network
  network->partition();

  // create factory
  gridpack::powerflow::PFFactory factory(network);
  factory.load();

  // set network components using factory
  factory.setComponents();

  // Set up bus data exchange buffers. Need to decide what data needs to be
  // exchanged
  factory.setExchange();

  // Create bus data exchange
  network->initBusUpdate();

  // set YBus components so that you can create Y matrix
  factory.setYBus();

  factory.setMode(YBus); 

  factory.setMode(S_Cal);
  gridpack::mapper::BusVectorMap<PFNetwork> vvMap(network);

  // make Sbus components to create S vector
  factory.setSBus();
  busIO.header("\nIteration 0\n");

  // Set PQ
  factory.setMode(RHS); 
  gridpack::mapper::BusVectorMap<PFNetwork> vMap(network);
  boost::shared_ptr<gridpack::math::Vector> PQ = vMap.mapToVector();
  factory.setMode(Jacobian);
  gridpack::mapper::FullMatrixMap<PFNetwork> jMap(network);
  boost::shared_ptr<gridpack::math::Matrix> J = jMap.mapToMatrix();

  // Create X vector by cloning PQ
  boost::shared_ptr<gridpack::math::Vector> X(PQ->clone());

  // Create linear solver
  gridpack::math::LinearSolver solver(*J);
  solver.configure(cursor);

  tol = 2.0*tolerance;
  int iter = 0;

  // First iteration
  X->zero(); //might not need to do this
  busIO.header("\nCalling solver\n");
  solver.solve(*PQ, *X);
  tol = PQ->normInfinity();


  while (real(tol) > tolerance && iter < max_iteration) {
    // Push current values in X vector back into network components
    // Need to implement setValues method in PFBus class in order for this to
    // work
    factory.setMode(RHS);
    vMap.mapToBus(X);

    // Exchange data between ghost buses (We don't need to exchange data
    // between branches)
    network->updateBuses();

    // Create new versions of Jacobian and PQ vector
    vMap.mapToVector(PQ);
    factory.setMode(Jacobian);
    jMap.mapToMatrix(J);

    // Create linear solver
    X->zero(); //might not need to do this
    solver.solve(*PQ, *X);

    tol = PQ->normInfinity();
    sprintf(ioBuf,"\nIteration %d Tol: %12.6e\n",iter+1,real(tol));
    busIO.header(ioBuf);
    iter++;
  }

  // Push final result back onto buses
  factory.setMode(RHS);
  vMap.mapToBus(X);

  // Make sure that ghost buses have up-to-date values before printing out
  // results
  network->updateBuses();

  gridpack::serial_io::SerialBranchIO<PFNetwork> branchIO(512,network);
  branchIO.header("\n   Branch Power Flow\n");
  branchIO.header("\n        Bus 1       Bus 2   CKT         P"
                  "                    Q\n");
  branchIO.write();


  busIO.header("\n   Bus Voltages and Phase Angles\n");
  busIO.header("\n   Bus Number      Phase Angle      Voltage Magnitude\n");
  busIO.write();
}
