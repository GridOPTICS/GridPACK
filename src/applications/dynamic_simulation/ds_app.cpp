/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_app.cpp
 * @author Shuangshuang Jin
 * @date   2015-03-06 14:48:54 d3g096
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/dynamic_simulation/ds_app.hpp"

// Calling program for dynamic simulation application

/**
 * Basic constructor
 */
gridpack::dynamic_simulation::DSApp::DSApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::DSApp::~DSApp(void)
{
}

/**
 * Execute application
 */
void gridpack::dynamic_simulation::DSApp::execute(int argc, char** argv)
{
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Total Application");
  timer->start(t_total);
  gridpack::parallel::Communicator world;
  boost::shared_ptr<DSNetwork> network(new DSNetwork(world));

  int t_setup = timer->createCategory("Read in Data");
  timer->start(t_setup);
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
  cursor = config->getCursor("Configuration.Dynamic_simulation");
  std::string filename = cursor->get("networkConfiguration",
      "No network configuration specified");
  double sim_time = cursor->get("simulationTime",0.0);
  if (sim_time == 0.0) {
    // TODO: some kind of error
  }
  double time_step = cursor->get("timeStep",0.0);
  if (time_step == 0.0) {
    // TODO: some kind of error
  }

  // Read in information about fault events and store them in internal data
  // structure
  cursor = config->getCursor("Configuration.Dynamic_simulation.faultEvents");
  gridpack::utility::Configuration::ChildCursors events;
  if (cursor) cursor->children(events);
  std::vector<gridpack::dynamic_simulation::DSBranch::Event>
     faults = setFaultEvents(events); 

  // load input file
  gridpack::parser::PTI23_parser<DSNetwork> parser(network);
  parser.parse(filename.c_str());
  cursor = config->getCursor("Configuration.Dynamic_simulation");
  filename = cursor->get("generatorParameters","");
#if 0
  if (filename.size() > 0) parser.parse(filename.c_str());
#endif
  timer->stop(t_setup);

  int t_part = timer->createCategory("Partition Network");
  timer->start(t_part);
  // partition network
  network->partition();
  timer->stop(t_part);
#if 1
  if (filename.size() > 0) parser.externalParse(filename.c_str());
#endif

  // Create serial IO object to export data from buses or branches
  gridpack::serial_io::SerialBusIO<DSNetwork> busIO(2048, network);
  gridpack::serial_io::SerialBranchIO<DSNetwork> branchIO(128, network);
  char ioBuf[128];

  int t_config = timer->createCategory("Configure Network");
  timer->start(t_config);
  // create factory
  gridpack::dynamic_simulation::DSFactory factory(network);
  factory.load();

  // set network components using factory
  factory.setComponents();
  timer->stop(t_config);

  // set YBus components so that you can create Y matrix  
  factory.setYBus();

  if (!factory.checkGen()) {
    busIO.header("Missing generators on at least one processor\n");
    return;
  }

  int t_matset = timer->createCategory("Matrix Setup");
  timer->start(t_matset);
  factory.setMode(YBUS);
  gridpack::mapper::FullMatrixMap<DSNetwork> ybusMap(network);
#if 0
  boost::shared_ptr<gridpack::math::Matrix> orgYbus = ybusMap.mapToMatrix();
  ///branchIO.header("\n=== orginal ybus: ============\n");
  ///orgYbus->print();

  // Form constant impedance load admittance yl for all buses and add it to 
  // system Y matrix: ybus = ybus + yl
  factory.setMode(YL);
  boost::shared_ptr<gridpack::math::Matrix> ybus = ybusMap.mapToMatrix();
  ///branchIO.header("\n=== ybus after added yl: ============\n");
  ///ybus->print();
 
  // Construct permutation matrix perm by checking for multiple generators at a bus
  factory.setMode(PERM);
  gridpack::mapper::FullMatrixMap<DSNetwork> permMap(network);
  boost::shared_ptr<gridpack::math::Matrix> perm = permMap.mapToMatrix();
#endif
  timer->stop(t_matset);
  ///busIO.header("\n=== perm: ============\n");
  ///perm->print(); 

  // Form a transposed matrix of perm
  int t_trans = timer->createCategory("Matrix Transpose");
#if 0
  timer->start(t_trans);
   boost::shared_ptr<gridpack::math::Matrix> permTrans(transpose(*perm));
  timer->stop(t_trans);
  ///busIO.header("\n=== permTrans: ============\n");
  ///permTrans->print();
#endif

  // Construct matrix Y_a using extracted xd and ra from gen data, 
  // and construct its diagonal matrix diagY_a
  timer->start(t_matset);
  factory.setMode(YA);
  gridpack::mapper::FullMatrixMap<DSNetwork> yaMap(network);
  boost::shared_ptr<gridpack::math::Matrix> diagY_a = yaMap.mapToMatrix();
  ///busIO.header("\n=== diagY_a: ============\n");
  ///diagY_a->print(); 
  // Convert diagY_a from sparse matrix to dense matrix Y_a so that SuperLU_DIST can solve
  gridpack::math::MatrixStorageType denseType = gridpack::math::Dense;
  boost::shared_ptr<gridpack::math::Matrix> Y_a(gridpack::math::storageType(*diagY_a, denseType));
  timer->stop(t_matset);

  // Construct matrix Ymod: Ymod = diagY_a * permTrans
  int t_matmul = timer->createCategory("Matrix Multiply");
#if 0
  timer->start(t_matmul);
  boost::shared_ptr<gridpack::math::Matrix> Ymod(multiply(*diagY_a, *permTrans));
  timer->stop(t_matmul);
#endif
  ///busIO.header("\n=== Ymod: ============\n");
  ///Ymod->print(); 
 
  // Form matrix Y_b: Y_b(1:ngen, jg) = -Ymod, where jg represents the 
  // corresponding index sets of buses that the generators are connected to. 
  // Then construct Y_b's transposed matrix Y_c: Y_c = Y_b'
  // This two steps can be done by using a P matrix to get Y_c directly.
  timer->start(t_matset);
#if 0
  factory.setMode(PMatrix);
  gridpack::mapper::FullMatrixMap<DSNetwork> pMap(network);
  boost::shared_ptr<gridpack::math::Matrix> P = pMap.mapToMatrix();
#endif
  ///busIO.header("\n=== P: ============\n");
  ///P->print();
  factory.setMode(YC);
  gridpack::mapper::FullMatrixMap<DSNetwork> cMap(network);
  boost::shared_ptr<gridpack::math::Matrix> Y_cDense = cMap.mapToMatrix(true);
  //Y_c->scale(-1.0);
  ///busIO.header("\n=== Y_c: ============\n");
  ///Y_c->print();
  factory.setMode(YB);
  gridpack::mapper::FullMatrixMap<DSNetwork> bMap(network);
  boost::shared_ptr<gridpack::math::Matrix> Y_b = bMap.mapToMatrix();
  timer->stop(t_matset);
  ///busIO.header("\n=== Y_b: ============\n");
  ///Y_b->print();
  
  //Y_c->scale(-1.0); // scale Y_c by -1.0 for linear solving
  // Convert Y_c from sparse matrix to dense matrix Y_cDense so that SuperLU_DIST can solve
  //gridpack::math::Matrix::StorageType denseType = gridpack::math::Matrix::Dense;
  timer->start(t_matset);
#if 0
  boost::shared_ptr<gridpack::math::Matrix> Y_cDense(gridpack::math::storageType(*Y_c, denseType));
   
  // Form matrix permYmod
  factory.setMode(permYMOD);
  gridpack::mapper::FullMatrixMap<DSNetwork> pymMap(network);
  boost::shared_ptr<gridpack::math::Matrix>  permYmod= pymMap.mapToMatrix();
  ///busIO.header("\n=== permYmod: ============\n");
  ///permYmod->print();

  // Update ybus: ybus = ybus+permYmod (different dimension) => prefy11ybus
#endif
  factory.setMode(updateYbus);
#if 0

  boost::shared_ptr<gridpack::math::Vector> vPermYmod(diagonal(*permYmod));
  ///busIO.header("\n=== vPermYmod: ============\n");
  ///vPermYmod->print();
  gridpack::mapper::BusVectorMap<DSNetwork> permYmodMap(network);
  permYmodMap.mapToBus(vPermYmod);
#endif

  boost::shared_ptr<gridpack::math::Matrix> prefy11ybus = ybusMap.mapToMatrix();
  timer->stop(t_matset);
  ///branchIO.header("\n=== prefy11ybus: ============\n");
  ///prefy11ybus->print();

  // Solve linear equations of ybus * X = Y_c
  //gridpack::math::LinearSolver solver1(*prefy11ybus);
  //solver1.configure(cursor);
  int t_solve = timer->createCategory("Solve Linear Equation");
  timer->start(t_solve);
  gridpack::math::LinearMatrixSolver solver1(*prefy11ybus);
  boost::shared_ptr<gridpack::math::Matrix> prefy11X(solver1.solve(*Y_cDense));
  timer->stop(t_solve);
  ///branchIO.header("\n=== prefy11X: ============\n");
  ///prefy11X->print(); 
  
  //-----------------------------------------------------------------------
  // Compute prefy11
  //-----------------------------------------------------------------------
  // Form reduced admittance matrix prefy11: prefy11 = Y_b * X
  timer->start(t_matmul);
  boost::shared_ptr<gridpack::math::Matrix> prefy11(multiply(*Y_b, *prefy11X)); 
  timer->stop(t_matmul);
  // Update prefy11: prefy11 = Y_a + prefy11
  prefy11->add(*Y_a);
  ///branchIO.header("\n=== Reduced Ybus: prefy11: ============\n");
  ///prefy11->print();

  //-----------------------------------------------------------------------
  // Compute fy11
  // Update ybus values at fault stage
  //-----------------------------------------------------------------------
  // Read the switch info from faults Event from input.xml
  int sw2_2 = faults[0].from_idx - 1;
  int sw3_2 = faults[0].to_idx - 1;

  boost::shared_ptr<gridpack::math::Matrix> fy11ybus(prefy11ybus->clone());
  ///branchIO.header("\n=== fy11ybus(original): ============\n");
  ///fy11ybus->print();
  gridpack::ComplexType x(0.0, -1e7);
#if 0
  fy11ybus->setElement(sw2_2, sw2_2, -x);
  fy11ybus->ready();
#else
  timer->start(t_matset);
  factory.setEvent(faults[0]);
  factory.setMode(onFY);
  ybusMap.overwriteMatrix(fy11ybus);
  timer->stop(t_matset);
#endif
  ///branchIO.header("\n=== fy11ybus: ============\n");
  ///fy11ybus->print();

  // Solve linear equations of fy11ybus * X = Y_c
  timer->start(t_solve);
  gridpack::math::LinearMatrixSolver solver2(*fy11ybus);
  boost::shared_ptr<gridpack::math::Matrix> fy11X(solver2.solve(*Y_cDense)); 
  timer->stop(t_solve);
  ///branchIO.header("\n=== fy11X: ============\n");
  ///fy11X->print();
  
  // Form reduced admittance matrix fy11: fy11 = Y_b * X
  timer->start(t_matmul);
  boost::shared_ptr<gridpack::math::Matrix> fy11(multiply(*Y_b, *fy11X)); 
  timer->stop(t_matmul);
  // Update fy11: fy11 = Y_a + fy11
  fy11->add(*Y_a);
  ///branchIO.header("\n=== Reduced Ybus: fy11: ============\n");
  ///fy11->print();

  //-----------------------------------------------------------------------
  // Compute posfy11
  // Update ybus values at clear fault stage
  //-----------------------------------------------------------------------
  // Get the updating factor for posfy11 stage ybus
#if 0
  gridpack::ComplexType myValue = factory.setFactor(sw2_2, sw3_2);
#endif
  timer->start(t_matset);
  boost::shared_ptr<gridpack::math::Matrix> posfy11ybus(prefy11ybus->clone());
  timer->stop(t_matset);
  ///branchIO.header("\n=== posfy11ybus (original): ============\n");
  ///posfy11ybus->print();
#if 0
  gridpack::ComplexType big(0.0, 1e7);
  gridpack::ComplexType x11 = big - myValue;
  posfy11ybus->addElement(sw2_2, sw2_2, x+x11);

  gridpack::ComplexType x12 = myValue;
  posfy11ybus->addElement(sw2_2, sw3_2, x12);

  gridpack::ComplexType x21 = myValue;
  posfy11ybus->addElement(sw3_2, sw2_2, x21);

  gridpack::ComplexType x22 = -myValue; 
  posfy11ybus->addElement(sw3_2, sw3_2, x22);
  
  posfy11ybus->ready(); 
#else
  timer->start(t_matset);
  factory.setMode(posFY);
  ybusMap.incrementMatrix(posfy11ybus);
  timer->stop(t_matset);
#endif
  ///branchIO.header("\n=== posfy11ybus: ============\n");
  ///posfy11ybus->print();
    
  // Solve linear equations of posfy11ybus * X = Y_c
  timer->start(t_solve);
  gridpack::math::LinearMatrixSolver solver3(*posfy11ybus);
  boost::shared_ptr<gridpack::math::Matrix> posfy11X(solver3.solve(*Y_cDense)); 
  timer->stop(t_solve);
  ///branchIO.header("\n=== posfy11X: ============\n");
  ///posfy11X->print();
  
  // Form reduced admittance matrix posfy11: posfy11 = Y_b * X
  timer->start(t_matmul);
  boost::shared_ptr<gridpack::math::Matrix> posfy11(multiply(*Y_b, *posfy11X)); 
  timer->stop(t_matmul);
  // Update posfy11: posfy11 = Y_a + posfy11
  posfy11->add(*Y_a);
  ///branchIO.header("\n=== Reduced Ybus: posfy11: ============\n");
  ///posfy11->print();

  //-----------------------------------------------------------------------
  // Integration implementation (Modified Euler Method)
  //-----------------------------------------------------------------------
 
  // Map to create vector pelect  
  int t_vecset = timer->createCategory("Setup Vector");
  timer->start(t_vecset);
  factory.setMode(init_pelect);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap1(network);
  boost::shared_ptr<gridpack::math::Vector> pelect = XMap1.mapToVector();
  ///busIO.header("\n=== pelect: ===\n");
  ///pelect->print();

  // Map to create vector eprime_s0 
  factory.setMode(init_eprime);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap2(network);
  boost::shared_ptr<gridpack::math::Vector> eprime_s0 = XMap2.mapToVector();
  ///busIO.header("\n=== eprime: ===\n");
  ///eprime_s0->print();

  // Map to create vector mac_ang_s0
  factory.setMode(init_mac_ang);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap3(network);
  boost::shared_ptr<gridpack::math::Vector> mac_ang_s0 = XMap3.mapToVector();
  ///busIO.header("\n=== mac_ang_s0: ===\n");
  ///mac_ang_s0->print();

  // Map to create vector mac_spd_s0
  factory.setMode(init_mac_spd);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap4(network);
  boost::shared_ptr<gridpack::math::Vector> mac_spd_s0 = XMap4.mapToVector();
  ///busIO.header("\n=== mac_spd_s0: ===\n");
  ///mac_spd_s0->print();

  // Map to create vector eqprime
  factory.setMode(init_eqprime);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap5(network);
  boost::shared_ptr<gridpack::math::Vector> eqprime = XMap5.mapToVector();
  ///busIO.header("\n=== eqprime: ===\n");
  ///eqprime->print();

  // Map to create vector pmech  
  factory.setMode(init_pmech);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap6(network);
  boost::shared_ptr<gridpack::math::Vector> pmech = XMap6.mapToVector();
  ///busIO.header("\n=== pmech: ===\n");
  ///pmech->print();

  // Map to create vector mva
  factory.setMode(init_mva);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap7(network);
  boost::shared_ptr<gridpack::math::Vector> mva = XMap7.mapToVector();
  ///busIO.header("\n=== mva: ===\n");
  ///mva->print();

  // Map to create vector d0
  factory.setMode(init_d0);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap8(network);
  boost::shared_ptr<gridpack::math::Vector> d0 = XMap8.mapToVector();
  ///busIO.header("\n=== d0: ===\n");
  ///d0->print();

  // Map to create vector h
  factory.setMode(init_h);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap9(network);
  boost::shared_ptr<gridpack::math::Vector> h = XMap9.mapToVector();
  ///busIO.header("\n=== h: ===\n");
  ///h->print();

  ///busIO.header("\n============Start of Simulation:=====================\n");

  // Declare vector mac_ang_s1, mac_spd_s1
  boost::shared_ptr<gridpack::math::Vector> mac_ang_s1(mac_ang_s0->clone()); 
  boost::shared_ptr<gridpack::math::Vector> mac_spd_s1(mac_spd_s0->clone()); 
 
  // Declare vector eprime_s1
  boost::shared_ptr<gridpack::math::Vector> eprime_s1(mac_ang_s0->clone());
 
  // Declare vector dmac_ang_s0, dmac_spd_s0, dmac_ang_s1, dmac_spd_s1
  boost::shared_ptr<gridpack::math::Vector> dmac_ang_s0(mac_ang_s0->clone()); 
  boost::shared_ptr<gridpack::math::Vector> dmac_ang_s1(mac_ang_s0->clone()); 
  boost::shared_ptr<gridpack::math::Vector> dmac_spd_s0(mac_spd_s0->clone()); 
  boost::shared_ptr<gridpack::math::Vector> dmac_spd_s1(mac_spd_s0->clone()); 

  // Declare vector curr
  boost::shared_ptr<gridpack::math::Vector> curr(mac_ang_s0->clone());
  timer->stop(t_vecset);

  timer->start(t_trans);
  boost::shared_ptr<gridpack::math::Matrix> trans_prefy11(transpose(*prefy11));
  boost::shared_ptr<gridpack::math::Matrix> trans_fy11(transpose(*fy11));
  boost::shared_ptr<gridpack::math::Matrix> trans_posfy11(transpose(*posfy11));
  timer->stop(t_trans);

  timer->start(t_vecset);
  boost::shared_ptr<gridpack::math::Vector> vecTemp(mac_ang_s0->clone());
  timer->stop(t_vecset);

  // Simulation related variables
  int simu_k;
  int t_step[20];
  double t_width[20];
  int flagF;
  int S_Steps;
  int last_S_Steps;
  int steps3, steps2, steps1;
  double h_sol1, h_sol2;
  int flagF1, flagF2;
  int I_Steps;

  const double sysFreq = 60.0;
  double pi = 4.0*atan(1.0);
  const double basrad = 2.0 * pi * sysFreq;
  gridpack::ComplexType jay(0.0, 1.0);

  // switch info is set up here 
  int nswtch = 4;
  static double sw1[4];
  static double sw7[4];
  sw1[0] = 0.0;
  sw1[1] = faults[0].start;
  sw1[2] = faults[0].end;
  sw1[3] = sim_time;
  sw7[0] = time_step;
  sw7[1] = faults[0].step;
  sw7[2] = time_step;
  sw7[3] = time_step;
  simu_k = 0; 
  for (int i = 0; i < nswtch-1; i++) {
    t_step[i] = (int) ((sw1[i+1] -sw1[i]) / sw7[i]);   
    t_width[i] = (sw1[i+1] - sw1[i]) / t_step[i];
    simu_k += t_step[i];
  }
  simu_k++;

  steps3 = t_step[0] + t_step[1] + t_step[2] - 1;
  steps2 = t_step[0] + t_step[1] - 1;
  steps1 = t_step[0] - 1;
  h_sol1 = t_width[0];
  h_sol2 = h_sol1;
  flagF1 = 0; 
  flagF2 = 0; 
  S_Steps = 1; 
  last_S_Steps = -1;

  for (I_Steps = 0; I_Steps < simu_k+1; I_Steps++) {
    if (I_Steps < steps1) {
      S_Steps = I_Steps;
      flagF1 = 0;
      flagF2 = 0;
    } else if (I_Steps == steps1) { 
      S_Steps = I_Steps;
      flagF1 = 0;
      flagF2 = 1;
    } else if (I_Steps == steps1+1) {
      S_Steps = I_Steps;
      flagF1 = 1;
      flagF2 = 1;
    } else if ((I_Steps>steps1+1) && (I_Steps<steps2+1)) {
      S_Steps = I_Steps - 1;
      flagF1 = 1;
      flagF2 = 1;
    } else if (I_Steps==steps2+1) {
      S_Steps = I_Steps - 1;
      flagF1 = 1;
      flagF2 = 2;
    } else if (I_Steps==steps2+2) {
      S_Steps = I_Steps - 1;
      flagF1 = 2;
      flagF2 = 2;
    } else if (I_Steps>steps2+2) {
      S_Steps = I_Steps - 2;
      flagF1 = 2;
      flagF2 = 2;
    }

    if (I_Steps !=0 && last_S_Steps != S_Steps) {
      mac_ang_s0->equate(*mac_ang_s1); 
      mac_spd_s0->equate(*mac_spd_s1); 
      eprime_s0->equate(*eprime_s1); 
    }

    vecTemp->equate(*mac_ang_s0); 
    vecTemp->scale(jay); 
    vecTemp->exp(); 
    vecTemp->elementMultiply(*eqprime); 
    eprime_s0->equate(*vecTemp); 
     
    // ---------- CALL i_simu_innerloop(k,S_Steps,flagF1): ----------
    int t_trnsmul = timer->createCategory("Transpose Multiply");
    timer->start(t_trnsmul);
    if (flagF1 == 0) {
      curr.reset(multiply(*trans_prefy11, *eprime_s0)); //MatMultTranspose(prefy11, eprime_s0, curr);
      //transposeMultiply(*prefy11,*eprime_s0,*curr);
    } else if (flagF1 == 1) {
      curr.reset(multiply(*trans_fy11, *eprime_s0)); //MatMultTranspose(fy11, eprime_s0, curr);
      //transposeMultiply(*fy11,*eprime_s0,*curr);
    } else if (flagF1 == 2) {
      curr.reset(multiply(*trans_posfy11, *eprime_s0)); //MatMultTranspose(posfy11, eprime_s0, curr);
      //transposeMultiply(*posfy11,*eprime_s0,*curr);
    } 
    timer->stop(t_trnsmul);
    
    // ---------- CALL mac_em2(k,S_Steps): ----------
    // ---------- pelect: ----------
    curr->conjugate(); 
    vecTemp->equate(*eprime_s0);
    vecTemp->elementMultiply(*curr);
    pelect->equate(*vecTemp); 
    pelect->real(); //Get the real part of pelect
    // ---------- dmac_ang: ----------
    vecTemp->equate(*mac_spd_s0); 
    vecTemp->add(-1.0); 
    vecTemp->scale(basrad); 
    dmac_ang_s0->equate(*vecTemp); 
    // ---------- dmac_spd: ----------
    vecTemp->equate(*pelect);
    vecTemp->elementMultiply(*mva); 
    dmac_spd_s0->equate(*pmech); 
    dmac_spd_s0->add(*vecTemp, -1.0); 
    vecTemp->equate(*mac_spd_s0); 
    vecTemp->add(-1.0); 
    vecTemp->elementMultiply(*d0); 
    dmac_spd_s0->add(*vecTemp, -1.0); 
    vecTemp->equate(*h); 
    vecTemp->scale(2.0); 
    dmac_spd_s0->elementDivide(*vecTemp); 

    mac_ang_s1->equate(*mac_ang_s0); 
    vecTemp->equate(*dmac_ang_s0);
    mac_ang_s1->add(*dmac_ang_s0, h_sol1); 
    mac_spd_s1->equate(*mac_spd_s0); 
    vecTemp->equate(*dmac_spd_s0);
    mac_spd_s1->add(*vecTemp, h_sol1); 

    vecTemp->equate(*mac_ang_s1); 
    vecTemp->scale(jay); 
    vecTemp->exp(); 
    vecTemp->elementMultiply(*eqprime);
    eprime_s1->equate(*vecTemp); 

    // ---------- CALL i_simu_innerloop2(k,S_Steps+1,flagF2): ----------
    timer->start(t_trnsmul);
    if (flagF2 == 0) {
      curr.reset(multiply(*trans_prefy11, *eprime_s1)); //MatMultTranspose(prefy11, eprime_s1, curr);
      //transposeMultiply(*prefy11,*eprime_s1,*curr);
    } else if (flagF2 == 1) {
      curr.reset(multiply(*trans_fy11, *eprime_s1)); //MatMultTranspose(fy11, eprime_s1, curr);
      //transposeMultiply(*fy11,*eprime_s1,*curr);
    } else if (flagF2 == 2) {
      curr.reset(multiply(*trans_posfy11, *eprime_s1)); //MatMultTranspose(posfy11, eprime_s1, curr);
      //transposeMultiply(*posfy11,*eprime_s1,*curr);
    }
    timer->stop(t_trnsmul);

    // ---------- CALL mac_em2(k,S_Steps+1): ---------- 
    // ---------- pelect: ----------
    curr->conjugate(); 
    vecTemp->equate(*eprime_s1);
    vecTemp->elementMultiply(*curr); 
    pelect->equate(*vecTemp); 
    pelect->real(); //Get the real part of pelect
    // ---------- dmac_ang: ----------
    vecTemp->equate(*mac_spd_s1); 
    vecTemp->add(-1.0); 
    vecTemp->scale(basrad); 
    dmac_ang_s1->equate(*vecTemp); 
    // ---------- dmac_spd: ----------
    vecTemp->equate(*pelect);
    vecTemp->elementMultiply(*mva); 
    dmac_spd_s1->equate(*pmech); 
    dmac_spd_s1->add(*vecTemp, -1.0); 
    vecTemp->equate(*mac_spd_s1); 
    vecTemp->add(-1.0); 
    vecTemp->elementMultiply(*d0); 
    dmac_spd_s1->add(*vecTemp, -1.0); 
    vecTemp->equate(*h); 
    vecTemp->scale(2.0); 
    dmac_spd_s1->elementDivide(*vecTemp); 

    vecTemp->equate(*dmac_ang_s0); 
    vecTemp->add(*dmac_ang_s1); 
    vecTemp->scale(0.5); 
    mac_ang_s1->equate(*mac_ang_s0); 
    mac_ang_s1->add(*vecTemp, h_sol2); 
    vecTemp->equate(*dmac_spd_s0); 
    vecTemp->add(*dmac_spd_s1); 
    vecTemp->scale(0.5); 
    mac_spd_s1->equate(*mac_spd_s0); 
    mac_spd_s1->add(*vecTemp, h_sol2); 

    // Print to screen
    if (last_S_Steps != S_Steps) {
      //sprintf(ioBuf, "\n========================S_Steps = %d=========================\n", S_Steps);
      //busIO.header(ioBuf);
      //mac_ang_s0->print();  
      //mac_spd_s0->print();  
      //pmech->print();
      //pelect->print();
      //sprintf(ioBuf, "========================End of S_Steps = %d=========================\n\n", S_Steps);
      //busIO.header(ioBuf);
    }
    if (I_Steps == simu_k) {
      factory.setMode(init_mac_ang);
      XMap3.mapToBus(mac_ang_s1);
      factory.setMode(init_mac_spd);
      XMap3.mapToBus(mac_spd_s1);
      factory.setMode(init_pmech);
      XMap6.mapToBus(pmech);
      factory.setMode(init_pelect);
      XMap1.mapToBus(pelect);
      sprintf(ioBuf, "\n========================S_Steps = %d=========================\n", S_Steps+1);
      busIO.header(ioBuf);
      sprintf(ioBuf, "\n         Bus ID     Generator ID"
          "    mac_ang         mac_spd         mech            elect\n\n");
      busIO.header(ioBuf);
      //mac_ang_s1->print();  
      //mac_spd_s1->print();  
      //pmech->print();
      //pelect->print();
      busIO.write();
      sprintf(ioBuf, "\n========================End of S_Steps = %d=========================\n\n", S_Steps+1);
      busIO.header(ioBuf);
    } // End of Print to screen 
    
    last_S_Steps = S_Steps;
  }
  timer->stop(t_total);
  timer->dump();
}

/**
 * Utility function to convert faults that are in event list into
 * internal data structure that can be used by code
 * @param cursors list of cursors pointing to individual events in input
 * deck
 * @return list of event data structures
 */
std::vector<gridpack::dynamic_simulation::DSBranch::Event>
   gridpack::dynamic_simulation::DSApp::setFaultEvents(
   std::vector<gridpack::utility::Configuration::CursorPtr > events)
{
  int size = events.size();
  int idx;
  std::vector<gridpack::dynamic_simulation::DSBranch::Event> faults;
  // Parse fault events
  for (idx=0; idx<size; idx++) {
    gridpack::dynamic_simulation::DSBranch::Event event;
    event.start = events[idx]->get("beginFault",0.0);
    event.end = events[idx]->get("endFault",0.0);
    std::string indices = events[idx]->get("faultBranch","0 0");
    //Parse indices to get from and to indices of branch
    int ntok1 = indices.find_first_not_of(' ',0);
    int ntok2 = indices.find(' ',ntok1);
    if (ntok2 - ntok1 > 0 && ntok1 != std::string::npos && ntok2 !=
        std::string::npos) {
      event.from_idx = atoi(indices.substr(ntok1,ntok2-ntok1).c_str());
      ntok1 = indices.find_first_not_of(' ',ntok2);
      ntok2 = indices.find(' ',ntok1);
      if (ntok1 != std::string::npos && ntok1 < indices.length()) {
        if (ntok2 == std::string::npos) {
          ntok2 = indices.length();
        }
        event.to_idx = atoi(indices.substr(ntok1,ntok2-ntok1).c_str());
      } else {
        event.from_idx = 0;
        event.to_idx = 0;
      }
    } else {
      event.from_idx = 0;
      event.to_idx = 0;
    }
    event.step = events[idx]->get("timeStep",0.0);
    if (event.step != 0.0 && event.end != 0.0 && event.from_idx != event.to_idx) {
      faults.push_back(event);
    }
  }
#if 0
  // Find local indices of branches on which faults occur. Start by constructing
  // a map object that maps to local branch indices using branch index pairs as the key
  std::map<std::pair<int, int>, int> pairMap;
  int numBranch = network->numBranches();
  for (idx = 0; idx<numBranch; idx++) {
    // Only set branch index for locally held branches
    if (network->getActiveBranch(idx)) {
      int idx1, idx2;
      network->getOriginalBranchEndpoints(idx, &idx1, &idx2);
      std::pair<int, int> branch_pair(idx1, idx2);
      pairMap.insert(std::pair<std::pair<int, int>, int>(branch_pair,idx));
    }
  }
  // run through all events and see if the branch exists on this processor. If
  // it does, then set the branch_idx member to the local branch index.
  size = faults.size();
  for (idx=0; idx<size; idx++) {
    std::map<std::pair<int, int>, int>::iterator it;
    std::pair<int, int> branch_pair(faults[idx].from_idx, faults[idx].to_idx);
    it = pairMap.find(branch_pair);
    if (it != pairMap.end()) {
      faults[idx].branch_idx = it->second;
    }
  }
#endif
  return faults;
}
