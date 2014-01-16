/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_app.cpp
 * @author Shuangshuang.Jin
 * @date   September 19, 2013
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
#include "gridpack/applications/dynamic_simulation/ds_app.hpp"
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/applications/dynamic_simulation/ds_factory.hpp"

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
  gridpack::parallel::Communicator world;
  boost::shared_ptr<DSNetwork> network(new DSNetwork(world));

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

  // Read in information about fault events and store them in internal data
  // structure
  cursor = config->getCursor("Configuration.Dynamic_simulation.faultEvents");
  gridpack::utility::Configuration::ChildCursors events;
  if (cursor) cursor->children(events);
  std::vector<Event> faults = setFaultEvents(events); 

  // load input file
  gridpack::parser::PTI23_parser<DSNetwork> parser(network);
  parser.parse(filename.c_str());

  // partition network
  network->partition();

  // create factory
  gridpack::dynamic_simulation::DSFactory factory(network);
  factory.load();

  // set network components using factory
  factory.setComponents();

  // set YBus components so that you can create Y matrix  
  factory.setYBus();

  factory.setMode(YBUS);
  gridpack::mapper::FullMatrixMap<DSNetwork> ybusMap(network);
  boost::shared_ptr<gridpack::math::Matrix> ybus = ybusMap.mapToMatrix();
  printf("\n=== ybus: ============\n");
  ybus->print();

  // Form constant impedance load admittance yl for all buses and add it to 
  // system Y matrix: ybus = ybus + yl
  factory.setMode(YL); 
  gridpack::mapper::FullMatrixMap<DSNetwork> newybusMap(network);
  ybus = newybusMap.mapToMatrix();
  printf("\n=== ybus after added yl: ============\n");
  ybus->print();
 
  // Construct permutation matrix perm by checking for multiple generators at a bus
  factory.setMode(PERM);
  gridpack::mapper::FullMatrixMap<DSNetwork> permMap(network);
  boost::shared_ptr<gridpack::math::Matrix> perm = permMap.mapToMatrix();
  printf("\n=== perm: ============\n");
  perm->print(); 

  // Form a transposed matrix of perm
  boost::shared_ptr<gridpack::math::Matrix> permTrans(transpose(*perm));
  printf("\n=== permTrans: ============\n");
  permTrans->print();

  // Construct matrix Y_a using extracted xd and ra from gen data, 
  // and construct its diagonal matrix diagY_a
  factory.setMode(YA);
  gridpack::mapper::FullMatrixMap<DSNetwork> yaMap(network);
  //boost::shared_ptr<gridpack::math::Matrix> Y_a = yaMap.mapToMatrix();
  //printf("\n=== Y_a: ============\n");
  //Y_a->print(); 
  boost::shared_ptr<gridpack::math::Matrix> diagY_a = yaMap.mapToMatrix();
  printf("\n=== diagY_a: ============\n");
  diagY_a->print(); 
  // Convert diagY_a from sparse matrix to dense matrix Y_a so that SuperLU_DIST can solve
  gridpack::math::Matrix::StorageType denseType = gridpack::math::Matrix::Dense;
  boost::shared_ptr<gridpack::math::Matrix> Y_a(gridpack::math::storageType(*diagY_a, denseType));

  // Construct matrix Ymod: Ymod = diagY_a * permTrans
  boost::shared_ptr<gridpack::math::Matrix> Ymod(multiply(*diagY_a, *permTrans));
  printf("\n=== Ymod: ============\n");
  Ymod->print(); 
 
  // Form matrix Y_b: Y_b(1:ngen, jg) = -Ymod, where jg represents the 
  // corresponding index sets of buses that the generators are connected to. 
  // Then construct Y_b's transposed matrix Y_c: Y_c = Y_b'
  // This two steps can be done by using a P matrix to get Y_c directly.
  factory.setMode(PMatrix);
  gridpack::mapper::FullMatrixMap<DSNetwork> pMap(network);
  boost::shared_ptr<gridpack::math::Matrix> P = pMap.mapToMatrix();
  printf("\n=== P: ============\n");
  P->print();
  boost::shared_ptr<gridpack::math::Matrix> Y_c(multiply(*P, *Ymod)); 
  printf("\n=== Y_c: ============\n");
  Y_c->print();
  boost::shared_ptr<gridpack::math::Matrix> Y_b(transpose(*Y_c));
  printf("\n=== Y_b: ============\n");
  Y_b->print();
  Y_c->scale(-1.0); // scale Y_c by -1.0 for linear solving
  // Convert Y_c from sparse matrix to dense matrix Y_cDense so that SuperLU_DIST can solve
  //gridpack::math::Matrix::StorageType denseType = gridpack::math::Matrix::Dense;
  boost::shared_ptr<gridpack::math::Matrix> Y_cDense(gridpack::math::storageType(*Y_c, denseType));
   
  // Form matrix permYmod
  boost::shared_ptr<gridpack::math::Matrix> permYmod(multiply(*perm, *Ymod)); 
  printf("\n=== permYmod: ============\n");
  permYmod->print();

  // Update ybus: ybus = ybus+permYmod (different dimension) => prefy11ybus
  factory.setMode(updateYbus);

  printf("\n=== vPermYmod: ============\n");
  boost::shared_ptr<gridpack::math::Vector> vPermYmod(diagonal(*permYmod));
  vPermYmod->print();
  gridpack::mapper::BusVectorMap<DSNetwork> permYmodMap(network);
  permYmodMap.mapToBus(vPermYmod);

  boost::shared_ptr<gridpack::math::Matrix> prefy11ybus = ybusMap.mapToMatrix();
  printf("\n=== prefy11ybus: ============\n");
  prefy11ybus->print();

  // Solve linear equations of ybus * X = Y_c
  //gridpack::math::LinearSolver solver1(*prefy11ybus);
  //solver1.configure(cursor);
  gridpack::math::LinearMatrixSolver solver1(*prefy11ybus);
  boost::shared_ptr<gridpack::math::Matrix> prefy11X(solver1.solve(*Y_cDense));
  printf("\n=== prefy11X: ============\n");
  prefy11X->print(); 
  
  //-----------------------------------------------------------------------
  // Compute prefy11
  //-----------------------------------------------------------------------
  // Form reduced admittance matrix prefy11: prefy11 = Y_b * X
  boost::shared_ptr<gridpack::math::Matrix> prefy11(multiply(*Y_b, *prefy11X)); 
  // Update prefy11: prefy11 = Y_a + prefy11
  prefy11->add(*Y_a);
  printf("\n=== Reduced Ybus: prefy11: ============\n");
  prefy11->print();

  //-----------------------------------------------------------------------
  // Compute fy11
  // Update ybus values at fault stage
  //-----------------------------------------------------------------------
  // assume switch info is set up here instead of reading from the input file
  int sw2_2 = 5; // 6-1
  int sw3_2 = 6; // 7-1

  //factory.setMode(FY); 
  //gridpack::mapper::FullMatrixMap<DSNetwork> fy11ybusMap(network);
  //boost::shared_ptr<gridpack::math::Matrix> fy11ybus = fy11ybusMap.mapToMatrix();
  //fy11ybus->print();

  boost::shared_ptr<gridpack::math::Matrix> fy11ybus(prefy11ybus->clone());
  gridpack::ComplexType x(0.0, -1e7);
  fy11ybus->setElement(sw2_2, sw2_2, -x);
  fy11ybus->ready();
  printf("\n=== fy11ybus: ============\n");
  fy11ybus->print();

  // Solve linear equations of fy11ybus * X = Y_c
  //gridpack::math::LinearSolver solver2(*fy11ybus);
  gridpack::math::LinearMatrixSolver solver2(*fy11ybus);
  boost::shared_ptr<gridpack::math::Matrix> fy11X(solver2.solve(*Y_cDense)); 
  printf("\n=== fy11X: ============\n");
  fy11X->print();
  
  // Form reduced admittance matrix fy11: fy11 = Y_b * X
  boost::shared_ptr<gridpack::math::Matrix> fy11(multiply(*Y_b, *fy11X)); 
  // Update fy11: fy11 = Y_a + fy11
  fy11->add(*Y_a);
  printf("\n=== Reduced Ybus: fy11: ============\n");
  fy11->print();

  //-----------------------------------------------------------------------
  // Compute posfy11
  // Update ybus values at clear fault stage
  //-----------------------------------------------------------------------
  //factory.setMode(POSFY); 
  //gridpack::mapper::FullMatrixMap<DSNetwork> posfy11ybusMap(network);
  //boost::shared_ptr<gridpack::math::Matrix> posfy11ybus = posfy11ybusMap.mapToMatrix();
  //posfy11ybus->print();

  // Get the updating factor for posfy11 stage ybus
  gridpack::ComplexType myValue = factory.setFactor(sw2_2, sw3_2);
  boost::shared_ptr<gridpack::math::Matrix> posfy11ybus(prefy11ybus->clone());

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
  printf("\n=== posfy11ybus: ============\n");
  posfy11ybus->print();
    
  // Solve linear equations of posfy11ybus * X = Y_c
  //gridpack::math::LinearSolver solver3(*posfy11ybus);
  gridpack::math::LinearMatrixSolver solver3(*posfy11ybus);
  boost::shared_ptr<gridpack::math::Matrix> posfy11X(solver3.solve(*Y_cDense)); 
  printf("\n=== posfy11X: ============\n");
  posfy11X->print();
  
  // Form reduced admittance matrix posfy11: posfy11 = Y_b * X
  boost::shared_ptr<gridpack::math::Matrix> posfy11(multiply(*Y_b, *posfy11X)); 
  //// Update posfy11: posfy11 = Y_a + posfy11
  posfy11->add(*Y_a);
  printf("\n=== Reduced Ybus: posfy11: ============\n");
  posfy11->print();

  //-----------------------------------------------------------------------
  // Integration implementation (Modified Euler Method)
  //-----------------------------------------------------------------------
  // Set initial conditions
  //factory.setMode(DAE_init);
 
  // Map to create vector pelect  
  factory.setMode(init_pelect);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap1(network);
  boost::shared_ptr<gridpack::math::Vector> pelect = XMap1.mapToVector();
  printf("\n=== pelect: ===\n");
  pelect->print();

  // Map to create vector eprime_s0 
  factory.setMode(init_eprime);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap2(network);
  boost::shared_ptr<gridpack::math::Vector> eprime_s0 = XMap2.mapToVector();
  printf("\n=== eprime: ===\n");
  eprime_s0->print();

  // Map to create vector mac_ang_s0
  factory.setMode(init_mac_ang);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap3(network);
  boost::shared_ptr<gridpack::math::Vector> mac_ang_s0 = XMap3.mapToVector();
  printf("\n=== mac_ang_s0: ===\n");
  mac_ang_s0->print();

  // Map to create vector mac_spd_s0
  factory.setMode(init_mac_spd);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap4(network);
  boost::shared_ptr<gridpack::math::Vector> mac_spd_s0 = XMap4.mapToVector();
  printf("\n=== mac_spd_s0: ===\n");
  mac_spd_s0->print();

  // Map to create vector eqprime
  factory.setMode(init_eqprime);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap5(network);
  boost::shared_ptr<gridpack::math::Vector> eqprime = XMap5.mapToVector();
  printf("\n=== eqprime: ===\n");
  eqprime->print();

  // Map to create vector pmech  
  factory.setMode(init_pmech);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap6(network);
  boost::shared_ptr<gridpack::math::Vector> pmech = XMap6.mapToVector();
  printf("\n=== pmech: ===\n");
  pmech->print();

  // Map to create vector mva
  factory.setMode(init_mva);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap7(network);
  boost::shared_ptr<gridpack::math::Vector> mva = XMap7.mapToVector();
  printf("\n=== mva: ===\n");
  mva->print();

  // Map to create vector d0
  factory.setMode(init_d0);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap8(network);
  boost::shared_ptr<gridpack::math::Vector> d0 = XMap8.mapToVector();
  printf("\n=== d0: ===\n");
  d0->print();

  // Map to create vector h
  factory.setMode(init_h);
  gridpack::mapper::BusVectorMap<DSNetwork> XMap9(network);
  boost::shared_ptr<gridpack::math::Vector> h = XMap9.mapToVector();
  printf("\n=== h: ===\n");
  h->print();
  printf("\n============Start of Simulation:=====================\n");

  // Declare vector mac_ang_s1, mac_spd_s1
  //boost::shared_ptr<gridpack::math::Vector> mac_ang_s1; 
  //boost::shared_ptr<gridpack::math::Vector> mac_spd_s1;
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
  boost::shared_ptr<gridpack::math::Vector> curr;

  boost::shared_ptr<gridpack::math::Matrix> trans_prefy11(transpose(*prefy11));
  boost::shared_ptr<gridpack::math::Matrix> trans_fy11(transpose(*fy11));
  boost::shared_ptr<gridpack::math::Matrix> trans_posfy11(transpose(*posfy11));

  boost::shared_ptr<gridpack::math::Vector> vecTemp(mac_ang_s0->clone());

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

  // assume switch info is set up here instead of reading from the input file
  int nswtch = 4;
  static double sw1[4] = {0.0, 0.03, 0.06, 0.1};
  static double sw7[4] = {0.01, 0.01, 0.01, 0.01}; 
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
      mac_ang_s0->equate(*mac_ang_s1); //VecCopy(mac_ang_s1, mac_ang_s0);
      mac_spd_s0->equate(*mac_spd_s1); //VecCopy(mac_spd_s1, mac_spd_s0);
      eprime_s0->equate(*eprime_s1); //VecCopy(eprime_s1, eprime_s0);
    }

    vecTemp->equate(*mac_ang_s0); //VecCopy(mac_ang_s0, vecTemp);
    vecTemp->scale(jay); //VecScale(vecTemp, PETSC_i);
    vecTemp->exp(); //VecExp(vecTemp);   
    vecTemp->elementMultiply(*eqprime); 
    eprime_s0->equate(*vecTemp); //VecPointwiseMult(eprime_s0, eqprime, vecTemp);
     
    // ---------- CALL i_simu_innerloop(k,S_Steps,flagF1): ----------
    if (flagF1 == 0) {
      curr.reset(multiply(*trans_prefy11, *eprime_s0)); //MatMultTranspose(prefy11, eprime_s0, curr);
    } else if (flagF1 == 1) {
      curr.reset(multiply(*trans_fy11, *eprime_s0)); //MatMultTranspose(fy11, eprime_s0, curr);
    } else if (flagF1 == 2) {
      curr.reset(multiply(*trans_posfy11, *eprime_s0)); //MatMultTranspose(posfy11, eprime_s0, curr);
    } 
    
    // ---------- CALL mac_em2(k,S_Steps): ----------
    // ---------- pelect: ----------
    curr->conjugate(); //VecConjugate(curr);
    vecTemp->equate(*eprime_s0);
    vecTemp->elementMultiply(*curr);
    pelect->equate(*vecTemp); //VecPointwiseMult(pelect, eprime_s0, curr);
    pelect->real(); //Get the real part of pelect
    // ---------- dmac_ang: ----------
    vecTemp->equate(*mac_spd_s0); //VecCopy(mac_spd_s0, vecTemp);
    vecTemp->add(-1.0); //VecShift(vecTemp, -1.0);
    vecTemp->scale(basrad); //VecScale(vecTemp, basrad);
    dmac_ang_s0->equate(*vecTemp); //VecCopy(vecTemp, dmac_ang_s0);
    // ---------- dmac_spd: ----------
    vecTemp->equate(*pelect);
    vecTemp->elementMultiply(*mva); //VecPointwiseMult(vecTemp, pelect, mva); // pelect * gen_mva
    dmac_spd_s0->equate(*pmech); //VecCopy(pmech, dmac_spd_s0);
    dmac_spd_s0->add(*vecTemp, -1.0); //VecAXPY(dmac_spd_s0, -1.0, vecTemp); // pmech - pelect * gen_mva
    vecTemp->equate(*mac_spd_s0); //VecCopy(mac_spd_s0, vecTemp);
    vecTemp->add(-1.0); //VecShift(vecTemp, -1.0);
    vecTemp->elementMultiply(*d0); //VecPointwiseMult(vecTemp1, d0, vecTemp); // gen_d0 * (mac_spd_s0 - 1.0)
    dmac_spd_s0->add(*vecTemp, -1.0); //VecAXPY(dmac_spd_s0, -1.0, vecTemp); // pmech - pelect * gen_mva - gen_d0 * (mac_spd_s0 - 1.0)
    vecTemp->equate(*h); //VecCopy(h, vecTemp);
    vecTemp->scale(2.0); //VecScale(vecTemp, 2); // 2 * gen_h
    dmac_spd_s0->elementDivide(*vecTemp); //VecPointwiseDivide(dmac_spd_s0, dmac_spd_s0, vecTemp); // (pmech-pelect*gen_mva-gen_d0*(mac_spd_s0-1.0) )/(2*gen_h)  

    mac_ang_s1->equate(*mac_ang_s0); //VecCopy(mac_ang_s0, mac_ang_s1);
    vecTemp->equate(*dmac_ang_s0);
    mac_ang_s1->add(*dmac_ang_s0, h_sol1); //VecAXPY(mac_ang_s1, h_sol1, dmac_ang_s0);
    mac_spd_s1->equate(*mac_spd_s0); //VecCopy(mac_spd_s0, mac_spd_s1);
    vecTemp->equate(*dmac_spd_s0);
    mac_spd_s1->add(*vecTemp, h_sol1); //VecAXPY(mac_spd_s1, h_sol1, dmac_spd_s0);

    vecTemp->equate(*mac_ang_s1); //VecCopy(mac_ang_s1, vecTemp);
    vecTemp->scale(jay); //VecScale(vecTemp, PETSC_i);
    vecTemp->exp(); //VecExp(vecTemp);
    vecTemp->elementMultiply(*eqprime);
    eprime_s1->equate(*vecTemp); //VecPointwiseMult(eprime_s1, eqprime, vecTemp);

    // ---------- CALL i_simu_innerloop2(k,S_Steps+1,flagF2): ----------
    if (flagF2 == 0) {
      curr.reset(multiply(*trans_prefy11, *eprime_s1)); //MatMultTranspose(prefy11, eprime_s1, curr)
    } else if (flagF2 == 1) {
      curr.reset(multiply(*trans_fy11, *eprime_s1)); //MatMultTranspose(fy11, eprime_s1, curr);
    } else if (flagF2 == 2) {
      curr.reset(multiply(*trans_posfy11, *eprime_s1)); //MatMultTranspose(posfy11, eprime_s1, curr);
    }

    // ---------- CALL mac_em2(k,S_Steps+1): ---------- 
    // ---------- pelect: ----------
    curr->conjugate(); //VecConjugate(curr);
    vecTemp->equate(*eprime_s1);
    vecTemp->elementMultiply(*curr); 
    pelect->equate(*vecTemp); //VecPointwiseMult(pelect, eprime_s1, curr);
    curr->print();
    eprime_s1->print();
    pelect->print();
    pelect->real(); //Get the real part of pelect
    // ---------- dmac_ang: ----------
    vecTemp->equate(*mac_spd_s1); //VecCopy(mac_spd_s1, vecTemp);
    vecTemp->add(-1.0); //VecShift(vecTemp, -1.0);
    vecTemp->scale(basrad); //VecScale(vecTemp, basrad);
    dmac_ang_s1->equate(*vecTemp); //VecCopy(vecTemp, dmac_ang_s1);
    // ---------- dmac_spd: ----------
    vecTemp->equate(*pelect);
    vecTemp->elementMultiply(*mva); //VecPointwiseMult(vecTemp, pelect, mva); // pelect * gen_mva
    dmac_spd_s1->equate(*pmech); //VecCopy(pmech, dmac_spd_s1);
    dmac_spd_s1->add(*vecTemp, -1.0); //VecAXPY(dmac_spd_s1, -1.0, vecTemp); // pmech - pelect * gen_mva
    vecTemp->equate(*mac_spd_s1); //VecCopy(mac_spd_s1, vecTemp);
    vecTemp->add(-1.0); //VecShift(vecTemp, -1.0);
    vecTemp->elementMultiply(*d0); //VecPointwiseMult(vecTemp1, d0, vecTemp); // gen_d0 * (mac_spd_s1 - 1.0)
    dmac_spd_s1->add(*vecTemp, -1.0); //VecAXPY(dmac_spd_s1, -1.0, vecTemp); // pmech - pelect * gen_mva - gen_d0 * (mac_spd_s1 - 1.0)
    vecTemp->equate(*h); //VecCopy(h, vecTemp);
    vecTemp->scale(2.0); //VecScale(vecTemp, 2); // 2 * gen_h
    dmac_spd_s1->elementDivide(*vecTemp); //VecPointwiseDivide(dmac_spd_s1, dmac_spd_s1, vecTemp); // (pmech-pelect*gen_mva-gen_d0*(mac_spd_s1-1.0) )/(2*gen_h) 

    vecTemp->equate(*dmac_ang_s0); //VecCopy(dmac_ang_s0, vecTemp);
    vecTemp->add(*dmac_ang_s1); //VecAXPY(vecTemp, 1.0, dmac_ang_s1);
    vecTemp->scale(0.5); //VecScale(vecTemp, 0.5);
    mac_ang_s1->equate(*mac_ang_s0); //VecCopy(mac_ang_s0, mac_ang_s1);
    mac_ang_s1->add(*vecTemp, h_sol2); //VecAXPY(mac_ang_s1, h_sol2, vecTemp);
    vecTemp->equate(*dmac_spd_s0); //VecCopy(dmac_spd_s0, vecTemp);
    vecTemp->add(*dmac_spd_s1); //VecAXPY(vecTemp, 1.0, dmac_spd_s1);
    vecTemp->scale(0.5); //VecScale(vecTemp, 0.5);
    mac_spd_s1->equate(*mac_spd_s0); //VecCopy(mac_spd_s0, mac_spd_s1);
    mac_spd_s1->add(*vecTemp, h_sol2); //VecAXPY(mac_spd_s1, h_sol2, vecTemp);

    // Print to screen
    if (last_S_Steps != S_Steps) {
      printf("\n========================S_Steps = %d=========================\n", S_Steps);
      mac_ang_s0->print();  
      mac_spd_s0->print();  
      pmech->print();
      pelect->print();
      printf("========================End of S_Steps = %d=========================\n\n", S_Steps);
    }
    if (I_Steps == simu_k) {
      printf("\n========================S_Steps = %d=========================\n", S_Steps+1);
      mac_ang_s1->print();  
      mac_spd_s1->print();  
      pmech->print();
      pelect->print();
      printf("========================End of S_Steps = %d=========================\n\n", S_Steps+1);
    } // End of Print to screen
    
    last_S_Steps = S_Steps;
  } 
}

/**
 * Utility function to convert faults that are in event list into
 * internal data structure that can be used by code
 * @param cursors list of cursors pointing to individual events in input
 * deck
 * @return list of event data structures
 */
std::vector<gridpack::dynamic_simulation::Event>
   gridpack::dynamic_simulation::DSApp::setFaultEvents(
   std::vector<gridpack::utility::Configuration::CursorPtr > events)
{
  int size = events.size();
  int idx;
  std::vector<gridpack::dynamic_simulation::Event> faults;
  for (idx=0; idx<size; idx++) {
    Event event;
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
    if (event.step != 0.0 && event.end != 0.0 && event.from_idx != event.to_idx)
    {
      faults.push_back(event);
    }
  }
  return faults;
}
