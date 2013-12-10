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
void gridpack::dynamic_simulation::DSApp::execute(void)
{
  gridpack::parallel::Communicator world;
  boost::shared_ptr<DSNetwork> network(new DSNetwork(world));

  // read configuration file
  gridpack::utility::Configuration *config = gridpack::utility::Configuration::configuration();
  config->open("input.xml",world);
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Dynamic_simulation");
  std::string filename = cursor->get("networkConfiguration",
      "No network configuration specified");

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
  // DAE implementation 
  //-----------------------------------------------------------------------
  // Set initial conditions
  factory.setMode(DAE_init);  
  gridpack::mapper::BusVectorMap<DSNetwork> XMap(network);
  boost::shared_ptr<gridpack::math::Vector> X = XMap.mapToVector();
  printf("\n=== DAE initial X: ============\n");
  X->print();

  //gridpack::dynamic_simulation::DSSimu simu(network);
  //simu.setIFunction();
  
}

