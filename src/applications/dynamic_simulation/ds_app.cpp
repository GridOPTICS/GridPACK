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
  gridpack::utility::Configuration::Cursor *cursor;
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
  boost::shared_ptr<gridpack::math::Matrix> Y_a = yaMap.mapToMatrix();
  printf("\n=== Y_a: ============\n");
  Y_a->print(); 

  /*boost::shared_ptr<gridpack::math::Matrix> diagY_a;
  diagY_a->gridpack::math::Matrix::identity();
  diagY_a->print();

  boost::shared_ptr<gridpack::math::Vector> vY_a(diagonal(*Y_a)); 
  vY_a->print();

  //diagY_a->gridpack::math::Matrix::multiplyDiagonal(*vY_a); 
  //diagY_a->muitiplyDiagonal(*vY_a); 
  printf("\n=== diagY_a: ============\n");
  diagY_a->print(); */

  boost::shared_ptr<gridpack::math::Matrix> diagY_a = Y_a;
  printf("\n=== diagY_a: ============\n");
  diagY_a->print(); 

  /************************************************************************
  !!! Compilation succeeded till this point                                 
  !!! Linear solver are not available in math yet 
  ************************************************************************/

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
  boost::shared_ptr<gridpack::math::Matrix> Y_c(multiply(*P, *Ymod)); 
  printf("\n=== Y_c: ============\n");
  Y_c->print();
   
  // Form matrix permYmod
  boost::shared_ptr<gridpack::math::Matrix> permYmod(multiply(*perm, *Ymod)); 
  printf("\n=== permYmod: ============\n");
  permYmod->print();

  // Update ybus: ybus = ybus+permYmod (different dimension)
  //ybus = ybus + permYmod;
  printf("\n=== ybus after added permYmod: ============\n");
  ybus->print();
/*
  // Solve linear equations of ybus * X = Y_c
  boost::shared_ptr<gridpack::math::Matrix> X; 

  //-----------------------------------------------------------------------
  // Compute prefy11
  //-----------------------------------------------------------------------
  // Form reduced admittance matrix prefy11: prefy11 = Y_b * X
  boost::shared_ptr<gridpack::math::Matrix> prefy11(multiply(*Y_b, *X)); 
  //prefy11->print();

  // Update prefy11: prefy11 = Y_a + prefy11
  prefy11->add(*Y_a);
*/
  //-----------------------------------------------------------------------
  // Compute fy11
  // Update ybus values at fault stage
  //-----------------------------------------------------------------------
  // assume switch info is set up here instead of reading from the input file
  int sw2_2 = 5; // 6-1
  int sw3_2 = 6; // 7-1

  /*factory.setMode(FY); 
  gridpack::mapper::FullMatrixMap<DSNetwork> fy11ybusMap(network);
  boost::shared_ptr<gridpack::math::Matrix> fy11ybus = fy11ybusMap.mapToMatrix();
  fy11ybus->print();*/

  boost::shared_ptr<gridpack::math::Matrix> fy11ybus(ybus->clone());
  gridpack::ComplexType x(0.0, -1e7);
  fy11ybus->setElement(sw2_2, sw2_2, -x);
  fy11ybus->ready();
  printf("\n=== fy11ybus: ============\n");
  fy11ybus->print();
/* 
  // Solve linear equations of fy11ybus * X = Y_c
  
  // Form reduced admittance matrix fy11: fy11 = Y_b * X
  boost::shared_ptr<gridpack::math::Matrix> fy11(multiply(*Y_b, *X)); 
  //fy11->print();

  // Update fy11: fy11 = Y_a + fy11
  fy11->add(*Y_a);
*/
  //-----------------------------------------------------------------------
  // Compute posfy11
  // Update ybus values at clear fault stage
  //-----------------------------------------------------------------------
  /*factory.setMode(POSFY); 
  gridpack::mapper::FullMatrixMap<DSNetwork> posfy11ybusMap(network);
  boost::shared_ptr<gridpack::math::Matrix> posfy11ybus = posfy11ybusMap.mapToMatrix();
  posfy11ybus->print();*/

  // Get the updating factor for posfy11 stage ybus
  gridpack::ComplexType myValue = factory.setFactor(sw2_2, sw3_2);
  boost::shared_ptr<gridpack::math::Matrix> posfy11ybus(ybus->clone());

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
/*    
  // Solve linear equations of posfy11ybus * X = Y_c
  
  // Form reduced admittance matrix posfy11: posfy11 = Y_b * X
  boost::shared_ptr<gridpack::math::Matrix> posfy11(multiply(*Y_b, *X)); 
  //fy11->print();

  // Update posfy11: posfy11 = Y_a + posfy11
  posfy11->add(*Y_a);
*/  
}

