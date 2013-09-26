// -------------------------------------------------------------
/**
 * @file   dynsim_app.cpp
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
#include "gridpack/parser/PTI23_parser.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/parser/Parser.hpp"
#include "gridpack/applications/dynamic_simulation/dynsim_app.hpp"
#include "gridpack/applications/dynamic_simulation/dynsim_factory.hpp"

// Calling program for dynamic simulation application

/**
 * Basic constructor
 */
gridpack::dynamic_simulation::DynSimApp::DynSimApp(void)
{
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::DynSimApp::~DynSimApp(void)
{
}

/**
 * Execute application
 */
void gridpack::dynamic_simulation::DynSimApp::execute(void)
{
  gridpack::parallel::Communicator world;
  boost::shared_ptr<DynSimNetwork> network(new DynSimNetwork(world));

  // read configuration file
  gridpack::utility::Configuration config;

  // load input file
  gridpack::parser::PTI23_parser<DynSimNetwork> parser(network);
  parser.getCase("SJin_3g9b.raw");
  parser.createNetwork();

  // partition network
  network->partition();

  // create factory
  gridpack::dynsim::DynSimFactory factory(network);
  factory.load();

  // set network components using factory
  factory.setComponents();

  // set YBus components so that you can create Y matrix  
  factory.setYBus();

  factory.setMode(YBus);
  gridpack::mapper::FullMatrixMap<DynSimNetwork> mMap(network);
  boost::shared_ptr<gridpack::math::Matrix> ybus = mMap.mapToMatrix();
  ybus->print();

  // assume switch info is set up here instead of reading from the input file
  int sw2_2 = 5; // 6-1
  int sw3_2 = 6; // 7-1

  // Form constant impedance load admittance yl for all buses and add it to 
  // system Y matrix: ybus = ybus + yl
  factory.setMode(YL); 
  gridpack::mapper::FullMatrixMap<DynSimNetwork> mMap(network);
  boost::shared_ptr<gridpack::math::Matrix> ybus = mMap.mapToMatrix();
  ybus->print();
 
  // Construct permutation matrix perm by checking for multiple generators at a bus
  factory.setMode(PERM);
  gridpack::mapper::FullMatrixMap<DynSimNetwork> mMap(network);
  boost::shared_ptr<gridpack::math::Matrix> perm = mMap.mapToMatrix();
  perm->print(); 

  // Form a transposed matrix of perm
  boost::shared_ptr<gridpack::math::Matrix> permTrans;
  gridpack::math::Matrix::transpose(perm, permTrans);
 
  // Construct matrix Y_a using extracted xd and ra from gen data, 
  // and construct its diagonal matrix diagY_a
  factory.setMode(YA);
  gridpack::mapper::FullMatrixMap<DynSimNetwork> mMap(network);
  boost::shared_ptr<gridpack::math::Matrix> Y_a = mMap.mapToMatrix();
  Y_a->print(); 
  boost::shared_ptr<gridpack::math::Matrix> diagY_a = identity();
  boost::shared_ptr<gridpack::math::Vector> vV_a = diagonal(V_a); 
  diagY_a.muitiply_diagonal(vV_a);
  diagY_a->print(); 

  // Construct matrix Ymod: Ymod = diagY_a * permTrans
  boost::shared_ptr<gridpack::math::Matrix> Ymod = diagY_a.multiply(permTrans);
  
  // Form matrix Y_b: Y_b(1:ngen, jg) = -Ymod, where jg represents the 
  // corresponding index sets of buses that the generators are connected to 
  /*factory.setMode(YB);
  gridpack::mapper::FullMatrixMap<DynSimNetwork> mMap(network);
  boost::shared_ptr<gridpack::math::Matrix> Y_b = mMap.mapToMatrix();
  Y_b->print(); */
  factory.setMode(PMatrix);
  gridpack::mapper::FullMatrixMap<DynSimNetwork> mMap(network);
  boost::shared_ptr<gridpack::math::Matrix> P = mMap.mapToMatrix();
  boost::shared_ptr<gridpack::math::Matrix> Y_b = P.multiply(Ymod);
  Y_b->print();
  
  // Construct Y_b's transposed matrix Y_c: Y_c = Y_b'
  boost::shared_ptr<gridpack::math::Matrix> Y_c;
  gridpack::math::Matrix::transpose(Y_b, Y_c);
  Y_c->print();
   
  // Form matrix permYmod
  boost::shared_ptr<gridpack::math::Matrix> permYmod;
  gridpack::math::Matrix::multiply(perm, Ymod, permYmod);

  // Update ybus: ybus = ybus+permYmod
  //ybus = ybus + permYmod;
  ybus.add(permYmod);

  // Solve linear equations of ybus * X = Y_c

  // Form reduced admittance matrix prefy11: prefy11 = Y_b * X

  // Update prefy11: prefy11 = Y_a + prefy11

}

