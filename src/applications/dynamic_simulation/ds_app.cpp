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

  // assume switch info is set up here instead of reading from the input file
  int sw2_2 = 5; // 6-1
  int sw3_2 = 6; // 7-1

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

/*  // Construct matrix Y_a using extracted xd and ra from gen data, 
  // and construct its diagonal matrix diagY_a
  factory.setMode(YA);
  gridpack::mapper::FullMatrixMap<DSNetwork> yaMap(network);
  boost::shared_ptr<gridpack::math::Matrix> Y_a = yaMap.mapToMatrix();
  Y_a->print(); 

  boost::shared_ptr<gridpack::math::Matrix> diagY_a;
  diagY_a->gridpack::math::Matrix::identity();

  boost::shared_ptr<gridpack::math::Vector> vY_a(diagonal(*Y_a)); 

  //diagY_a->gridpack::math::Matrix::muitiplyDiagonal(*vY_a); //???
  //diagY_a->muitiplyDiagonal(*vY_a); //???
  diagY_a->print(); 
*/
/*  // Construct matrix Ymod: Ymod = diagY_a * permTrans
  //boost::shared_ptr<gridpack::math::Matrix> Ymod(multiply(*diagY_a, *permTrans)); //???
 
  // Form matrix Y_b: Y_b(1:ngen, jg) = -Ymod, where jg represents the 
  // corresponding index sets of buses that the generators are connected to 
  factory.setMode(PMatrix);
  gridpack::mapper::FullMatrixMap<DSNetwork> pMap(network);
  boost::shared_ptr<gridpack::math::Matrix> P = pMap.mapToMatrix();
  //boost::shared_ptr<gridpack::math::Matrix> Y_b(multiply(*P, *Ymod)); //???
  Y_b->print();
  
  // Construct Y_b's transposed matrix Y_c: Y_c = Y_b'
  boost::shared_ptr<gridpack::math::Matrix> Y_c(transpose(*Y_b));
  Y_c->print();
   
  // Form matrix permYmod
  //boost::shared_ptr<gridpack::math::Matrix> permYmod(multiply(*perm, *Ymod)); //???

  // Update ybus: ybus = ybus+permYmod
  //ybus = ybus + permYmod;
  ybus->add(*permYmod);
*/
  // Solve linear equations of ybus * X = Y_c

  // Form reduced admittance matrix prefy11: prefy11 = Y_b * X

  // Update prefy11: prefy11 = Y_a + prefy11

}

