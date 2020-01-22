// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   optimizer_test.cpp
 * @author William A. Perkins
 * @date   2017-03-22 08:49:29 d3g096
 * 
 * @brief  Unit tests for gridpack::optimization::Optimizer class
 * 
 * 
 */
// -------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <vector>

#include "gridpack/environment/environment.hpp"
#include "optimizer.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

namespace go = gridpack::optimization;
namespace gp = gridpack::parallel;

/// The configuration used for these tests
static gridpack::utility::Configuration::CursorPtr test_config;

BOOST_AUTO_TEST_SUITE( Optimization )

// -------------------------------------------------------------
// UNIT TEST: flow
// This is from a GLPK example. A simple network flow optimization. 
// -------------------------------------------------------------
BOOST_AUTO_TEST_CASE( flow )
{
  gp::Communicator world;
  int nproc(world.size());
  int me(world.rank());

  static const int netnodes(9);
  static const int netedges(14);

  BOOST_REQUIRE(nproc <= netnodes);
  BOOST_REQUIRE(test_config);

  gridpack::utility::Configuration::CursorPtr
    flow_config(test_config->getCursor("FlowTest"));

  go::Optimizer opt(world);
  opt.configure(flow_config);

  std::vector<bool> iownnode(netnodes+1);
  for (int i = 1; i < netnodes+1; ++i) {
    iownnode[i] = ((i % nproc ) == me);
  }

  std::vector<go::VariablePtr> vars(netedges+2);

  go::VariablePtr v;

  // flux 0. total in flux (everybody makes this so it can be checked)
  v.reset(new go::RealVariable(0.0)); 
  v->name("TotalInFlux"); 
  vars[0] = v;

  // flux 1. from 1->2 
  if (iownnode[1] || iownnode[2]) {
    v.reset(new go::RealVariable(0, 0.0, 14.0)); 
    v->name("Flux01");
    vars[1] = v;
  }
  // flux 2. from 1->4 
  if (iownnode[1] || iownnode[4]) {
    v.reset(new go::RealVariable(0, 0.0, 23.0)); 
    v->name("Flux02");
    vars[2] = v;
  }
  // 3. from 5->2 
  if (iownnode[5] || iownnode[2]) {
    v.reset(new go::RealVariable(0, 0.0, 11.0)); 
    v->name("Flux03");
    vars[3] = v;
  }
  // 4. from 2->3 
  if (iownnode[2] || iownnode[3]) {
    v.reset(new go::RealVariable(0, 0.0, 10.0)); 
    v->name("Flux04");
    vars[4] = v;
  }
  // 5. from 2->4 
  if (iownnode[2] || iownnode[4]) {
    v.reset(new go::RealVariable(0, 0.0,  9.0)); 
    v->name("Flux05");
    vars[5] = v;
  }
  // 6. from 3->5 
  if (iownnode[3] || iownnode[5]) {
    v.reset(new go::RealVariable(0, 0.0, 12.0)); 
    v->name("Flux06");
    vars[6] = v;
  }
  // 7. from 3->8 
  if (iownnode[3] || iownnode[8]) {
    v.reset(new go::RealVariable(0, 0.0, 18.0)); 
    v->name("Flux07");
    vars[7] = v;
  }
  // 8. from 4->5 
  if (iownnode[4] || iownnode[5]) {
    v.reset(new go::RealVariable(0, 0.0, 26.0)); 
    v->name("Flux08");
    vars[8] = v;
  }
  // 9. from 5->6 
  if (iownnode[5] || iownnode[6]) {
    v.reset(new go::RealVariable(0, 0.0, 25.0)); 
    v->name("Flux09");
    vars[9] = v;
  }
  //10. from 5->7 
  if (iownnode[5] || iownnode[7]) {
    v.reset(new go::RealVariable(0, 0.0,  4.0)); 
    v->name("Flux10");
    vars[10] = v;
  }
  //11. from 6->7 
  if (iownnode[6] || iownnode[7]) {
    v.reset(new go::RealVariable(0, 0.0,  7.0)); 
    v->name("Flux11");
    vars[11] = v;
  }
  //12. from 6->8 
  if (iownnode[6] || iownnode[8]) {
    v.reset(new go::RealVariable(0, 0.0,  8.0)); 
    v->name("Flux12");
    vars[12] = v;
  }
  //13. from 7->9 
  if (iownnode[7] || iownnode[9]) {
    v.reset(new go::RealVariable(0, 0.0, 15.0)); 
    v->name("Flux13");
    vars[13] = v;
  }
  //14. from 8->9 
  if (iownnode[8] || iownnode[9]) {
    v.reset(new go::RealVariable(0, 0.0, 20.0)); 
    v->name("Flux14");
    vars[14] = v;
  }

  // 15. total out flux
  if (iownnode[9]) {
    v.reset(new go::RealVariable(0.0)); 
    v->name("TotalOutFlux");
    vars[15] = v;
  }

  for (std::vector<go::VariablePtr>::iterator i = vars.begin();
       i != vars.end(); ++i) {
    if (*i) { opt.addVariable(*i); }
  }

  // add a global constraint (inflow = outflow)
  go::ExpressionPtr empty;
  go::ConstraintPtr c( empty == 0.0 );
  c->name("inout");
  opt.createGlobalConstraint(c->name(), c);
  c.reset();

  if ( iownnode[1] ) {
    c = ( - vars[1] - vars[2] + vars[0] == 0); 
    c->name("Node1");  
    opt.addConstraint( c );
    go::ExpressionPtr v(new go::VariableExpression(vars[0]));
    opt.addToObjective(v);
    opt.addToGlobalConstraint("inout", v);
  }
  if ( iownnode[2] ) {
    c = ( + vars[1] + vars[3] - vars[4] - vars[5] == 0);  
    c->name("Node2");  
    opt.addConstraint( c );
  }
  if ( iownnode[3] ) {
    c = ( + vars[4] - vars[6] - vars[7] == 0 );
    c->name("Node3");
    opt.addConstraint( c );
  }
  if ( iownnode[4] ) {
    c = ( + vars[2] + vars[5] - vars[8] == 0 );           
    c->name("Node4");  
    opt.addConstraint( c );
  }
  if ( iownnode[5] ) {
    c = ( - vars[3] + vars[6] + vars[8] - vars[9] - vars[10] == 0 ); 
    c->name("Node5");  
    opt.addConstraint( c );
  }
  if ( iownnode[6] ) {
    c = ( + vars[9] - vars[11] - vars[12] == 0 );         
    c->name("Node6");  
    opt.addConstraint( c );
  }
  if ( iownnode[7] ) {
    c = ( + vars[10]+ vars[11] - vars[13] == 0 ); 
    c->name("Node7");  
    opt.addConstraint( c );
  }
  if ( iownnode[8] ) {
    c = ( + vars[7] + vars[12] - vars[14] == 0 ); 
    c->name("Node8");  
    opt.addConstraint( c );
  }
  if ( iownnode[9] ) {
    c = ( + vars[13]+ vars[14] + vars[15] == 0 );
    c->name("Node9");
    opt.addConstraint( c ); 
    go::ExpressionPtr v(new go::VariableExpression(vars[15]));
    opt.addToGlobalConstraint("inout", v);
  }

  opt.maximize();

#if defined(HAVE_GLPK) || defined(HAVE_CPLEX)
  // if an optimizer can be run, the answer is 29
  go::GetVariableInitial g;
  vars[0]->accept(g);
  BOOST_CHECK_EQUAL(g.value(), 29.0);
#endif

  world.barrier();

}

// -------------------------------------------------------------
// UNIT TEST: uc
// -------------------------------------------------------------
BOOST_AUTO_TEST_CASE( uc )
{
  gp::Communicator world;
  gp::Communicator self(world.self());

  gridpack::utility::Configuration::CursorPtr
    uc_config(test_config->getCursor("UCTest"));

  go::Optimizer opt(self);
  opt.configure(uc_config);

  std::vector<go::VariablePtr> vars;

  vars.push_back(go::VariablePtr(new go::RealVariable(150, 0, 250))); // V0
  vars.push_back(go::VariablePtr(new go::RealVariable(0,   0, 100))); 
  vars.push_back(go::VariablePtr(new go::RealVariable(0,   0,  50))); 
  vars.push_back(go::VariablePtr(new go::RealVariable(0,   0, 250))); 
  vars.push_back(go::VariablePtr(new go::RealVariable(0,   0, 100))); 

  vars.push_back(go::VariablePtr(new go::RealVariable(0,   0,  50))); // V5 
  vars.push_back(go::VariablePtr(new go::RealVariable(0,   0, 250)));
  vars.push_back(go::VariablePtr(new go::RealVariable(0,   0, 100))); 
  vars.push_back(go::VariablePtr(new go::RealVariable(0,   0,  50))); 
  vars.push_back(go::VariablePtr(new go::IntegerVariable(0, 0,  1)));

  vars.push_back(go::VariablePtr(new go::IntegerVariable(0, 0,  1))); // V10
  vars.push_back(go::VariablePtr(new go::IntegerVariable(0, 0,  1)));
  vars.push_back(go::VariablePtr(new go::IntegerVariable(0, 0,  1)));
  vars.push_back(go::VariablePtr(new go::IntegerVariable(0, 0,  1)));
  vars.push_back(go::VariablePtr(new go::IntegerVariable(0, 0,  1)));

  vars.push_back(go::VariablePtr(new go::IntegerVariable(0, 0,  1))); // V15
  vars.push_back(go::VariablePtr(new go::IntegerVariable(0, 0,  1)));
  vars.push_back(go::VariablePtr(new go::IntegerVariable(0, 0,  1)));

 
  for (std::vector<go::VariablePtr>::iterator i = vars.begin();
       i != vars.end(); ++i) {
    opt.addVariable(*i);
  }

  opt.addConstraint( vars[0] == 150.0 ); // C0
  opt.addConstraint( vars[9] == 1.0);
  opt.addConstraint( vars[3] - 10000*vars[12] <= 0 );
  opt.addConstraint( vars[3] - 150*vars[12] >= 0 );
  opt.addConstraint( vars[12] - vars[9] - 10000*vars[12] <= 0 );

  opt.addConstraint( vars[12] - vars[9] - 10000*vars[15] <= 0 ); // C5
  opt.addConstraint( -vars[9] + vars[12] - 10000*vars[12] <= -1 );
  opt.addConstraint( -vars[9] + vars[12] - 10000*vars[15] <= -1 );
  opt.addConstraint( vars[4] - 10000*vars[13] <= 0 );
  opt.addConstraint( vars[4] - 50*vars[13] >= 0 );

  opt.addConstraint( vars[13] - vars[10] - 10000*vars[13] <= 0 ); // C10
  opt.addConstraint( -vars[10] + vars[13] - 10000*vars[13] <= -1 );
  opt.addConstraint( vars[5] - 10000*vars[14] <= 0 );
  opt.addConstraint( vars[5] - 10*vars[14] >= 0 );
  opt.addConstraint( vars[14] - vars[11] - 10000*vars[14] <= 0 );

  opt.addConstraint( -vars[11] + vars[14] - 10000*vars[14] <= -1 ); // C15
  opt.addConstraint( -vars[11] + vars[14] - 10000*vars[17] <= -1 );
  opt.addConstraint( vars[3] + vars[4] + vars[5] == 300 );
  opt.addConstraint( vars[6] - 10000*vars[15] <= 0 );
  opt.addConstraint( vars[6] - 150*vars[15] >= 0 );

  opt.addConstraint( vars[15] - vars[12] - 10000*vars[15] <= 0 ); // C20
  opt.addConstraint( -vars[12] + vars[15] - 10000*vars[15] <= -1 );
  opt.addConstraint( vars[7] - 10000*vars[16] <= 0 );
  opt.addConstraint( vars[7] - 50*vars[16] >= 0 );
  opt.addConstraint( vars[16] - vars[13] - 10000*vars[16] <= 0 );

  opt.addConstraint( -vars[13] + vars[16] - 10000*vars[16] <= -1 ); // C25
  opt.addConstraint( vars[8] - 10000*vars[17] <= 0 );
  opt.addConstraint( vars[8] - 10*vars[17] >= 0 );
  opt.addConstraint( vars[17] - vars[14] - 10000*vars[17] <= 0 );
  opt.addConstraint( -vars[14] + vars[17] - 10000*vars[17] <= -1 );

  opt.addConstraint( vars[6] + vars[7] + vars[8] == 200 ); // C30


  opt.addToObjective(510.0*vars[9]);
  opt.addToObjective( 7.9*vars[0] );
  opt.addToObjective( 0.00344*(vars[0]^2) );
  opt.addToObjective( 310*vars[10] );
  opt.addToObjective(7.85*vars[1] );
  opt.addToObjective(0.00388*(vars[1]^2) );
  opt.addToObjective(78*vars[11] );
  opt.addToObjective(9.56*vars[2] );
  opt.addToObjective(0.01388*(vars[2]^2) );
  opt.addToObjective(510*vars[12] );
  opt.addToObjective(7.9*vars[3] );
  opt.addToObjective(0.00344*(vars[3]^2) );
  opt.addToObjective(310*vars[13] );
  opt.addToObjective(7.85*vars[4] );
  opt.addToObjective(0.00388*(vars[4]^2) );
  opt.addToObjective(78*vars[14] );
  opt.addToObjective(9.56*vars[5] );
  opt.addToObjective(0.01388*(vars[5]^2) );
  opt.addToObjective(510*vars[15] );
  opt.addToObjective(7.9*vars[6] );
  opt.addToObjective(0.00344*(vars[6]^2));
  opt.addToObjective(310*vars[16] );
  opt.addToObjective(7.85*vars[7] );
  opt.addToObjective(0.00388*(vars[7]^2));
  opt.addToObjective(78*vars[17] );
  opt.addToObjective(9.56*vars[8] );
  opt.addToObjective(0.01388*(vars[8]^2));

  // CPLEX should be able to handle this. I don't know what the answer
  // should be though.
  opt.minimize();

  world.barrier();

}

BOOST_AUTO_TEST_CASE (function)
{
  gp::Communicator world;
  gp::Communicator self(world.self());

  gridpack::utility::Configuration::CursorPtr
    uc_config(test_config->getCursor("FunctionTest"));

  go::Optimizer opt(self);
  opt.configure(uc_config);

  go::VariablePtr A(new go::RealVariable(M_PI, 0, 2*M_PI));
  go::VariablePtr B(new go::RealVariable(M_PI, 0, 2*M_PI));

  opt.addVariable(A);
  opt.addVariable(B);

  opt.addConstraint( go::sin(A) + go::cos(B) > 0.0 );
  opt.addToObjective( go::sin(A) + go::sin(B) );
  
  opt.minimize();

  world.barrier();
}  

BOOST_AUTO_TEST_SUITE_END()

// -------------------------------------------------------------
// init_function
// -------------------------------------------------------------
bool init_function()
{
  return true;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gridpack::Environment env(argc, argv);
  gridpack::parallel::Communicator world;

  boost::scoped_ptr<gridpack::utility::Configuration> 
    config(gridpack::utility::Configuration::configuration());
  
  config->enableLogging();
  config->open("gridpack.xml", world);

  test_config = config->getCursor("GridPACK.OptimizerTests");

  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );

   return result;
}

