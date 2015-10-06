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
 * @date   2015-10-06 11:29:34 d3g096
 * 
 * @brief  Unit tests for gridpack::optimization::Optimizer class
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <vector>

#include "optimizer.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

namespace go = gridpack::optimization;
namespace gp = gridpack::parallel;

BOOST_AUTO_TEST_SUITE( Optimization )

// -------------------------------------------------------------
// UNIT TEST: flow
// This is from a GLPK example. A simple network flow optimization. 
// -------------------------------------------------------------
BOOST_AUTO_TEST_CASE( flow )
{
  gp::Communicator world;
  gp::Communicator self(world.self());
  go::Optimizer opt(self);

  std::vector<go::VariablePtr> vars;

  vars.push_back(go::VariablePtr(new go::RealVariable(0))); // 0. total flow
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 14.0))); // 1. flow cap from 1->2 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 23.0))); // 2. flow cap from 1->4 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 11.0))); // 3. flow cap from 5->2 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 10.0))); // 4. flow cap from 2->3 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0,  9.0))); // 5. flow cap from 2->4 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 12.0))); // 6. flow cap from 3->5 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 18.0))); // 7. flow cap from 3->8 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 26.0))); // 8. flow cap from 4->5 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 25.0))); // 9. flow cap from 5->6 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0,  4.0))); //10. flow cap from 5->7 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0,  7.0))); //11. flow cap from 6->7 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0,  8.0))); //12. flow cap from 6->8 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 15.0))); //13. flow cap from 7->9 
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 20.0))); //14. flow cap from 8->9 

  for (std::vector<go::VariablePtr>::iterator i = vars.begin();
       i != vars.end(); ++i) {
    opt.addVariable(*i);
  }

  opt.addConstraint( -1.0 * vars[1] - vars[2] + vars[0] == -0 );
  opt.addConstraint( + vars[1] + vars[3] - vars[4] - vars[5] == -0 );
  opt.addConstraint( + vars[4] - vars[6] - vars[7] == -0 );
  opt.addConstraint( + vars[2] + vars[5] - vars[8] == -0 );
  opt.addConstraint( - vars[3] + vars[6] + vars[8] - vars[9] - vars[10] == -0 );
  opt.addConstraint( + vars[9] - vars[11] - vars[12] == -0 );
  opt.addConstraint( + vars[10]+ vars[11] - vars[13] == -0 );
  opt.addConstraint( + vars[7] + vars[12] - vars[14] == -0 );
  opt.addConstraint( + vars[13]+ vars[14] - vars[0] == -0 );
  
  go::ExpressionPtr obj(new go::VariableExpression(vars[0]));
  opt.addToObjective(obj);

  opt.maximize();

  // if an optimizer can be run, the answer is 29

#if defined(HAVE_CPLEX) || defined(HAVE_GLPK)
  go::GetVariableInitial g;
  vars[0]->accept(g);
  BOOST_CHECK_EQUAL(g.value(), 29.0);
#endif
}

// -------------------------------------------------------------
// UNIT TEST: uc
// -------------------------------------------------------------
BOOST_AUTO_TEST_CASE( uc )
{
  gp::Communicator world;
  gp::Communicator self(world.self());
  go::Optimizer opt(self);

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

#if defined(HAVE_GLPK) 
  // This is a quadratic problem. It will not be understood by GLPK.
  BOOST_CHECK_THROW(opt.minimize(), gridpack::Exception);
#elif defined(HAVE_CPLEX)
  // CPLEX should be able to handle this. I don't know what the answer
  // should be though.
  opt.minimize();

#else

  // This should just print the problem
  opt.minimize();

#endif
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
  gridpack::parallel::Environment env(argc, argv);
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  return result;
}

