// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   uc_test.cpp
 * @author William A. Perkins
 * @date   2015-09-28 15:53:43 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <vector>
#include "optimizer.hpp"

namespace go = gridpack::optimization;
namespace gp = gridpack::parallel;

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  gp::Environment env(argc, argv);
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

  opt.minimize();
  
  return 0;
}


