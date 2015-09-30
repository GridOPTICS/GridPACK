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
 * @date   2015-09-28 14:59:06 d3g096
 * 
 * @brief  Unit tests for gridpack::optimization::Optimizer class
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

  return 0;
}

