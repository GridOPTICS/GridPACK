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
 * @date   2015-08-28 16:50:04 d3g096
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

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  go::Optimizer opt;

  std::vector<go::VariablePtr> vars;

  // id0-8
  for (int i = 0; i < 9; ++i) {
    go::VariablePtr v(new go::BinaryVariable(1));
    vars.push_back(v);
  }

  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 250.0))); // id9
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 100.0)));
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0,  50.0)));
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 250.0)));
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 100.0)));
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0,  50.0)));
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 250.0)));
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0, 100.0)));
  vars.push_back(go::VariablePtr(new go::RealVariable(0, 0.0,  50.0))); // id17

  for (std::vector<go::VariablePtr>::iterator i = vars.begin();
       i != vars.end(); ++i) {
    opt.addVariable(*i);
  }

  opt.addConstraint( vars[0] == 1 );
  opt.addConstraint( vars[9] == 150 );
  opt.addConstraint( -10000*vars[3] + vars[12] <= 0 );
  opt.addConstraint( -1*vars[3] + vars[12] >= 0 );
  opt.addConstraint( -1*vars[0] + -9999*vars[3] <= 0 );
  opt.addConstraint( -1*vars[0] + vars[3] + -10000*vars[6] <= 0 );
  opt.addConstraint( -1*vars[0] + -9999*vars[3] <= -1 );

  opt.solve();

  return 0;
}

