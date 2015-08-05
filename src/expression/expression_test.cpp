// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   expression_test.cpp
 * @author William A. Perkins
 * @date   2015-07-31 13:38:55 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <list>
#include "expression.hpp"

namespace go = gridpack::optimization;

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  // 4*6 + 2
  go::ExpressionPtr four(new go::IntegerConstant(4));
  go::ExpressionPtr six(new go::IntegerConstant(6));
  go::ExpressionPtr two(new go::IntegerConstant(2));
  go::VariablePtr A(new go::RealVariable(13.0));
  go::VariablePtr B(new go::RealVariable(0.0, -1.0, 1.0));
  go::VariablePtr C(new go::IntegerVariable(0, -1, 1));
  go::ExpressionPtr junk;
  junk = A*4.0;
  junk->evaluate();
  std::cout << std::endl;
  junk = four*six + two*C;
  junk->evaluate();
  std::cout << std::endl;
  junk = four * 6 + 2;
  junk->evaluate();
  std::cout << std::endl;
  return 0;
}
