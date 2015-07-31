// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   variable_test.cpp
 * @author William A. Perkins
 * @date   2015-07-30 15:23:25 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <list>
#include "variable.hpp"

namespace go = gridpack::optimization;

// -------------------------------------------------------------
//  class VariableNamePrinter
// -------------------------------------------------------------
class VariablePrinter 
  : public go::VariableVisitor
{
public:

  /// Default constructor.
  VariablePrinter(void)
    : VariableVisitor()
  {}

  /// Destructor
  ~VariablePrinter(void)
  {}

  void visit(const go::Variable& var)
  {
    std::cout << var.name() << std::endl;
  }
  
};

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  std::list<go::VariablePtr> vlist;
  vlist.push_back(go::VariablePtr(new go::RealVariable(13.0)));
  vlist.push_back(go::VariablePtr(new go::RealVariable(0.0, -1.0, 1.0)));
  vlist.push_back(go::VariablePtr(new go::IntegerVariable(0, -1, 1)));
  vlist.push_back(go::VariablePtr(new go::BinaryVariable(1)));

  VariablePrinter vp;
  for (std::list<go::VariablePtr>::iterator i = vlist.begin();
       i != vlist.end(); ++i) {
    (*i)->accept(vp);
  }
  return 0;
}
