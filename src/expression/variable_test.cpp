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
 * @date   2015-08-28 10:04:47 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <list>
#include <boost/bind.hpp>
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
  void visit(const go::RealVariable& var)
  {
    std::cout << var.name() << ": real"
              << ": (" << var.lowerBound() << ":" << var.upperBound() << ")"
              << std::endl;
  }

  void visit(const go::IntegerVariable& var)
  {
    std::cout << var.name()  << ": integer"
              << ": (" << var.lowerBound() << ":" << var.upperBound() << ")"
              << std::endl;
  }

  void visit(const go::BinaryVariable& var)
  {
    std::cout << var.name() << ": binary"
              << ": (" << var.lowerBound() << ":" << var.upperBound() << ")"
              << std::endl;
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
    std::cout << (*i)->name() << std::endl;
  }

  for_each(vlist.begin(), vlist.end(),
           boost::bind(&go::Variable::accept, _1, boost::ref(vp)));
  
  return 0;
}
