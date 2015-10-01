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
 * @date   2015-10-01 07:26:47 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <list>
#include <boost/bind.hpp>
#include "gridpack/utilities/exception.hpp"
#include "gridpack/parallel/parallel.hpp"
#include "variable.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

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

  void visit(go::Variable& var)
  {
    std::cout << var.name() << std::endl;
  }
  void visit(go::RealVariable& var)
  {
    std::cout << var.name() << ": real"
              << ": (" << var.lowerBound() << ":" << var.upperBound() << ")"
              << std::endl;
  }

  void visit(go::IntegerVariable& var)
  {
    std::cout << var.name()  << ": integer"
              << ": (" << var.lowerBound() << ":" << var.upperBound() << ")"
              << std::endl;
  }

  void visit(go::BinaryVariable& var)
  {
    std::cout << var.name() << ": binary"
              << ": (" << var.lowerBound() << ":" << var.upperBound() << ")"
              << std::endl;
  }
};

BOOST_AUTO_TEST_SUITE( VectorIOTest )

BOOST_AUTO_TEST_CASE( create )
{
  std::list<go::VariablePtr> vlist;
  vlist.push_back(go::VariablePtr(new go::RealVariable(13.0)));
  vlist.push_back(go::VariablePtr(new go::RealVariable(0.0, -1.0, 1.0)));
  vlist.push_back(go::VariablePtr(new go::IntegerVariable(0, -1, 1)));
  vlist.push_back(go::VariablePtr(new go::BinaryVariable(1)));

  go::VariableCounter cnt;
  for_each(vlist.begin(), vlist.end(),
           boost::bind(&go::Variable::accept, _1, boost::ref(cnt)));
  BOOST_CHECK_EQUAL(cnt.numVar, 4);
  BOOST_CHECK_EQUAL(cnt.numReal, 2);
  BOOST_CHECK_EQUAL(cnt.numInt, 1);
  BOOST_CHECK_EQUAL(cnt.numBin, 1);

  {
    go::VariableTable vt(std::cout);
    for_each(vlist.begin(), vlist.end(),
             boost::bind(&go::Variable::accept, _1, boost::ref(vt)));
  }
}

// BOOST_AUTO_TEST_CASE( serialization )
// {

// }


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
