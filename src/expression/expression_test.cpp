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
 * @date   2015-10-07 14:45:05 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/bind.hpp>
#include "expression.hpp"

#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace go = gridpack::optimization;

// -------------------------------------------------------------
//  class ExpressionVariableChecker
// -------------------------------------------------------------
class ExpressionVariableChecker 
  : public go::ExpressionVisitor
{
public:

  /// Default constructor.
  ExpressionVariableChecker(void)
    : go::ExpressionVisitor()
  {}

  /// Destructor
  ~ExpressionVariableChecker(void)
  {}

  void visit(go::VariableExpression& e)
  {
    go::VariablePtr v(e.var());
    void *vaddr = v.get();
    var[vaddr] = v->name();
  } 

  std::map<void *, std::string> var;

};



// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{

  std::vector<go::VariablePtr> vars;
  go::VariablePtr A(new go::RealVariable(13.0));
  go::VariablePtr B(new go::RealVariable(0.0, -1.0, 1.0));
  go::VariablePtr C(new go::IntegerVariable(0, -1, 1));
  vars.push_back(A);
  vars.push_back(B);
  vars.push_back(C);

  go::ExpressionPtr four(new go::IntegerConstant(4));
  go::ExpressionPtr six(new go::IntegerConstant(6));
  go::ExpressionPtr two(new go::IntegerConstant(2));

  std::vector<go::ExpressionPtr> exprs;

  exprs.push_back( A*4.0 ); 
  exprs.push_back( four*six + two*C );
  exprs.push_back( four*(6*C + 2*A) );
  exprs.push_back( exprs[0] + exprs[1] + exprs[2] );
  exprs.push_back( 6*C + A/4 );
  exprs.push_back( (6.0*C + A)/4.6 ); 
  exprs.push_back( 6*B + 2*(A^2) ); 
  exprs.push_back( 6*C + ((2*A)^2) ); 

  std::for_each(exprs.begin(), exprs.end(), 
                boost::bind(&go::Expression::evaluate, _1));
  std::cout << std::endl;

  std::vector<go::ConstraintPtr> cons;
  cons.push_back( A < 4 );
  cons.push_back( exprs.back() >= 6 );
  cons.push_back( C == 3 );

  std::for_each(cons.begin(), cons.end(), 
                boost::bind(&go::Expression::evaluate, _1));
  std::cout << std::endl;

  std::string buf;
  {
    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa & vars & exprs & cons;
    buf = oss.str();
  }
  {
    std::vector<go::VariablePtr> invars;
    std::vector<go::ExpressionPtr> inexprs;
    std::vector<go::ConstraintPtr> incons;
    std::istringstream iss(buf);
    boost::archive::binary_iarchive ia(iss);
    ia & invars & inexprs & incons;

    std::for_each(inexprs.begin(), inexprs.end(), 
                  boost::bind(&go::Expression::evaluate, _1));
    std::cout << std::endl;

    std::for_each(incons.begin(), incons.end(), 
                  boost::bind(&go::Expression::evaluate, _1));

    ExpressionVariableChecker vcheck;
    std::for_each(inexprs.begin(), inexprs.end(), 
                  boost::bind(&go::Expression::accept, _1, boost::ref(vcheck)));
    std::for_each(incons.begin(), incons.end(), 
                  boost::bind(&go::Expression::accept, _1, boost::ref(vcheck)));
 
    std::cout << vcheck.var.size() << std::endl;
    for (std::map<void *, std::string>::iterator i = vcheck.var.begin();
         i != vcheck.var.end(); ++i) {
      std::cout << i->first << ": " << i->second 
                << ": " << ((go::Variable *)(i->first))->name()
                << std::endl;
    }
  }

  return 0;
}
