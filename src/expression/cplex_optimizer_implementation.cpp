// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   cplex_optimizer_implementation.cpp
 * @author William A. Perkins
 * @date   2015-08-28 16:44:31 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <boost/bind.hpp>
#include "cplex_optimizer_implementation.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class CPlexVarVisitor
// -------------------------------------------------------------
class CPlexVarVisitor 
  : public VariableVisitor
{
public:

  /// Default constructor.
  CPlexVarVisitor(std::ostream& o)
    : VariableVisitor(),
      p_stream(o)
  {
  }

  /// Destructor
  ~CPlexVarVisitor(void)
  {
  }

protected:

  std::ostream& p_stream;

};

// -------------------------------------------------------------
// class CPlexVarBoundsLister
// -------------------------------------------------------------
class CPlexVarBoundsLister 
  : public CPlexVarVisitor
{
public:

  /// Default constructor.
  CPlexVarBoundsLister(std::ostream& o)
    : CPlexVarVisitor(o)
  {}

  /// Destructor
  ~CPlexVarBoundsLister(void)
  {}

  void visit(const Variable& var)
  {
    BOOST_ASSERT(false);
  }

  void visit(const RealVariable& var)
  {
    std::string s;
    s += "";
    if (var.lowerBound() > var.veryLowValue) {
      s += boost::str(boost::format("%8.4g") % var.lowerBound());
    } else {
      s += "        ";
    }
    s += " <= ";
    s += boost::str(boost::format("%-4.4s") % var.name());
    s += " <= ";
    if (var.lowerBound() < var.veryHighValue) {
      s += boost::str(boost::format("%8.4g") % var.upperBound());
    } else {
      s += "        ";
    }
    this->p_stream << s << std::endl;
  }
    
  void visit(const IntegerVariable& var)
  {
    std::string s;
    s = "";
    if (var.lowerBound() > var.veryLowValue) {
      s += boost::str(boost::format("%8d") % var.lowerBound());
    } else {
      s += "        ";
    }
    s += " <= ";
    s += boost::str(boost::format("%-4.4s") % var.name());
    s += " <= ";
    if (var.lowerBound() < var.veryHighValue) {
      s += boost::str(boost::format("%8d") % var.upperBound());
    } else {
      s += "        ";
    }
    this->p_stream << s << std::endl;
  }
};

// -------------------------------------------------------------
//  class CPlexBinVarLister
// -------------------------------------------------------------
class CPlexBinVarLister 
  : public CPlexVarVisitor
{
public:

  /// Default constructor.
  CPlexBinVarLister(std::ostream& o)
    : CPlexVarVisitor(o)
  {
  }

  /// Destructor
  ~CPlexBinVarLister(void)
  {
  }

  void visit(const BinaryVariable& var)
  {
    this->p_stream << " " << var.name();
  };
};


// -------------------------------------------------------------
//  class CPlexOptimizerImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// CPlexOptimizerImplementation::p_solve
// -------------------------------------------------------------
void
CPlexOptimizerImplementation::p_solve(void)
{
  std::cout << "\\* Problem name: Test\\*" << std::endl;

  std::cout << "Subject To" << std::endl;
  for (std::vector<ConstraintPtr>:: iterator i = p_constraints.begin();
       i != p_constraints.end(); ++i) {
    std::cout << (*i)->render() << std::endl;
  }
  std::cout << std::endl;   


  std::cout << "Bounds" << std::endl;
  {
    CPlexVarBoundsLister v(std::cout);
    for_each(p_variables.begin(), p_variables.end(),
             boost::bind(&Variable::accept, _1, boost::ref(v)));

  }
  std::cout << std::endl;

  std::cout << "Binaries" << std::endl;
  std::cout << "    ";
  {
    CPlexBinVarLister v(std::cout);
    for_each(p_variables.begin(), p_variables.end(),
             boost::bind(&Variable::accept, _1, boost::ref(v)));

  }
  std::cout << std::endl << std::endl;


  std::cout << "End" << std::endl;
}
    
  


} // namespace optimization
} // namespace gridpack

