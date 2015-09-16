// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   lpfile_optimizer_implementation.cpp
 * @author William A. Perkins
 * @date   2015-09-16 08:20:56 d3g096
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
#include "lpfile_optimizer_implementation.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class LPFileVarVisitor
// -------------------------------------------------------------
class LPFileVarVisitor 
  : public VariableVisitor
{
public:

  /// Default constructor.
  LPFileVarVisitor(std::ostream& o)
    : VariableVisitor(),
      p_stream(o)
  {
  }

  /// Destructor
  ~LPFileVarVisitor(void)
  {
  }

protected:

  std::ostream& p_stream;

};

// -------------------------------------------------------------
// class LPFileVarBoundsLister
// -------------------------------------------------------------
class LPFileVarBoundsLister 
  : public LPFileVarVisitor
{
public:

  /// Default constructor.
  LPFileVarBoundsLister(std::ostream& o)
    : LPFileVarVisitor(o)
  {}

  /// Destructor
  ~LPFileVarBoundsLister(void)
  {}

  void visit(const Variable& var)
  {
    BOOST_ASSERT(false);
  }

  void visit(const RealVariable& var)
  {
    std::string s;
    s += "";
    if (var.bounded()) {
      if (var.lowerBound() > var.veryLowValue) {
        s += boost::str(boost::format("%8.4g") % var.lowerBound());
        s += " <= ";
      } else {
        s += "        ";
        s += "    ";
      }
      s += boost::str(boost::format("%-4.4s") % var.name());
      if (var.upperBound() < var.veryHighValue) {
        s += " <= ";
        s += boost::str(boost::format("%8.4g") % var.upperBound());
      } else {
        s += "        ";
        s += "    ";
      }
    } else {
      s += "        ";
      s += "    ";
      s += boost::str(boost::format("%-4.4s") % var.name());
      s += " free";
    }
    this->p_stream << s << std::endl;
  }
    
  void visit(const IntegerVariable& var)
  {
    std::string s;
    s = "";
    if (var.lowerBound() > var.veryLowValue) {
      s += boost::str(boost::format("%8d") % var.lowerBound());
      s += " <= ";
    } else {
      s += "        ";
      s += "    ";
    }
    s += boost::str(boost::format("%-4.4s") % var.name());
    if (var.upperBound() < var.veryHighValue) {
      s += " <= ";
      s += boost::str(boost::format("%8d") % var.upperBound());
    } else {
      s += "    ";
      s += "        ";
    }
    this->p_stream << s << std::endl;
  }
};

// -------------------------------------------------------------
//  class LPFileGenVarLister
// -------------------------------------------------------------
class LPFileGenVarLister 
  : public LPFileVarVisitor
{
public:

  /// Default constructor.
  LPFileGenVarLister(std::ostream& o)
    : LPFileVarVisitor(o)
  {
  }

  /// Destructor
  ~LPFileGenVarLister(void)
  {
  }

  void visit(const RealVariable& var)
  {
    this->p_stream << " " << var.name();
  };
};

// -------------------------------------------------------------
//  class LPFileIntVarLister
// -------------------------------------------------------------
class LPFileIntVarLister 
  : public LPFileVarVisitor
{
public:

  /// Default constructor.
  LPFileIntVarLister(std::ostream& o)
    : LPFileVarVisitor(o)
  {
  }

  /// Destructor
  ~LPFileIntVarLister(void)
  {
  }

  void visit(const IntegerVariable& var)
  {
    this->p_stream << " " << var.name();
  };
};

// -------------------------------------------------------------
//  class LPFileBinVarLister
// -------------------------------------------------------------
class LPFileBinVarLister 
  : public LPFileVarVisitor
{
public:

  /// Default constructor.
  LPFileBinVarLister(std::ostream& o)
    : LPFileVarVisitor(o)
  {
  }

  /// Destructor
  ~LPFileBinVarLister(void)
  {
  }

  void visit(const BinaryVariable& var)
  {
    this->p_stream << " " << var.name();
  };
};

// -------------------------------------------------------------
//  class LPFileOptimizerImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// LPFileOptimizerImplementation::p_solve
// -------------------------------------------------------------
void
LPFileOptimizerImplementation::p_solve(const p_optimizeMethod& method)
{
  std::cout << "\\* Problem name: Test\\*" << std::endl << std::endl;

  switch (method) {
  case Maximize:
    std::cout << "Maximize" << std::endl;
    break;
  case Minimize:
    std::cout << "Minimize" << std::endl;
    break;
  default:
    BOOST_ASSERT(false);
  }
  std::cout << "    "
            << p_objective->render() 
            << std::endl << std::endl;

  std::cout << "Subject To" << std::endl;
  for (std::vector<ConstraintPtr>:: iterator i = p_constraints.begin();
       i != p_constraints.end(); ++i) {
    std::cout << " " << (*i)->name() << ": " 
              <<(*i)->render() << std::endl;
  }
  std::cout << std::endl;   


  std::cout << "Bounds" << std::endl;
  {
    LPFileVarBoundsLister v(std::cout);
    for_each(p_variables.begin(), p_variables.end(),
             boost::bind(&Variable::accept, _1, boost::ref(v)));

  }
  std::cout << std::endl;

  std::cout << "General" << std::endl;
  std::cout << "    ";
  {
    LPFileGenVarLister v(std::cout);
    for_each(p_variables.begin(), p_variables.end(),
             boost::bind(&Variable::accept, _1, boost::ref(v)));

  }
  std::cout << std::endl << std::endl;

  std::cout << "Integer" << std::endl;
  std::cout << "    ";
  {
    LPFileIntVarLister v(std::cout);
    for_each(p_variables.begin(), p_variables.end(),
             boost::bind(&Variable::accept, _1, boost::ref(v)));

  }
  std::cout << std::endl << std::endl;

  std::cout << "Binary" << std::endl;
  std::cout << "    ";
  {
    LPFileBinVarLister v(std::cout);
    for_each(p_variables.begin(), p_variables.end(),
             boost::bind(&Variable::accept, _1, boost::ref(v)));

  }
  std::cout << std::endl << std::endl;


  std::cout << "End" << std::endl;
}
    
  


} // namespace optimization
} // namespace gridpack

