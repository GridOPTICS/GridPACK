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
 * @date   2015-09-16 11:31:17 d3g096
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
#include <boost/filesystem.hpp>
#include "gridpack/utilities/exception.hpp"
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

  void visit(Variable& var)
  {
    BOOST_ASSERT(false);
  }

  void visit(RealVariable& var)
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
    
  void visit(IntegerVariable& var)
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

  void visit(RealVariable& var)
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

  void visit(IntegerVariable& var)
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

  void visit(BinaryVariable& var)
  {
    this->p_stream << " " << var.name();
  };
};

// -------------------------------------------------------------
//  class LPFileOptimizerImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// LPFileOptimizerImplementation::p_temporaryFile
// -------------------------------------------------------------
std::string
LPFileOptimizerImplementation::p_temporaryFileName(void)
{
  using namespace boost::filesystem;
  path model("gridpack%%%%.lp");
  path tmp(temp_directory_path());
  tmp /= unique_path(model);

  boost::system::error_code ec;
  file_status istat = status(tmp);
  std::string result(tmp.c_str());
  return result;
  
}

// -------------------------------------------------------------
// LPFileOptimizerImplementation::p_write
// -------------------------------------------------------------
void
LPFileOptimizerImplementation::p_write(const p_optimizeMethod& method, std::ostream& out)
{
  out << "\\* Problem name: Test\\*" << std::endl << std::endl;

  switch (method) {
  case Maximize:
    out << "Maximize" << std::endl;
    break;
  case Minimize:
    out << "Minimize" << std::endl;
    break;
  default:
    BOOST_ASSERT(false);
  }
  out << "    "
            << p_objective->render() 
            << std::endl << std::endl;

  out << "Subject To" << std::endl;
  for (std::vector<ConstraintPtr>:: iterator i = p_constraints.begin();
       i != p_constraints.end(); ++i) {
    out << " " << (*i)->name() << ": " 
              <<(*i)->render() << std::endl;
  }
  out << std::endl;   


  out << "Bounds" << std::endl;
  {
    LPFileVarBoundsLister v(out);
    for_each(p_variables.begin(), p_variables.end(),
             boost::bind(&Variable::accept, _1, boost::ref(v)));

  }
  out << std::endl;

  out << "General" << std::endl;
  out << "    ";
  {
    LPFileGenVarLister v(out);
    for_each(p_variables.begin(), p_variables.end(),
             boost::bind(&Variable::accept, _1, boost::ref(v)));

  }
  out << std::endl << std::endl;

  out << "Integer" << std::endl;
  out << "    ";
  {
    LPFileIntVarLister v(out);
    for_each(p_variables.begin(), p_variables.end(),
             boost::bind(&Variable::accept, _1, boost::ref(v)));

  }
  out << std::endl << std::endl;

  out << "Binary" << std::endl;
  out << "    ";
  {
    LPFileBinVarLister v(out);
    for_each(p_variables.begin(), p_variables.end(),
             boost::bind(&Variable::accept, _1, boost::ref(v)));

  }
  out << std::endl << std::endl;


  out << "End" << std::endl;
}


// -------------------------------------------------------------
// LPFileOptimizerImplementation::p_solve
// -------------------------------------------------------------
void
LPFileOptimizerImplementation::p_solve(const p_optimizeMethod& method)
{
  p_write(method, std::cout);
}
    
  


} // namespace optimization
} // namespace gridpack

