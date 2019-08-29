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
 * @date   2016-12-08 14:58:38 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/bind.hpp>

#include "gridpack/utilities/exception.hpp"
#include "line_wrapping_output_filter.hpp"
#include "constraint_renderer.hpp"
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
      s += boost::str(boost::format("%s") % var.name());
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
      s += boost::str(boost::format("%s") % var.name());
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
    s += boost::str(boost::format("%s") % var.name());
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
//  class LPFileConstraintRenderer
// -------------------------------------------------------------
class LPFileConstraintRenderer 
  : public ConstraintRenderer
{
public:

  /// Default constructor.
  LPFileConstraintRenderer(std::ostream& out)
    : ConstraintRenderer(out)
  {}

  /// Destructor
  ~LPFileConstraintRenderer(void)
  {}

  void visit(Exponentiation& e)
  {
    ExpressionPtr lhs(e.lhs());
    ExpressionPtr rhs(e.rhs());
    ExpressionChecker rcheck;
    rhs->accept(rcheck);

    // Only integer constants are allowed as exponents -- check that
    if (!rcheck.isInteger) {
      BOOST_ASSERT_MSG(false, "LPFileConstraintRenderer: Only integer exponents allowed");
    }

    ConstraintRenderer::visit(e);
  }

  void visit(Multiplication& e)
  {
    ExpressionPtr lhs(e.lhs());
    ExpressionChecker lcheck;
    lhs->accept(lcheck);

    ExpressionPtr rhs(e.rhs());
    ExpressionChecker rcheck;
    rhs->accept(rcheck);

    // check to see that the LHS is a constant, otherwise it's a
    // nonlinear expression that is not handled by the LP language
    if (!lcheck.isConstant) { 
      BOOST_ASSERT_MSG(false, "LPFileConstraintRenderer: Invalid LHS to Multiplication");
    }

    if (rcheck.isExponentiation) {
      // SPECIAL CASE: if the RHS is
      // exponentiation, the whole
      // expression needs to be grouped
      p_out << "[";
      lhs->accept(*this);
      p_out << " ";
      rhs->accept(*this);
      p_out << "]";
    } else {
      // consider switching the LHS and RHS
      // if the RHS is constant; for now,
      // the user should write correctly.

      lhs->accept(*this);         // it's a constant, right?
      p_out << " ";
      if (rhs->precedence() > e.precedence()) {
        p_group(rhs);
      } else {
        rhs->accept(*this);
      }
    }
  }

  void visit(Constraint& e)
  {
    p_out << e.name() << ": ";
    ConstraintRenderer::visit(e);
    p_out << std::endl;
  }

protected:

  /// How to group an expression with higher precedence
  void p_group(ExpressionPtr e)
  {
    p_out << "[";
    e->accept(*this);
    p_out << "]";
  }
};

// -------------------------------------------------------------
//  class LPFileOptimizerImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// LPFileOptimizerImplementation::p_configure
// -------------------------------------------------------------
void
LPFileOptimizerImplementation::p_configure(utility::Configuration::CursorPtr props)
{
  FileOptimizerImplementation::p_configure(props);
  p_outputName += ".lp";
}

// -------------------------------------------------------------
// LPFileOptimizerImplementation::p_setFilename
// -------------------------------------------------------------
void LPFileOptimizerImplementation::p_setFilename(std::string file)
{
  FileOptimizerImplementation::p_setFilename(file);
  p_outputName += ".lp";
}

// -------------------------------------------------------------
// LPFileOptimizerImplementation::p_write
// -------------------------------------------------------------
void
LPFileOptimizerImplementation::p_write(const p_optimizeMethod& method, std::ostream& output)
{
  p_gatherProblem();

  io::filtering_stream<io::output> out;
  out.push(line_wrapping_output_filter());
  out.push(output);

  out << "\\* Problem name: GridPACK \\*" << std::endl << std::endl;

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
  {
    LPFileConstraintRenderer r(out);
    p_fullObjective->accept(r);
  }
  out << std::endl << std::endl;

  out << "Subject To" << std::endl;
  {
    LPFileConstraintRenderer r(out);
    BOOST_FOREACH(ConstraintPtr c, p_allConstraints) {
      c->accept(r);
    }
  }
  out << std::endl;   


  out << "Bounds" << std::endl;
  {
    LPFileVarBoundsLister v(out);
    BOOST_FOREACH(VarMap::value_type& i, p_allVariables) {
      i.second->accept(v);
    }
  }
  out << std::endl;

  VariableCounter cnt;
  BOOST_FOREACH(VarMap::value_type& i, p_allVariables) {
    i.second->accept(cnt);
  }
/**
  if (cnt.numReal > 0) {
    out << "General" << std::endl;
    {
      LPFileGenVarLister v(out);
      BOOST_FOREACH(VarMap::value_type& i, p_allVariables) {
        i.second->accept(v);
      }
    }
    out << std::endl << std::endl;
  }
**/
  if (cnt.numInt > 0) {
    out << "Integer" << std::endl;
    {
      LPFileIntVarLister v(out);
      BOOST_FOREACH(VarMap::value_type& i, p_allVariables) {
        i.second->accept(v);
      }
    }
    out << std::endl << std::endl;
  }

  if (cnt.numBin > 0) {
    out << "Binary" << std::endl;
    {
      LPFileBinVarLister v(out);
      BOOST_FOREACH(VarMap::value_type& i, p_allVariables) {
        i.second->accept(v);
      }
    }
    out << std::endl << std::endl;
  }

  out << "End" << std::endl;
}


    
  


} // namespace optimization
} // namespace gridpack

