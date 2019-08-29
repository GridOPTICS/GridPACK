// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   julia_optimizer_implementation.cpp
 * @author William A. Perkins
 * @date   2016-12-08 14:58:51 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December  7, 2016 by William A. Perkins
// Last Change: 2016-12-07 14:35:30 d3g096
// -------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/bind.hpp>

#include "gridpack/utilities/exception.hpp"
#include "line_wrapping_output_filter.hpp"
#include "constraint_renderer.hpp"
#include "julia_optimizer_implementation.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class JuliaVarLister
// -------------------------------------------------------------
class JuliaVarLister
  : public VariableVisitor
{
public:

  /// Default constructor.
  JuliaVarLister(const std::string& mname, std::ostream& o)
    : VariableVisitor(),
      p_model(mname),
      p_stream(o)
  {
  }

  /// Destructor
  ~JuliaVarLister(void)
  {
  }

  void visit(Variable& var)
  {
    BOOST_ASSERT(false);
  }

  void visit(RealVariable& var)
  {
    if (var.getNoInit()) return;
    std::string s;
    s = var.name();
    if (var.bounded()) {
      if (var.lowerBound() > var.veryLowValue) {
        s += ", lowerbound = ";
        s += boost::str(boost::format("%8.4g") % var.lowerBound());
      } 
      if (var.upperBound() < var.veryHighValue) {
        s += ", upperbound = ";
        s += boost::str(boost::format("%8.4g") % var.upperBound());
      } 
    }
    this->p_stream << "@variable(" << p_model << ", " << s << ")" << std::endl;
  }
    
  void visit(IntegerVariable& var)
  {
    if (var.getNoInit()) return;
    std::string s;
    s = var.name();
    if (var.bounded()) {
      if (var.lowerBound() > var.veryLowValue) {
        s += ", lowerbound = ";
        s += boost::str(boost::format("%8d") % var.lowerBound());
      } 
      if (var.upperBound() < var.veryHighValue) {
        s += ", upperbound = ";
        s += boost::str(boost::format("%8d") % var.upperBound());
      } 
    }
    this->p_stream << "@variable(" << p_model << ", " << s << ", Int)" << std::endl;
  }

  void visit(BinaryVariable& var)
  {
    if (var.getNoInit()) return;
    this->p_stream << "@variable(" << p_model << ", " << var.name() << ", Bin)" << std::endl;
  }

protected:

  std::string p_model;
  std::ostream& p_stream;

};

// -------------------------------------------------------------
//  class JuliaVarValueLister
// -------------------------------------------------------------
class JuliaVarValueLister
  : public VariableVisitor
{
public:

  /// Default constructor.
  JuliaVarValueLister(const std::string& mname, std::ostream& o)
    : VariableVisitor(),
      p_model(mname),
      p_stream(o)
  {
  }

  /// Destructor
  ~JuliaVarValueLister(void)
  {
  }

  void visit(Variable& var)
  {
    BOOST_ASSERT(false);
  }

  void visit(RealVariable& var)
  {
    if (var.getNoInit()) return;
    std::string s;
    s = var.name();
    s += ", ";
    s += boost::str(boost::format("%8.4g") %var.initial());
    this->p_stream << "setvalue(" << s << ")" << std::endl;
  }
    
  void visit(IntegerVariable& var)
  {
    if (var.getNoInit()) return;
    std::string s;
    s = var.name();
    s += ", ";
    s += boost::str(boost::format("%8d") %var.initial());
    this->p_stream << "setvalue(" << s << ")" << std::endl;
  }

  void visit(BinaryVariable& var)
  {
    if (var.getNoInit()) return;
    std::string s;
    s = var.name();
    s += ", ";
    s += boost::str(boost::format("%8d") %var.initial());
    this->p_stream << "setvalue(" << s << ")" << std::endl;
  }

protected:

  std::string p_model;
  std::ostream& p_stream;

};

// -------------------------------------------------------------
//  class JuliaVarPrintLister
// -------------------------------------------------------------
class JuliaVarPrintLister
  : public VariableVisitor
{
public:

  /// Default constructor.
  JuliaVarPrintLister(const std::string& mname, std::ostream& o)
    : VariableVisitor(),
      p_model(mname),
      p_stream(o)
  {
  }

  /// Destructor
  ~JuliaVarPrintLister(void)
  {
  }

  void visit(Variable& var)
  {
    BOOST_ASSERT(false);
  }

  void visit(RealVariable& var)
  {
    if (var.getNoInit()) return;
    std::string s;
    s = var.name();
    this->p_stream << "println(\""<<s<<" value: \",getvalue(" << s << "))" << std::endl;
  }
    
  void visit(IntegerVariable& var)
  {
    if (var.getNoInit()) return;
    std::string s;
    s = var.name();
    this->p_stream << "println(\""<<s<<" value: \",getvalue(" << s << "))" << std::endl;
  }

  void visit(BinaryVariable& var)
  {
    if (var.getNoInit()) return;
    std::string s;
    this->p_stream << "println(\""<<s<<" value: \",getvalue(" << s << "))" << std::endl;
  }

protected:

  std::string p_model;
  std::ostream& p_stream;

};

// -------------------------------------------------------------
//  class JuliaConstraintRenderer
// -------------------------------------------------------------
class JuliaConstraintRenderer 
  : public ConstraintRenderer
{
public:

  /// Default constructor.
  JuliaConstraintRenderer(const std::string& mname, std::ostream& out)
    : ConstraintRenderer(out), p_model(mname)
  {
  }

  /// Destructor
  ~JuliaConstraintRenderer(void)
  {}

  void visit(Constraint& e)
  {
    p_out << "@NLconstraint(" << p_model << ", ";
    ConstraintRenderer::visit(e);
    p_out << ")" << std::endl;
  }

  void visit(Equal& e)
  {
    ExpressionPtr lhs(e.lhs());
    ExpressionPtr rhs(e.rhs());

    p_out << "@NLconstraint(" << p_model << ", ";
    if (lhs->precedence() > e.precedence()) {
      p_group(lhs);
    } else {
      lhs->accept(*this);
    }
    p_out << " " << "==" << " ";
    if (rhs->precedence() > e.precedence()) {
      p_group(rhs);
    } else {
      rhs->accept(*this);
    }
    p_out << ")" << std::endl;
  }

  void visit(Exponentiation& e)
  {
    ExpressionPtr lhs(e.lhs());
    ExpressionPtr rhs(e.rhs());
    ExpressionChecker rcheck;
    rhs->accept(rcheck);

    if (lhs->precedence() > e.precedence()) {
      p_group(lhs);
    } else {
      lhs->accept(*this);
    }
    p_out << "^";
    rhs->accept(*this);         // it's a constant, right?
  }

protected:

  /// Name of the Model
  std::string p_model;
};

// -------------------------------------------------------------
//  class JuliaOptimizerImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// JuliaOptimizerImplementation:: constructors / destructor
// -------------------------------------------------------------
JuliaOptimizerImplementation::JuliaOptimizerImplementation(const parallel::Communicator& comm)
  : FileOptimizerImplementation(comm)
{
  
}

JuliaOptimizerImplementation::~JuliaOptimizerImplementation(void)
{
}

// -------------------------------------------------------------
// JuliaOptimizerImplementation::p_configure
// -------------------------------------------------------------
void
JuliaOptimizerImplementation::p_configure(utility::Configuration::CursorPtr props)
{
  FileOptimizerImplementation::p_configure(props);
  p_outputName += ".jl";
}

// -------------------------------------------------------------
// JuliaOptimizerImplementation::p_setFilename
// -------------------------------------------------------------
void JuliaOptimizerImplementation::p_setFilename(std::string file)
{
  FileOptimizerImplementation::p_setFilename(file);
  p_outputName += ".jl";
}

// -------------------------------------------------------------
// JuliaOptimizerImplementation::p_write
// -------------------------------------------------------------
void 
JuliaOptimizerImplementation::p_write(const p_optimizeMethod& m, std::ostream& output)
{
  p_gatherProblem();

  io::filtering_stream<io::output> out;
  out.push(line_wrapping_output_filter());
  out.push(output);

  std::string mname("gpm");
  std::string solver;

  out << "using JuMP" << std::endl;

  // FIXME: configurable
  // out << "using GLPKMathProgInterface" << std::endl;
  // solver = "GLPKSolverLP()");

  out << "using Ipopt" << std::endl;
  solver = "IpoptSolver()";

  out << mname << " = Model(solver=" << solver << ")" << std::endl;
  // List all variables
  { 
    JuliaVarLister v(mname, out);
    BOOST_FOREACH(VarMap::value_type& i, p_exportVariables) {
      i.second->accept(v);
    }
  }
  // List initial value of all variables
  { 
    JuliaVarValueLister v(mname, out);
    BOOST_FOREACH(VarMap::value_type& i, p_exportVariables) {
      i.second->accept(v);
    }
  }
  {
    JuliaConstraintRenderer r(mname, out);
    BOOST_FOREACH(ConstraintPtr c, p_allConstraints) {
      c->accept(r);
    }
  }
  out << "@objective(" << mname << ", ";
  switch (m) {
  case Maximize:
    out << ":Max";
    break;
  case Minimize:
    out << ":Min";
    break;
  default:
    BOOST_ASSERT(false);
    break;
  }
  out << ", ";
  {
    ConstraintRenderer r(out);
    if (p_fullObjective) {
      p_fullObjective->accept(r);
    }
  }
  out << ")" << std::endl;

  out << "print(" << mname << ")" << std::endl;
  out << "status = solve(" << mname << ")" << std::endl;
  out << "println(\"Objective value: \", getobjectivevalue(" << mname << "))" << std::endl;

  // Print values of optimized variables
  { 
    JuliaVarPrintLister v(mname, out);
    BOOST_FOREACH(VarMap::value_type& i, p_exportVariables) {
      i.second->accept(v);
    }
  }

}




} // namespace optimization
} // namespace gridpack
