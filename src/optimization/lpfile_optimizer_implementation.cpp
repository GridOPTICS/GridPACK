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
 * @date   2015-11-03 13:10:41 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace io = boost::iostreams;

#include "gridpack/utilities/exception.hpp"
#include "lpfile_optimizer_implementation.hpp"

// -------------------------------------------------------------
// line_wrapping_output_filter
//
// slightly modified from the Boost iostreams example
// -------------------------------------------------------------
class line_wrapping_output_filter : public io::output_filter {
public:
  explicit line_wrapping_output_filter(int line_length = 80, int margin = 8)
    : do_new_line_(false), line_length_(line_length), end_margin_(margin), col_no_(0) 
  { }

  template<typename Sink>
  bool put(Sink& dest, int c)
  {
    // if the column is past the goal end of line, find the next white space, eat it up, then
    if (do_new_line_) {
      if (c == '\n') {
        do_new_line_ = false;
      } else if (!std::isspace(c)) {
        if (!put_char(dest, '\n')) return false;
        do_new_line_ = false;
      } else {
        // eating white space
        return true;
      }
    } else if (col_no_ >= (line_length_ - end_margin_)) {
      // wait until we see some white space before doing a new line
      if (std::isspace(c)) {
        do_new_line_ = true;
      }
    }
    return put_char(dest, c);
  }

  template<typename Sink>
  void close(Sink&)
  { col_no_ = 0; }

private:
  template<typename Sink>
  bool put_char(Sink& dest, int c)
  {
    if (!io::put(dest, c))
      return false;
    if (c != '\n')
      ++col_no_;
    else
      col_no_ = 0;
    return true;
  }
  int  do_new_line_;
  int  line_length_;
  int  end_margin_;
  int  col_no_;
};


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
  : public ExpressionVisitor
{
public:

  /// Default constructor.
  LPFileConstraintRenderer(std::ostream& out)
    : ExpressionVisitor(), p_out(out)
  {}

  /// Destructor
  ~LPFileConstraintRenderer(void)
  {}

  void visit(IntegerConstant& e)
  { 
    p_out << boost::str(boost::format("%d") % e.value());
  }

  void visit(RealConstant& e)
  { 
    p_out << boost::str(boost::format("%g") % e.value());
  }

  void visit(VariableExpression& e)
  { 
    p_out << e.name();
  }

  void visit(UnaryExpression& e)
  {
    // should not be called
    BOOST_ASSERT(false);
  }

  void visit(UnaryMinus& e)
  {
    p_out << "-";
    if (e.rhs()->precedence() > e.precedence()) {
      p_group(e.rhs());
    } else {
      e.rhs()->accept(*this);
    }
  }

  void visit(UnaryPlus& e)
  {
    p_out << "+";
    if (e.rhs()->precedence() > e.precedence()) {
      p_group(e.rhs());
    } else {
      e.rhs()->accept(*this);
    }
  }

  void visit(BinaryExpression& e)
  {
    ExpressionPtr lhs(e.lhs());
    ExpressionPtr rhs(e.rhs());

    if (lhs->precedence() > e.precedence()) {
      p_group(lhs);
    } else {
      lhs->accept(*this);
    }
    p_out << " " << e.op() << " ";
    if (rhs->precedence() > e.precedence()) {
      p_group(rhs);
    } else {
      rhs->accept(*this);
    }
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

  void visit(Division& e)
  {
    ExpressionPtr lhs(e.lhs());
    ExpressionChecker lcheck;
    lhs->accept(lcheck);

    ExpressionPtr rhs(e.rhs());
    ExpressionChecker rcheck;
    rhs->accept(rcheck);

    // check to see that the RHS is a constant, otherwise it's a
    // nonlinear expression that is not handled by the LP language
    if (!lcheck.isConstant) { 
      BOOST_ASSERT_MSG(false, "LPFileConstraintRenderer: Invalid RHS to Division");
    }

    // consider switching the LHS and RHS if the RHS is constant; for
    // now, the user should write correctly.

    if (lhs->precedence() > e.precedence()) {
      p_group(lhs);
    } else {
      lhs->accept(*this);
    }
    p_out << "/";
    rhs->accept(*this);         // it's a constant, right?
  }  

  void visit(Addition& e)
  {
    this->visit(static_cast<BinaryExpression&>(e));
  }

  void visit(Subtraction& e)
  {
    ExpressionPtr lhs(e.lhs());
    ExpressionPtr rhs(e.rhs());

    if (lhs->precedence() > e.precedence()) {
      p_group(lhs);
    } else {
      lhs->accept(*this);
    }
    p_out << " - ";
    if (rhs->precedence() >= e.precedence()) {
      p_group(rhs);
    } else {
      rhs->accept(*this);
    }
  }
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

    if (lhs->precedence() > e.precedence()) {
      p_group(lhs);
    } else {
      lhs->accept(*this);
    }
    p_out << " ^ ";
    rhs->accept(*this);         // it's a constant, right?
  }

  void visit(Constraint& e)
  {
    ExpressionPtr lhs(e.lhs());
    ExpressionPtr rhs(e.rhs());

    p_out << e.name() << ": ";
    if (lhs->precedence() > e.precedence()) {
      p_group(lhs);
    } else {
      lhs->accept(*this);
    }
    p_out << " " << e.op() << " ";
    if (rhs->precedence() > e.precedence()) {
      p_group(rhs);
    } else {
      rhs->accept(*this);
    }
    p_out << std::endl;
  }

  void visit(LessThan& e)
  {
    this->visit(static_cast<Constraint&>(e));
  }

  void visit(LessThanOrEqual& e)
  {
    this->visit(static_cast<Constraint&>(e));
  }

  void visit(GreaterThan& e)
  {
    this->visit(static_cast<Constraint&>(e));
  }

  void visit(GreaterThanOrEqual& e)
  {
    this->visit(static_cast<Constraint&>(e));
  }

  void visit(Equal& e)
  {
    this->visit(static_cast<Constraint&>(e));
  }
  

protected:

  /// The stream to send renderings
  std::ostream& p_out;

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
// LPFileOptimizerImplementation::p_temporaryFile
// -------------------------------------------------------------
std::string
LPFileOptimizerImplementation::p_temporaryFileName(void)
{
  using namespace boost::filesystem;
  std::string pattern = 
    boost::str(boost::format("%s_%d.lp") % "gridpack%%%%" % this->processor_rank());
  path model(pattern);
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


// -------------------------------------------------------------
// LPFileOptimizerImplementation::p_solve
// -------------------------------------------------------------
void
LPFileOptimizerImplementation::p_solve(const p_optimizeMethod& method)
{
  std::string tmpname(p_temporaryFileName());
  std::ofstream tmp;
  tmp.open(tmpname.c_str());
  if (!tmp) {
    std::string msg("Cannot open temporary file: ");
    msg += tmpname.c_str();
    throw gridpack::Exception(msg);
  }
  p_write(method, tmp);
  tmp.close();
}
    
  


} // namespace optimization
} // namespace gridpack

