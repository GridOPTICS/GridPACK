// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   constraint_renderer.hpp
 * @author William A. Perkins
 * @date   2017-03-22 09:15:48 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December  8, 2016 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _constraint_renderer_hpp_
#define _constraint_renderer_hpp_

#include <boost/assert.hpp>
#include <boost/format.hpp>

#include "gridpack/expression/expression.hpp"
#include "gridpack/expression/functions.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class ConstraintRenderer
// -------------------------------------------------------------
class ConstraintRenderer 
  : public ExpressionVisitor
{
public:

  /// Default constructor.
  ConstraintRenderer(std::ostream& out)
    : ExpressionVisitor(), p_out(out)
  {}

  /// Destructor
  ~ConstraintRenderer(void)
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
    this->visit(static_cast<BinaryExpression&>(e));
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
      BOOST_ASSERT_MSG(false, "ConstraintRenderer: Invalid RHS to Division");
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

  void visit(Function& f)
  {
    bool first(true);
    p_out << f.name(); 
    p_out << "(";
    BOOST_FOREACH(const ExpressionPtr a, f.args()) {
      if (!first) {
        p_out << ", ";
      }
      first = false;
      a->accept(*this);
    }
    p_out << ")";
  }

protected:

  /// The stream to send renderings
  std::ostream& p_out;

  /// How to group an expression with higher precedence
  virtual void p_group(ExpressionPtr e)
  {
    p_out << "(";
    e->accept(*this);
    p_out << ")";
  }

  
};



} // namespace optimization
} // namespace gridpack

#endif
