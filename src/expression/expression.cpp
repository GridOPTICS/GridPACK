// -------------------------------------------------------------
// file: expression.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 28, 2015 by William A. Perkins
// Last Change: 2015-09-21 11:46:37 d3g096
// -------------------------------------------------------------

#include <boost/assert.hpp>
#include <boost/format.hpp>
#include "expression.hpp"

namespace gridpack {
namespace optimization {

template <typename T>
std::string
ConstantExpression<T>::p_render(void) const
{
  BOOST_ASSERT(false);
}

template <>
std::string
ConstantExpression<int>::p_render(void) const
{
  return boost::str(boost::format("%d") % p_value);
}

template <>
std::string
ConstantExpression<double>::p_render(void) const
{
  return boost::str(boost::format("%g") % p_value);
}

// -------------------------------------------------------------
// Constraint
// -------------------------------------------------------------
Constraint::Constraint(const int& prec, const std::string& op, 
                       ExpressionPtr lhs, ExpressionPtr rhs)
  : BinaryExpression(prec, op, lhs, rhs), 
    utility::Named()
  {
    Named::name(boost::str(boost::format("C%d") % p_nextID++));
  }

int Constraint::p_nextID(0);


// -------------------------------------------------------------
//  class ExpressionVisitor
// -------------------------------------------------------------

// -------------------------------------------------------------
// ExpressionVisitor:: constructors / destructor
// -------------------------------------------------------------
ExpressionVisitor::ExpressionVisitor()
  : utility::Uncopyable()
{
  
}

ExpressionVisitor::~ExpressionVisitor(void)
{
}

// -------------------------------------------------------------
// ExpressionVisitor::visit
//
// The default behavior is to do nothing. 
// -------------------------------------------------------------

void ExpressionVisitor::visit(IntegerConstant& e)  { return; };
void ExpressionVisitor::visit(RealConstant& e) { return; };
void ExpressionVisitor::visit(VariableExpression& e) { return; };

void ExpressionVisitor::visit(UnaryExpression& e) { return; };
void ExpressionVisitor::visit(UnaryMinus& e) { return; };
void ExpressionVisitor::visit(UnaryPlus& e) { return; };

void ExpressionVisitor::visit(BinaryExpression& e) { return; };
void ExpressionVisitor::visit(Multiplication& e) { return; };
void ExpressionVisitor::visit(Addition& e) { return; };
void ExpressionVisitor::visit(Subtraction& e) { return; };
void ExpressionVisitor::visit(Exponentiation& e) { return; };

void ExpressionVisitor::visit(Constraint& e) { return; };
void ExpressionVisitor::visit(LessThan& e) { return; };
void ExpressionVisitor::visit(LessThanOrEqual& e) { return; };
void ExpressionVisitor::visit(GreaterThan& e) { return; };
void ExpressionVisitor::visit(GreaterThanOrEqual& e) { return; };
void ExpressionVisitor::visit(Equal& e) { return; };


} // namespace optimization
} // namespace gridpack

