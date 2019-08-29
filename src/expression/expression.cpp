// -------------------------------------------------------------
// file: expression.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 28, 2015 by William A. Perkins
// Last Change: 2016-12-20 07:20:15 d3g096
// -------------------------------------------------------------

#include <boost/assert.hpp>
#include <boost/format.hpp>
#include "expression.hpp"

// Cannot do this because these classes have a boost::shared_ptr member
// #include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/weak_ptr.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::ConstantExpression<int>)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::ConstantExpression<double>)

BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::VariableExpression)

BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::UnaryMinus)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::UnaryPlus)

BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Multiplication)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Division)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Addition)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Subtraction)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Exponentiation)

BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Constraint)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::LessThan)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::LessThanOrEqual)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::GreaterThan)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::GreaterThanOrEqual)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Equal)

namespace gridpack {
namespace optimization {

template <typename T>
std::string
ConstantExpression<T>::p_render(void) const
{
  BOOST_ASSERT(false);
  return "";
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
    utility::Named(boost::str(boost::format("C%d") % p_nextID++))
{
  // Constraints can have a null LHS
  // BOOST_ASSERT_MSG(lhs, "Constraint: LHS null");
  BOOST_ASSERT_MSG(rhs, "Constraint: RHS null");
  // BOOST_ASSERT_MSG(!lhs->null(), "Constraint: LHS null contents");
  BOOST_ASSERT_MSG(!rhs->null(), "Constraint: RHS null contests");
}

Constraint::Constraint(void)
  : BinaryExpression(),
    utility::Named()
{}

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
// -------------------------------------------------------------

// For constants and variable do nothing
void ExpressionVisitor::visit(IntegerConstant& e)  { return; }
void ExpressionVisitor::visit(RealConstant& e) { return; }
void ExpressionVisitor::visit(VariableExpression& e) { return; }

// Send the visitor to the rhs
void ExpressionVisitor::visit(UnaryExpression& e)
{
  e.rhs()->accept(*this);
}
void ExpressionVisitor::visit(UnaryMinus& e)
{
  visit(static_cast<UnaryExpression&>(e));
}
void ExpressionVisitor::visit(UnaryPlus& e)
{
  visit(static_cast<UnaryExpression&>(e));
}

// Send the visitor to the rhs and the rhs
void ExpressionVisitor::visit(BinaryExpression& e) 
{
  e.lhs()->accept(*this);
  e.rhs()->accept(*this);
}
void ExpressionVisitor::visit(Multiplication& e)
{
  visit(static_cast<BinaryExpression&>(e));
}
void ExpressionVisitor::visit(Division& e)
{
  visit(static_cast<BinaryExpression&>(e));
}
void ExpressionVisitor::visit(Addition& e)
{
  visit(static_cast<BinaryExpression&>(e));
}
void ExpressionVisitor::visit(Subtraction& e)
{
  visit(static_cast<BinaryExpression&>(e));
}
void ExpressionVisitor::visit(Exponentiation& e)
{
  visit(static_cast<BinaryExpression&>(e));
}

// Send the visitor to the rhs and the rhs
void ExpressionVisitor::visit(Constraint& e)
{
  e.lhs()->accept(*this);
  e.rhs()->accept(*this);
}
void ExpressionVisitor::visit(LessThan& e)
{
  visit(static_cast<Constraint&>(e));
}
void ExpressionVisitor::visit(LessThanOrEqual& e)
{
  visit(static_cast<Constraint&>(e));
}
void ExpressionVisitor::visit(GreaterThan& e)
{
  visit(static_cast<Constraint&>(e));
}
void ExpressionVisitor::visit(GreaterThanOrEqual& e)
{
  visit(static_cast<Constraint&>(e));
}
void ExpressionVisitor::visit(Equal& e)
{
  visit(static_cast<Constraint&>(e));
}


} // namespace optimization
} // namespace gridpack

