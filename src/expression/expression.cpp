// -------------------------------------------------------------
// file: expression.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created August 28, 2015 by William A. Perkins
// Last Change: 2015-10-07 13:11:40 d3g096
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

BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::ConstantExpression<int>);
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::ConstantExpression<double>);

BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::VariableExpression);

BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::UnaryMinus);
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::UnaryPlus);

BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Multiplication);
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Division);
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Addition);
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Subtraction);
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Exponentiation);

BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::LessThan);
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::LessThanOrEqual);
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::GreaterThan);
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::GreaterThanOrEqual);
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Equal);

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
    utility::Named()
{
  Named::name(boost::str(boost::format("C%d") % p_nextID++));
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
void ExpressionVisitor::visit(Division& e) { return; };
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

