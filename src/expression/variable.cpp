// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   variable.cpp
 * @author William A. Perkins
 * @date   2015-08-28 15:01:59 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#include <climits>
#include <boost/format.hpp>
#include "variable.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class Variable
// -------------------------------------------------------------

int Variable::p_nextID(0);

// -------------------------------------------------------------
// Variable:: constructors / destructor
// -------------------------------------------------------------
Variable::Variable()
    : utility::Named(),
      utility::Uncopyable(),
      p_id(p_nextID++)
{
  Named::name(boost::str(boost::format("V%d") % p_id));
}

Variable::~Variable(void)
{
}

// -------------------------------------------------------------
// template class BoundedVariableT<>
// -------------------------------------------------------------

template <> const double BoundedVariableT<double>::veryLowValue(-1.0E100);
template <> const double BoundedVariableT<double>::veryHighValue(1.0E100);

template <> const int BoundedVariableT<int>::veryLowValue(INT_MIN);
template <> const int BoundedVariableT<int>::veryHighValue(INT_MAX);


// -------------------------------------------------------------
//  class VariableVisitor
// -------------------------------------------------------------

// -------------------------------------------------------------
// VariableVisitor:: constructors / destructor
// -------------------------------------------------------------
VariableVisitor::VariableVisitor(void)
{
  // empty
}

VariableVisitor::~VariableVisitor(void)
{
  // empty
}

// -------------------------------------------------------------
// VariableVisitor::visit
// -------------------------------------------------------------
void
VariableVisitor::visit(const Variable& var)
{
  // do nothing
}

void
VariableVisitor::visit(const RealVariable& var)
{
  this->visit(static_cast<const Variable&>(var));
}

void
VariableVisitor::visit(const IntegerVariable& var)
{
  this->visit(static_cast<const Variable&>(var));
}

void
VariableVisitor::visit(const BinaryVariable& var)
{
  this->visit(static_cast<const IntegerVariable&>(var));
}


} // namespace optimization
} // namespace gridpack

