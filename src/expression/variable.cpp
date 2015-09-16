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
 * @date   2015-09-16 12:16:22 d3g096
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
VariableVisitor::visit(Variable& var)
{
  // do nothing
}

void
VariableVisitor::visit(RealVariable& var)
{
  this->visit(static_cast<Variable&>(var));
}

void
VariableVisitor::visit(IntegerVariable& var)
{
  this->visit(static_cast<Variable&>(var));
}

void
VariableVisitor::visit(BinaryVariable& var)
{
  this->visit(static_cast<IntegerVariable&>(var));
}


// -------------------------------------------------------------
//  class VariableTable
// -------------------------------------------------------------

// -------------------------------------------------------------
// VariableTable:: constructors / destructor
// -------------------------------------------------------------
VariableTable::VariableTable(std::ostream& out)
  : VariableVisitor(), p_out(out), p_first(true)
{
  
}

VariableTable::~VariableTable(void)
{
}

// -------------------------------------------------------------
// VariableTable::visit
// -------------------------------------------------------------
void 
VariableTable::visit(RealVariable& var)
{
  if (p_first) {
    p_header();
    p_first = false;
  }
  char fmt[] = "%-10.10s (%c) %10.4g";
  p_out << boost::str(boost::format(fmt) % var.name() % 'R' % var.initial())
        << std::endl;
}
 
void 
VariableTable::visit(IntegerVariable& var)
{
  if (p_first) {
    p_header();
    p_first = false;
  }
  char fmt[] = "%-10.10s (%c) %10d";
  p_out << boost::str(boost::format(fmt) % var.name() % 'I' % var.initial())
        << std::endl;
}

void 
VariableTable::visit(BinaryVariable& var)
{
  if (p_first) {
    p_header();
    p_first = false;
  }
  char fmt[] = "%-10.10s (%c) %10d";
  p_out << boost::str(boost::format(fmt) % var.name() % 'B' % var.initial())
        << std::endl;
}


static const std::string bar(25, '-');

// -------------------------------------------------------------
// VariableTable::p_header
// -------------------------------------------------------------
void
VariableTable::p_header(void) const
{
  p_out << bar << std::endl
        << "Variable                 " << std::endl
        << "Name      Type      Value" << std::endl
        << bar << std::endl;
}


} // namespace optimization
} // namespace gridpack

