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
 * @date   2017-02-10 08:18:20 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#include <climits>
#include <boost/format.hpp>
#include "variable.hpp"

#include <boost/mpi.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#if defined(GRIDPACK_AVOID_APPLECLANG_MPI_PROBLEMS)
// For some reason this is required for CLang on Mac for the MPI
// serialization is implemented. 
namespace boost
{
namespace mpi
{
template<> struct is_mpi_datatype<std::string> : public mpl::true_ { };
}
}
#endif

BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::Variable)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::BoundedVariableT<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::BoundedVariableT<int>)
BOOST_CLASS_EXPORT_IMPLEMENT(gridpack::optimization::BinaryVariable)


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
      p_id(p_nextID++),
      p_no_init(false)
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
  char fmt[] = "%-10.10s (%c) %10.4g %10.4g %10.4g";
  p_out << boost::str(boost::format(fmt) % 
                      var.name() % 'R' % var.initial() % 
                      var.lowerBound() % var.upperBound())
        << std::endl;
}
 
void 
VariableTable::visit(IntegerVariable& var)
{
  if (p_first) {
    p_header();
    p_first = false;
  }
  char fmt[] = "%-10.10s (%c) %10d %10d %10d";
  p_out << boost::str(boost::format(fmt) % 
                      var.name() % 'I' % var.initial() %
                      var.lowerBound() % var.upperBound())
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


static const std::string bar(47, '-');

// -------------------------------------------------------------
// VariableTable::p_header
// -------------------------------------------------------------
void
VariableTable::p_header(void) const
{
  p_out << bar << std::endl
        << "Variable                                       " << std::endl
        << "Name      Type      Value      Lower      Upper" << std::endl
        << bar << std::endl;
}

// -------------------------------------------------------------
//  class VariableCounter
// -------------------------------------------------------------

// -------------------------------------------------------------
// VariableCounter:: constructors / destructor
// -------------------------------------------------------------
VariableCounter::VariableCounter()
  : VariableVisitor(), 
    numVar(0), numReal(0), numInt(0), numBin(0)
{
  
}

VariableCounter::~VariableCounter(void)
{
}

// -------------------------------------------------------------
// VariableCounter::visit
// -------------------------------------------------------------
void 
VariableCounter::visit(Variable& var)
{
  BOOST_ASSERT_MSG(false, "VariableCounter: error: variable type not handled");
}

void 
VariableCounter::visit(RealVariable& var)
{
  numVar += 1;
  numReal += 1;
}

void 
VariableCounter::visit(IntegerVariable& var)
{
  numVar += 1;
  numInt += 1;
}

void 
VariableCounter::visit(BinaryVariable& var)
{
  numVar += 1;
  numBin += 1;
}



} // namespace optimization
} // namespace gridpack

