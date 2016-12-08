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
 * @date   2016-12-07 15:30:05 d3g096
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
    this->p_stream << "@variable(" << p_model << ", " << var.name() << ", Bin)" << std::endl;
  }

protected:

  std::string p_model;
  std::ostream& p_stream;

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
// JuliaOptimizerImplementation::p_temporaryFileName
// -------------------------------------------------------------
std::string
JuliaOptimizerImplementation::p_temporaryFileName(void)
{
  std::string s(FileOptimizerImplementation::p_temporaryFileName());
  s += ".jl";
  return s;
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

  out << "using JuMP" << std::endl;
  std::string mname("gpm");
  out << mname << " = Model()" << std::endl;
  { 
    JuliaVarLister v(mname, out);
    BOOST_FOREACH(VarMap::value_type& i, p_allVariables) {
      i.second->accept(v);
    }
  }
}




} // namespace optimization
} // namespace gridpack
