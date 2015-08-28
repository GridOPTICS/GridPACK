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
 * @file   optimizer.hpp
 * @author William A. Perkins
 * @date   2015-08-28 15:47:25 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _optimizer_hpp_
#define _optimizer_hpp_

#include <vector>
#include <boost/scoped_ptr.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/expression/expression.hpp>

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class OptimizerInterface
// -------------------------------------------------------------
class OptimizerInterface 
  : private utility::Uncopyable
{
public:

  /// Default constructor.
  OptimizerInterface(void)
  {}

  /// Destructor
  virtual ~OptimizerInterface(void)
  {}

  /// Add a (local) variable to be optimized
  void addVariable(VariablePtr v)
  {
    this->p_addVariable(v);
  }

  /// Add a (local) constraint 
  void addConstraint(ConstraintPtr c)
  { 
    this->p_addConstraint(c);
  }

  /// Add the local part of a global constraint (added to other parts)
  void addToGlobalConstraint(const std::string& name, ExpressionPtr expr)
  {
    this->p_addToGlobalConstraint(name, expr);
  }

  /// Add to the local part of the global objective function (added to other parts)
  void addToObjective(ExpressionPtr expr)
  {
    this->p_addToObjective(expr);
  }

  /// Do the problem
  void solve(void)
  {
    this->p_solve();
  }

protected:
  
  /// Add a (local) variable to be optimized (specialized)
  virtual void p_addVariable(VariablePtr v) = 0;

  /// Add a (local) constraint (specialized)
  virtual void p_addConstraint(ConstraintPtr c) = 0;

  /// Add the local part of a global constraint (specialized)
  virtual void p_addToGlobalConstraint(const std::string& name, ExpressionPtr expr) = 0;

  /// Add to the local part of the global objective function (added to other parts)
  virtual void p_addToObjective(ExpressionPtr expr) = 0;

  /// Do the problem (specialized)
  virtual void p_solve(void) = 0;
  
};

// -------------------------------------------------------------
//  class OptimizerImplementation
// -------------------------------------------------------------
class OptimizerImplementation 
  : public OptimizerInterface
{
public:

  /// Default constructor.
  OptimizerImplementation(void)
  {}

  /// Destructor
  ~OptimizerImplementation(void)
  {}

protected:

  /// The variables involved
  std::vector<VariablePtr> p_variables;

  /// The collection of constraints involved
  std::vector<ConstraintPtr> p_constraints;

  /// The objective fuction
  ExpressionPtr p_objective;

  /// Add a (local) variable to be optimized (specialized)
  void p_addVariable(VariablePtr v)
  {
    p_variables.push_back(v);
  }

  /// Add a (local) constraint (specialized)
  void p_addConstraint(ConstraintPtr c)
  {
    p_constraints.push_back(c);
  }

  /// Add the local part of a global constraint (specialized)
  void p_addToGlobalConstraint(const std::string& name, ExpressionPtr expr)
  {
  }    

  /// Add to the local part of the global objective function (added to other parts)
  void p_addToObjective(ExpressionPtr expr)
  {
    if (p_objective) {
      p_objective = p_objective + expr;
    } else {
      p_objective = expr;
    }
  }
};

// -------------------------------------------------------------
//  class Optimizer
// -------------------------------------------------------------
class Optimizer 
  : public OptimizerInterface
{
public:

  /// Default constructor.
  Optimizer(void);

  /// Destructor
  ~Optimizer(void);

protected:
  
  /// Where the work happens
  boost::scoped_ptr<OptimizerImplementation> p_impl;

  /// Add a (local) variable to be optimized (specialized)
  void p_addVariable(VariablePtr v)
  {
    p_impl->addVariable(v);
  }

  /// Add a (local) constraint (specialized)
  void p_addConstraint(ConstraintPtr c)
  {
    p_impl->addConstraint(c);
  }

  /// Add the local part of a global constraint (specialized)
  void p_addToGlobalConstraint(const std::string& name, ExpressionPtr expr)
  {
    p_impl->addToGlobalConstraint(name, expr);
  }    

  /// Add to the local part of the global objective function (added to other parts)
  void p_addToObjective(ExpressionPtr expr)
  {
    p_impl->addToObjective(expr);
  }

  /// Solve the problem
  void p_solve(void)
  {
    p_impl->solve();
  }
};


} // namespace optimization
} // namespace gridpack


#endif
