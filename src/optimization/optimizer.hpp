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
 * @date   2015-10-13 12:31:02 d3g096
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
#include <gridpack/utilities/exception.hpp>
#include <gridpack/configuration/configurable.hpp>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/optimization/expression.hpp>

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

  /// Maximize the objective
  void maximize(void)
  {
    this->p_solve(Maximize);
  }

  /// Minimize the objective
  void minimize(void)
  {
    this->p_solve(Minimize);
  }

protected:

  /// Ways in which the objective can be optimized
  enum p_optimizeMethod { Maximize, Minimize };
  
  /// Add a (local) variable to be optimized (specialized)
  virtual void p_addVariable(VariablePtr v) = 0;

  /// Add a (local) constraint (specialized)
  virtual void p_addConstraint(ConstraintPtr c) = 0;

  /// Add the local part of a global constraint (specialized)
  virtual void p_addToGlobalConstraint(const std::string& name, ExpressionPtr expr) = 0;

  /// Add to the local part of the global objective function (added to other parts)
  virtual void p_addToObjective(ExpressionPtr expr) = 0;

  /// Do the problem (specialized)
  virtual void p_solve(const p_optimizeMethod& m) = 0;
  
};

// -------------------------------------------------------------
//  class OptimizerImplementation
// -------------------------------------------------------------
class OptimizerImplementation 
  : public OptimizerInterface,
    public parallel::Distributed,
    public utility::Configurable
{
public:

  /// A map of variable name to variable
  typedef std::map<std::string, VariablePtr> VarMap;

  /// Default constructor.
  OptimizerImplementation(const parallel::Communicator& comm)
    : OptimizerInterface(),
      parallel::Distributed(comm),
      utility::Configurable("Optimizer")
  {}

  /// Destructor
  ~OptimizerImplementation(void)
  {}

protected:

  /// The (local) variables involved
  std::vector<VariablePtr> p_variables;

  /// The collection of (local) constraints involved
  std::vector<ConstraintPtr> p_constraints;

  /// The (local part of the) objective fuction
  ExpressionPtr p_objective;

  /// Variables involved from all processes (unique instances)
  VarMap p_allVariables;

  /// Constraints involved from all processes 
  std::vector<ConstraintPtr> p_allConstraints;

  /// The entire objective function
  ExpressionPtr p_fullObjective;

  /// The global constraint expressions (local parts)
  std::map<std::string, ExpressionPtr> p_globalConstraints;

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
    if (p_globalConstraints[name]) {
      p_globalConstraints[name] = p_globalConstraints[name] + expr;
    } else {
      p_globalConstraints[name] = expr;
    }
  }  

  /// Get the global constraint expression
  ExpressionPtr p_getGlobalConstraint(const std::string& name)
  {
    ExpressionPtr result;
    try {
      result = p_globalConstraints.at(name);
    } catch (const std::out_of_range& e) {
      std::string msg(name);
    }
    return result;
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

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props)
  {
    if (props) {
    }
  }

  /// Gather the problem to all processors
  void p_gatherProblem(void);
};

// -------------------------------------------------------------
//  class Optimizer
// -------------------------------------------------------------
class Optimizer 
  : public OptimizerInterface,
    public parallel::WrappedDistributed,
    public utility::WrappedConfigurable
{
public:

  /// Default constructor.
  Optimizer(const parallel::Communicator& comm);

  /// Destructor
  ~Optimizer(void);

protected:
  
  /// Where the work happens
  boost::scoped_ptr<OptimizerImplementation> p_impl;

  /// Set the implementation
  /** 
   * Does what is necessary to set the \ref
   * OptimizerImplementation "implementation".  Subclasses are
   * required to call this at construction.
   * 
   * @param impl specific nonlinear solver implementation to use
   */
  void p_setImpl(OptimizerImplementation *impl)
  {
    p_impl.reset(impl);
    p_setDistributed(p_impl.get());
    p_setConfigurable(p_impl.get());
  }
  

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
  void p_solve(const p_optimizeMethod& m)
  {
    switch (m) {
    case (Maximize):
      p_impl->maximize();
      break;
    case (Minimize):
      p_impl->minimize();
      break;
    default:
      BOOST_ASSERT(false);
    }
  }
};


} // namespace optimization
} // namespace gridpack


#endif
