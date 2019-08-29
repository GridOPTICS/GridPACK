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
 * @date   2017-03-22 07:53:39 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _optimizer_hpp_
#define _optimizer_hpp_

#include <vector>
#include <map>
#include <boost/scoped_ptr.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/utilities/exception.hpp>
#include <gridpack/configuration/configurable.hpp>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/expression/expression.hpp>
#include <gridpack/expression/functions.hpp>

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

  /// Set filename for output
  virtual void setFilename(std::string file)
  {
    this->p_setFilename(file);
  }

  /// Add a (local) variable to be optimized
  void addVariable(VariablePtr v)
  {
    this->p_addVariable(v);
  }

  /// Add a (local) auxiliary variable. This is used for
  /// parallel applications to create ghost copies of variables
  /// defined on other processors
  void addAuxVariable(VariablePtr v)
  {
    this->p_addAuxVariable(v);
  }

  /// Add a (local) constraint 
  void addConstraint(ConstraintPtr c)
  { 
    this->p_addConstraint(c);
  }

  /// Create a global constraint (with possibly empty LHS)
  /** 
   * It's assumed that all processes that contribute to a global
   * constraint will call this routine with an equivalent Constraint
   * instance, i.e., the operator and RHS will be the same on all
   * proceses.  
   * 
   * @param name 
   * @param cons 
   */
  void createGlobalConstraint(const std::string& name, ConstraintPtr cons)
  {
    this->p_createGlobalConstraint(name, cons);
  }

  /// Add to the LHS of a global constraint, maybe
  void addToGlobalConstraint(const std::string& name, ExpressionPtr expr)
  {
    if (expr) this->p_addToGlobalConstraint(name, expr);
  }

  /// Add to the local part of the global objective function, maybe
  void addToObjective(ExpressionPtr expr)
  {
    if (expr) this->p_addToObjective(expr);
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

  /// Set file name for results (specialized)
  virtual void p_setFilename(std::string file) = 0;
  
  /// Add a (local) variable to be optimized (specialized)
  virtual void p_addVariable(VariablePtr v) = 0;

  /// Add a (local) auxiliary variable. This is used for
  /// parallel applications to create ghost copies of variables
  /// defined on other processors
  virtual void p_addAuxVariable(VariablePtr v) = 0;

  /// Add a (local) constraint (specialized)
  virtual void p_addConstraint(ConstraintPtr c) = 0;

  /// Create a global constraint (with possibly empty LHS)
  virtual void p_createGlobalConstraint(const std::string& name, ConstraintPtr cons) = 0;

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

  /// The auxiliary variables involved
  std::vector<VariablePtr> p_aux_variables;

  /// The collection of (local) constraints involved
  std::vector<ConstraintPtr> p_constraints;

  /// The (local part of the) objective fuction
  ExpressionPtr p_objective;

  /// Variables involved from all processes (unique instances)
  VarMap p_allVariables;

  /// Auxiliary variables from all processes
  VarMap p_exportVariables;

  /// Auxiliary variables from all processes
  VarMap p_auxVariables;

  /// Constraints involved from all processes 
  std::vector<ConstraintPtr> p_allConstraints;

  /// The entire objective function
  ExpressionPtr p_fullObjective;

  /// A thing to hold named constraints
  typedef std::map<std::string, ConstraintPtr> ConstraintMap;

  /// The global constraint expressions (local parts)
  ConstraintMap p_globalConstraints;

  /// The global constraints from all processes
  ConstraintMap p_allGlobalConstraints;

  /// Add a (local) variable to be optimized (specialized)
  void p_addVariable(VariablePtr v)
  {
    p_variables.push_back(v);
  }

  /// Add a (local) auxiliary variable. This is used for
  /// parallel applications to create ghost copies of variables
  /// defined on other processors
  virtual void p_addAuxVariable(VariablePtr v)
  {
    p_aux_variables.push_back(v);
  }

  /// Add a (local) constraint (specialized)
  void p_addConstraint(ConstraintPtr c)
  {
    p_constraints.push_back(c);
  }

  /// Create a global constraint (with possibly empty LHS) (specialized)
  void p_createGlobalConstraint(const std::string& name, ConstraintPtr cons);

  /// Add the local part of a global constraint (specialized)
  void p_addToGlobalConstraint(const std::string& name, ExpressionPtr expr);

  /// Get the global constraint expression
  ConstraintPtr p_getGlobalConstraint(const std::string& name);

  /// Add to the local part of the global objective function (added to other parts)
  void p_addToObjective(ExpressionPtr expr)
  {
    if (p_objective) {
      p_objective = p_objective + expr;
    } else {
      p_objective = expr;
    }
  }

  /// Gather the global constraints from a processor
  void p_gatherGlobalConstraints(const ConstraintMap& tmpglobal);

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

  void p_preconfigure(utility::Configuration::CursorPtr theprops);

  void p_setFilename(std::string file)
  {
    p_impl->setFilename(file);
  }

  /// Add a (local) variable to be optimized (specialized)
  void p_addVariable(VariablePtr v)
  {
    p_impl->addVariable(v);
  }

  /// Add a (local) auxiliary variable. This is used for
  /// parallel applications to create ghost copies of variables
  /// defined on other processors
  void p_addAuxVariable(VariablePtr v)
  {
    p_impl->addAuxVariable(v);
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

  /// Create a global constraint (with possibly empty LHS)
  void p_createGlobalConstraint(const std::string& name, ConstraintPtr cons)
  {
    p_impl->createGlobalConstraint(name, cons);
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
