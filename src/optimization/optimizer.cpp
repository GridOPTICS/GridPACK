// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   optimizer.cpp
 * @author William A. Perkins
 * @date   2016-12-13 12:01:01 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <sstream>
#include <boost/bind.hpp>
#include <boost/format.hpp>

// These two includes are needed for Boost 1.56

#include <boost/serialization/singleton.hpp>
#include <boost/serialization/extended_type_info.hpp>

#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "optimizer.hpp"
#if defined(HAVE_CPLEX)
#include "cplex_optimizer_implementation.hpp"
#endif
#if defined(HAVE_GLPK)
#include "glpk_optimizer_implementation.hpp"
#endif
#include "lpfile_optimizer_implementation.hpp"
#include "julia_optimizer_implementation.hpp"

namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class VariableSubstituter
// -------------------------------------------------------------
class VariableSubstituter 
  : public ExpressionVisitor
{
public:

  /// Default constructor.
  VariableSubstituter(const OptimizerImplementation::VarMap& vmap)
    : ExpressionVisitor(), p_vmap(vmap)
  {}

  /// Destructor
  ~VariableSubstituter(void)
  {}

  /// Replace variables as needed
  void visit(VariableExpression& e)
  { 
    VariablePtr vold(e.var());
    OptimizerImplementation::VarMap::const_iterator v = p_vmap.find(vold->name());
    BOOST_ASSERT(v != p_vmap.end());
    VariablePtr vnew(v->second);

    BOOST_ASSERT(vold);
    BOOST_ASSERT(vnew);

    if (vold != vnew) {
      e.var(vnew);
    }
  }
protected:

  /// One variable per name
  const OptimizerImplementation::VarMap& p_vmap;

};


// -------------------------------------------------------------
//  class ConstraintRenamer
// -------------------------------------------------------------
class ConstraintRenamer 
  : public ExpressionVisitor
{
public:

  /// Default constructor.
  ConstraintRenamer(void)
    : ExpressionVisitor(), p_nextID(0)
  {}

  /// Destructor
  ~ConstraintRenamer(void)
  {}

protected:

  /// The constraint id
  int p_nextID;

  /// Rename constraints
  void visit(Constraint& e)
  {
    std::string nname = 
      boost::str(boost::format("C%d") % p_nextID++);
    e.name(nname);
  }
  
};


// -------------------------------------------------------------
//  class OptimizerImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// OptimizerImplementation::p_createGlobalConstraint
// -------------------------------------------------------------
void 
OptimizerImplementation::p_createGlobalConstraint(const std::string& name, ConstraintPtr cons)
{
  ConstraintPtr c;
  try {
    c = p_globalConstraints.at(name);
    std::string msg("p_createGlobalConstraint: Global constraint \"");
    msg += name;
    msg += "\" already exists";
    throw gridpack::Exception(msg);
  } catch (const std::out_of_range& e) {
    // do nothing
  }
  p_globalConstraints[name] = cons;
  p_globalConstraints[name]->name() = name;
}

// -------------------------------------------------------------
// OptimizerImplementation::p_addToGlobalConstraint
// -------------------------------------------------------------
void 
OptimizerImplementation::p_addToGlobalConstraint(const std::string& name, ExpressionPtr expr)
{
  ConstraintPtr c;
  try {
    c = p_globalConstraints.at(name);
  } catch (const std::out_of_range& e) {
    std::string msg("p_addToGlobalConstraint: Global constraint \"");
    msg += name;
    msg += "\" not defined";
    throw gridpack::Exception(msg);
  }
  c->addToLHS(expr);
}  

// -------------------------------------------------------------
// OptimizerImplementation::p_getGlobalConstraint
// -------------------------------------------------------------
ConstraintPtr 
OptimizerImplementation::p_getGlobalConstraint(const std::string& name)
{
  ConstraintPtr result;
  try {
    result = p_globalConstraints.at(name);
  } catch (const std::out_of_range& e) {
    std::string msg("p_getGlobalConstraint: Global constraint \"");
    msg += name;
    msg += "\" not defined";
    throw gridpack::Exception(msg);
  }
  return result;
}

// -------------------------------------------------------------
// OptimizerImplementation::p_gatherGlobalConstraints
// -------------------------------------------------------------
void
OptimizerImplementation::p_gatherGlobalConstraints(const ConstraintMap& tmpglobal)
{
  ConstraintMap::const_iterator c;
  for (c = tmpglobal.begin(); c != tmpglobal.end(); ++c) {
    std::string name(c->first);
    ConstraintPtr cons(c->second);
    if (cons->lhs()) {
      ConstraintMap::iterator gc(p_allGlobalConstraints.find(name));
      if (gc != p_allGlobalConstraints.end()) {
        gc->second->addToLHS(cons->lhs());
      } else {
        p_allGlobalConstraints[name] = cons;
      }
    }
  }
}


// -------------------------------------------------------------
// OptimizerImplementation::p_gatherProblem
// -------------------------------------------------------------
void
OptimizerImplementation::p_gatherProblem(void)
{
  parallel::Communicator comm(this->communicator());
  int nproc(comm.size());
  int me(comm.rank());

  // package up the local part of the problem into a string buffer;
  // MPI serialization cannot be used directly because VariablePtr's
  // are used in Expression's

  std::string lbuf;
  {
    std::ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa & p_variables;
    oa & p_constraints;
    oa & p_objective;
    oa & p_globalConstraints;
    lbuf = oss.str();
  }

  // all processes get a copy of the other processor parts of the problem 

  std::vector<std::string> gbuf(nproc);
  boost::mpi::all_gather(comm, lbuf, gbuf);

  // extract the problem and put it in local variables as required

  for (int p = 0; p < nproc; ++p) {
    if (p != me) {
      std::istringstream iss(gbuf[p]);
      boost::archive::binary_iarchive ia(iss);
      std::vector<VariablePtr> tmpvars;
      std::vector<ConstraintPtr> tmpcons;
      ConstraintMap tmpglobal;
      ExpressionPtr tmpobj;
      ia & tmpvars;
      ia & tmpcons;
      ia & tmpobj;
      ia & tmpglobal;

      for (std::vector<VariablePtr>::iterator v = tmpvars.begin();
           v != tmpvars.end(); ++v) {
        p_allVariables[(*v)->name()] = *v;
      }
      
      std::copy(tmpcons.begin(), tmpcons.end(), 
                std::back_inserter(p_allConstraints));

      // ConstraintMap::const_iterator c;
      // for (c = tmpglobal.begin(); c != tmpglobal.end(); ++c) {
      //   std::cout << me << ": " << p << ": " << c->second->name() << ": "
      //             << c->second->render() << std::endl;
      // }

      p_gatherGlobalConstraints(tmpglobal);

      if (tmpobj) {
        if (!p_fullObjective) {
          p_fullObjective = tmpobj;
        } else {
          p_fullObjective = p_fullObjective + tmpobj;
        }
      }
      
    } else {
      std::copy(p_constraints.begin(), p_constraints.end(), 
                std::back_inserter(p_allConstraints));

      // ConstraintMap::const_iterator c;
      // for (c = p_globalConstraints.begin(); c != p_globalConstraints.end(); ++c) {
      //   std::cout << me << ": " << p << ": " << c->second->name() << ": "
      //             << c->second->render() << std::endl;
      // }

      p_gatherGlobalConstraints(p_globalConstraints);

      if (p_objective) {
        if (!p_fullObjective) {
          p_fullObjective = p_objective;
        } else {
          p_fullObjective = p_fullObjective + p_objective;
        }
      }
    }
    comm.barrier();
  }


  // add global constraints to the constraints list
  ConstraintMap::const_iterator c;
  for (c = p_allGlobalConstraints.begin(); c != p_allGlobalConstraints.end(); ++c) {
    // std::cout << me << ": after: " << c->second->name() << ": "
    //           << c->second->render() << std::endl;
    p_allConstraints.push_back(c->second);
  }
  
  // make sure the locally defined variables are used on this processor

  for (std::vector<VariablePtr>::iterator v = p_variables.begin();
       v != p_variables.end(); ++v) {
    p_allVariables[(*v)->name()] = *v;
    p_exportVariables[(*v)->name()] = *v;
  }
  for (std::vector<VariablePtr>::iterator v = p_aux_variables.begin();
       v != p_aux_variables.end(); ++v) {
    p_allVariables[(*v)->name()] = *v;
  }

  // subsitute variables in constraints and objective so they are unique

  VariableSubstituter vs(p_allVariables);
  std::for_each(p_allConstraints.begin(), p_allConstraints.end(),
                boost::bind(&Constraint::accept, _1, boost::ref(vs)));
  if (p_fullObjective) {
    p_fullObjective->accept(vs);
  }

  // uniquely name all constraints in parallel
  if (nproc > 1) {
    ConstraintRenamer r;
    std::for_each(p_allConstraints.begin(), p_allConstraints.end(),
                  boost::bind(&Constraint::accept, _1, boost::ref(r)));
  }
}


// -------------------------------------------------------------
//  class Optimizer
// -------------------------------------------------------------

// -------------------------------------------------------------
// Optimizer:: constructors / destructor
// -------------------------------------------------------------
Optimizer::Optimizer(const parallel::Communicator& comm)
  : OptimizerInterface(),
    parallel::WrappedDistributed(),
    utility::WrappedConfigurable(),
    p_impl()
{
  p_setImpl(new LPFileOptimizerImplementation(comm));
}

Optimizer::~Optimizer(void)
{
}

// -------------------------------------------------------------
// Optimizer::p_preconfigure
// -------------------------------------------------------------
void
Optimizer::p_preconfigure(utility::Configuration::CursorPtr theprops)
{
  parallel::Communicator comm(p_impl->communicator());
  std::string key(p_impl->configurationKey());
  utility::Configuration::CursorPtr p = theprops->getCursor(key);
  std::string solver;
  solver = p->get("Solver", solver);

  if (!solver.empty()) {
    p_setImpl(NULL);
    if (solver == "GLPK") {
#if defined(HAVE_GLPK)
      p_setImpl(new GLPKOptimizerImplementation(comm));
#else
      throw gridpack::Exception("GLPK Optimizer not supported"); 
#endif
    } 
    if (solver == "CPLEX") {
#if defined(HAVE_CPLEX)
      p_setImpl(new CPlexOptimizerImplementation(comm));
#else
      throw gridpack::Exception("CPLEX Optimizer not supported"); 
#endif
    } 
    if (solver == "Julia") {
      p_setImpl(new JuliaOptimizerImplementation(comm));
    }
    if (solver == "LPFile") {
      p_setImpl(new LPFileOptimizerImplementation(comm));
    }
    if (!p_impl) {
      std::string s("Unknown ptimizer solver type \"");
      s += solver;
      s += "\"";
      throw gridpack::Exception(s); 
    }
  }
}

} // namespace optimization
} // namespace gridpack




