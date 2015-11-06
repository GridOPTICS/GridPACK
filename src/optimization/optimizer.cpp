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
 * @date   2015-11-04 09:00:53 d3g096
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
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "optimizer.hpp"
#if defined(HAVE_CPLEX)
#include "cplex_optimizer_implementation.hpp"
#elif defined(HAVE_GLPK)
#include "glpk_optimizer_implementation.hpp"
#else
#include "lpfile_optimizer_implementation.hpp"
#endif

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
      ExpressionPtr tmpobj;
      ia & tmpvars;
      ia & tmpcons;
      ia & tmpobj;

      for (std::vector<VariablePtr>::iterator v = tmpvars.begin();
           v != tmpvars.end(); ++v) {
        p_allVariables[(*v)->name()] = *v;
      }
      
      std::copy(tmpcons.begin(), tmpcons.end(), 
                std::back_inserter(p_allConstraints));

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
      if (p_objective) {
        if (!p_fullObjective) {
          p_fullObjective = p_objective;
        } else {
          p_fullObjective = p_fullObjective + p_objective;
        }
      }
    }
  }

  comm.barrier();
  // make sure the locally defined variables are used on this processor

  for (std::vector<VariablePtr>::iterator v = p_variables.begin();
       v != p_variables.end(); ++v) {
    p_allVariables[(*v)->name()] = *v;
  }

  // subsitute variables in constraints and objective so they are unique

  VariableSubstituter vs(p_allVariables);
  std::for_each(p_allConstraints.begin(), p_allConstraints.end(),
                boost::bind(&Constraint::accept, _1, boost::ref(vs)));
  p_fullObjective->accept(vs);

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
  p_setImpl(
#if defined(HAVE_CPLEX)
            new CPlexOptimizerImplementation(comm)
#elif defined(HAVE_GLPK)
            new GLPKOptimizerImplementation(comm)
#else
            new LPFileOptimizerImplementation(comm)
#endif
            );
}

Optimizer::~Optimizer(void)
{
}

} // namespace optimization
} // namespace gridpack




