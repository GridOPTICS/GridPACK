/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   optimization_ifc.hpp
 * @author Yilin Fang and Bruce Palmer
 * @date   2015-06-12
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/smart_ptr/shared_ptr.hpp>
#include "gridpack/component/optimization_ifc.hpp"

namespace gridpack {
namespace component {
// base implementation for the optimization interface
/**
 * Constructor
 */
OptimizationInterface::OptimizationInterface(void)
{
}

/**
 * Constructor
 */
OptimizationInterface::~OptimizationInterface(void)
{
}

/**
 * Return a vector of optimization variables associated witht this
 * interface
 * @return list of variables
 */
std::vector<boost::shared_ptr<gridpack::optimization::Variable> >
  OptimizationInterface::getVariables()
{
  std::vector<boost::shared_ptr<gridpack::optimization::Variable> > ret;
  return ret;
}

/**
 * Return a vector of auxiliary variables associated with this interface.
 * These are variables that are used in expressions but may not be
 * defined by this network
 * @return list of variables
 */
std::vector<boost::shared_ptr<gridpack::optimization::Variable> >
          OptimizationInterface::getAuxVariables()
{
  std::vector<boost::shared_ptr<gridpack::optimization::Variable> > ret;
  return ret;
}

/**
 * Return contribution from bus to a global constraint
 * @param tag string that can be parsed by bus to determine which constraint
 * contribution is being requested
 * @return contribution to global constraint. If no contribution, return null
 * pointer
 */
boost::shared_ptr<gridpack::optimization::Expression>
  OptimizationInterface::getGlobalConstraint(const char* tag)
{
  boost::shared_ptr<gridpack::optimization::Expression> ret;
  return ret;
}

/**
 * Return a list of local constraints from component
 * @return list of constraints
 */
std::vector<boost::shared_ptr<gridpack::optimization::Constraint> >
  OptimizationInterface::getLocalConstraints()
{
  std::vector<boost::shared_ptr<gridpack::optimization::Constraint> > ret;
  return ret;
}

/**
 * Return contribution to objective function
 * @return expression representing contribution to objective function. If no
 * contribution, return null pointer
 */
boost::shared_ptr<gridpack::optimization::Expression>
  OptimizationInterface::getObjectiveFunction()
{
  boost::shared_ptr<gridpack::optimization::Expression> ret;
  return ret;
}

}  // component
}  // gridpack
