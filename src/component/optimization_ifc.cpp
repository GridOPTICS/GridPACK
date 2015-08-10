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
std::vector<gridpack::optimization::Variable*> OptimizationInterface::getVariables()
{
  std::vector<gridpack::optimization::Variable*> ret;
  return ret;
}

/**
 * Return contribution from bus to a global constraint
 * @param tag string that can be parsed by bus to determine which constraint
 * contribution is being requested
 * @param flag bool that returns false if there is no contribution to constaint
 * @return contribution to global constraint
 */
gridpack::optimization::Expression*
  OptimizationInterface::getGlobalConstraint(const char* tag, bool *flag)
{
  gridpack::optimization::Expression *ret = NULL;
  *flag = false;
  return ret;
}

}  // component
}  // gridpack
