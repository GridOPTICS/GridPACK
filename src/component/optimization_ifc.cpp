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
 * Get the number of optimization variables contributed from component
 * @return number of optimization variables
 */
int OptimizationInterface::numOptVariables(void)
{
  return 0;
}

/**
 * Return the data type of the optimization variable
 * @param idx index of variable inside network component
 * @return data type of network variable
 */
int OptimizationInterface::optVariableType(int idx)
{
  return OPT_DBL;
}

/**
 * Return bounds on variable. If the lower bound is a null pointer the lower
 * bound is -inf, if the upper bound is a null pointer, the upper bound is
 * inf
 * @param idx index of variable inside network component
 * @param lo lower bound
 * @param hi upper bound
 */
void OptimizationInterface::optVariableBounds(int idx, double *lo, double *hi)
{
  lo = NULL;
  hi = NULL;
}

void OptimizationInterface::optVariableBounds(int idx, int *lo, int *hi)
{
  lo = NULL;
  hi = NULL;
}


/**
 * Return the objective function contributed from the bus 
 */
double OptimizationInterface::objectiveFunction(void)
{
  double value;
  value = 0.0;
  return value;
}

/**
 * Return the solution from the bus 
 */
bool OptimizationInterface::solution(void)
{
  return false;
}

}  // component
}  // gridpack
