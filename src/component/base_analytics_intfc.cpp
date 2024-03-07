/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

#include "gridpack/component/base_analytics_intfc.hpp"
#include "gridpack/parser/dictionary.hpp"

// Simple outline of data collection object

namespace gridpack{
namespace component{
/**
 * Simple constructor
 */
BaseBusAnalyticsInterface::BaseBusAnalyticsInterface(void)
{
  p_data = NULL;
}

/**
 * Simple destructor
 */
BaseBusAnalyticsInterface::~BaseBusAnalyticsInterface(void)
{
}

/**
 * Set data collection object inside interface
 * @param data_collection pointer to data collection object
 */
void BaseBusAnalyticsInterface::setData(
    boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  p_data = data.get();
}

/**
 * return number of generators on bus
 */
int BaseBusAnalyticsInterface::numGenerators()
{
  int ret = 0;
  if (p_data != NULL) {
    if (!p_data->getValue(GENERATOR_NUMBER,&ret)) ret = 0;
  }
  return ret;
}

/**
 * return number of loads on bus
 */
int BaseBusAnalyticsInterface::numLoads()
{
  int ret = 0;
  if (p_data != NULL) {
    if (!p_data->getValue(LOAD_NUMBER,&ret)) ret = 0;
  }
  return ret;
}

/**
 * Simple constructor
 */
BaseBranchAnalyticsInterface::BaseBranchAnalyticsInterface(void)
{
}

/**
 * Simple destructor
 */
BaseBranchAnalyticsInterface::~BaseBranchAnalyticsInterface(void)
{
}

/**
 * Set data collection object inside interface
 * data_collection pointer to data collection object
 */
void BaseBranchAnalyticsInterface::setData(
    boost::shared_ptr<gridpack::component::DataCollection> &data)
{
}

/**
 * return number of lines on branch
 */
int BaseBranchAnalyticsInterface::numLines()
{
  int ret = 0;
  if (p_data != NULL) {
    if (!p_data->getValue(BRANCH_NUM_ELEMENTS,&ret)) ret = 0;
  }
  return ret;
}

}  // component
}  // gridpack
