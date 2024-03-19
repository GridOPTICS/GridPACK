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
    boost::shared_ptr<gridpack::component::DataCollection> data)
{
  p_data = data.get();
}

/**
 * Return pointer to data collection object stored in the base interface
 */
gridpack::component::DataCollection* BaseBusAnalyticsInterface::getData()
{
  return p_data;
}


/**
 * return number of generators on bus
 */
int BaseBusAnalyticsInterface::numGenerators()
{
  int ret = 0;
  if (p_data != NULL) {
    int i, ngen, gcnt, istat;
    if (!p_data->getValue(GENERATOR_NUMBER,&ngen)) ngen = 0;
    gcnt = 0;
    for (i=0; i<ngen; i++) {
      if (p_data->getValue(GENERATOR_STAT,&istat,i)) {
        if (istat > 0) gcnt++;
      } else {
        gcnt++;
      }
    }
    ret = gcnt;
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
    int i, nload, lcnt, istat;
    if (!p_data->getValue(LOAD_NUMBER,&nload)) nload = 0;
    lcnt = 0;
    for (i=0; i<nload; i++) {
      if (p_data->getValue(LOAD_STATUS,&istat,i)) {
        if (istat > 0) lcnt++;
      } else {
        lcnt++;
      }
    }
    ret = lcnt;
  }
  return ret;
}

/**
 * return number of storage units on bus
 */
int BaseBusAnalyticsInterface::numStorage()
{
  return 0;
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
    boost::shared_ptr<gridpack::component::DataCollection> data)
{
  p_data = data.get();
}

/**
 * Return pointer to data collection object stored in the base interface
 */
gridpack::component::DataCollection* BaseBranchAnalyticsInterface::getData()
{
  return p_data;
}

/**
 * return number of lines on branch
 */
int BaseBranchAnalyticsInterface::numLines()
{
  int ret = 0;
  if (p_data != NULL) {
    int i, nline, lcnt, istat;
    if (!p_data->getValue(BRANCH_NUM_ELEMENTS,&nline)) nline = 0;
    for (i=0; i<nline; i++) {
      if (p_data->getValue(BRANCH_STATUS,&istat,i)) {
        if (istat > 0) lcnt++;
      } else {
        lcnt++;
      }
    }
    ret = lcnt;
  }
  return ret;
}

}  // component
}  // gridpack
