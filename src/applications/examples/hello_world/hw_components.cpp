/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_components.cpp
 * @author Bruce Palmer
 * @date   2013-11-19 13:32:05 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/applications/examples/hello_world/hw_components.hpp"
#include "gridpack/parser/dictionary.hpp"

//#define LARGE_MATRIX

/**
 *  Simple constructor
 */
gridpack::hello_world::HWBus::HWBus(void)
{
  p_original_idx = 0;
}

/**
 *  Simple destructor
 */
gridpack::hello_world::HWBus::~HWBus(void)
{
}


/**
 * Load values stored in DataCollection object into HWBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::hello_world::HWBus::load(const
         boost::shared_ptr<gridpack::component::DataCollection> &data)
{
   data->getValue(BUS_NUMBER,&p_original_idx);
}

/**
 * Write output from buses to standard out
 * @param string (output) string with information to be printed out
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool gridpack::hello_world::HWBus::serialWrite(char *string, const char *signal)
{
  sprintf(string,"Hello world from bus %d\n",p_original_idx);
  return true;
}

/**
 *  Simple constructor
 */
gridpack::hello_world::HWBranch::HWBranch(void)
{
  p_original_idx1 = 0;
  p_original_idx2 = 0;
}

/**
 *  Simple destructor
 */
gridpack::hello_world::HWBranch::~HWBranch(void)
{
}

/**
 * Load values stored in DataCollection object into HWBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::hello_world::HWBranch::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  data->getValue(BRANCH_FROMBUS,&p_original_idx1);
  data->getValue(BRANCH_TOBUS,&p_original_idx2);
}

/**
 * Write output from branches to standard out
 * @param string (output) string with information to be printed out
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if branch is contributing string to output, false otherwise
 */
bool gridpack::hello_world::HWBranch::serialWrite(char *string, const char *signal)
{
  sprintf(string,"Hello world from the branch connecting bus %d to bus %d\n",
      p_original_idx1, p_original_idx2);
  return true;
}
