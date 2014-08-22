/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   data_collection_c.cpp
 * @author Bruce Palmer
 * @date   2014-08-15 11:05:08 d3g293
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
#include <boost/smart_ptr/shared_ptr.hpp>
#include "gridpack/component/data_collection.hpp"

typedef gridpack::component::DataCollection FortranData;

struct dataWrapper {
  boost::shared_ptr<FortranData> data;
};

/**
 * Retrieve current value of existing data element in
 * data_collection object
 * @param data pointer to GridPACK data collection object
 * @param name name of data element
 * @param value current value of data element
 * @return false if no element of the correct name and type exists in
 * data_collection object
 */
extern "C" bool data_collection_get_int_value(dataWrapper *wdata, char *name,
      int *value)
{
  return wdata->data->getValue(name,value);
}
