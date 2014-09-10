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
  FortranData *data;
};

/**
 * Retrieve current value of existing data element in
 * data_collection object
 * @param data pointer to GridPACK data collection wrapper object
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
extern "C" bool data_collection_get_logical_value(dataWrapper *wdata, char *name,
      bool *value)
{
  return wdata->data->getValue(name,value);
}
extern "C" bool data_collection_get_string_value(dataWrapper *wdata, char *name,
      char *value)
{
  return wdata->data->getValue(name,value);
}
extern "C" bool data_collection_get_real_value(dataWrapper *wdata, char *name,
      float *value)
{
  return wdata->data->getValue(name,value);
}
extern "C" bool data_collection_get_double_value(dataWrapper *wdata, char *name,
      double *value)
{
  return wdata->data->getValue(name,value);
}
/**
 * Retrieve current value of existing data element in
 * data_collection object. Assume that the item appears in data_collection with
 * tag "name:idx"
 * @param p_data pointer to GridPACK data collection wrapper object
 * @param name name of data element
 * @param value current value of data element
 * @param idx index of value
 * @return false if no element of the correct name and type exists in
 * data_collection object
 */
extern "C" bool data_collection_get_int_indexed_value(dataWrapper *wdata, char *name,
      int *value, int idx)
{
  return wdata->data->getValue(name,value,idx);
}
extern "C" bool data_collection_get_logical_indexed_value(dataWrapper *wdata, char *name,
      bool *value, int idx)
{
  return wdata->data->getValue(name,value,idx);
}
extern "C" bool data_collection_get_string_indexed_value(dataWrapper *wdata, char *name,
      char *value, int idx)
{
  return wdata->data->getValue(name,value,idx);
}
extern "C" bool data_collection_get_real_indexed_value(dataWrapper *wdata, char *name,
      float *value, int idx)
{
  return wdata->data->getValue(name,value,idx);
}
extern "C" bool data_collection_get_double_indexed_value(dataWrapper *wdata, char *name,
      double *value, int idx)
{
  return wdata->data->getValue(name,value,idx);
}
