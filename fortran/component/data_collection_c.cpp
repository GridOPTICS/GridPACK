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
#include <string.h>
#include <boost/smart_ptr/shared_ptr.hpp>
#include "gridpack/component/data_collection.hpp"
#include <stdio.h>

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
extern "C" bool data_collection_get_int_value(FortranData *data, char *name,
      int *value)
{
  return data->getValue(name,value);
}
extern "C" bool data_collection_get_logical_value(FortranData *data, char *name,
      bool *value)
{
  return data->getValue(name,value);
}
extern "C" bool data_collection_get_string_value(FortranData *data, char *name,
      char *value)
{
  std::string s_value;
  bool ret = data->getValue(name,&s_value);
  strcpy(value,s_value.c_str());
  return ret;
}
extern "C" bool data_collection_get_real_value(FortranData *data, char *name,
      float *value)
{
  return data->getValue(name,value);
}
extern "C" bool data_collection_get_double_value(FortranData *data, char *name,
      double *value)
{
  return data->getValue(name,value);
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
extern "C" bool data_collection_get_int_indexed_value(FortranData *data, char *name,
      int *value, int idx)
{
  idx--;
  return data->getValue(name,value,idx);
}
extern "C" bool data_collection_get_logical_indexed_value(FortranData *data, char *name,
      bool *value, int idx)
{
  idx--;
  return data->getValue(name,value,idx);
}
extern "C" bool data_collection_get_string_indexed_value(FortranData *data, char *name,
      char *value, int idx)
{
  idx--;
  std::string s_value;
  bool ret = data->getValue(name,&s_value,idx);
  strcpy(value,s_value.c_str());
  return ret;
}
extern "C" bool data_collection_get_real_indexed_value(FortranData *data, char *name,
      float *value, int idx)
{
  idx--;
  return data->getValue(name,value,idx);
}
extern "C" bool data_collection_get_double_indexed_value(FortranData *data, char *name,
      double *value, int idx)
{
  idx--;
  return data->getValue(name,value,idx);
}
