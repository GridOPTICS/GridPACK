/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#ifndef _data_collection_h
#define _data_collection_h

#define OLD_MAP

#ifndef OLD_MAP
#include <boost/unordered_map.hpp>
#endif
#include <map>
#include <string>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>

#include "gridpack/utilities/complex.hpp"

#include <boost/serialization/export.hpp>

// Simple outline of data collection object

namespace gridpack{
namespace component{

class DataCollection {
public:
  /**
   * Simple constructor
   */
  DataCollection(void);

  /**
   * Simple destructor
   */
  ~DataCollection(void);

  /**
   * Assignment operator
   */
  DataCollection & operator= (const DataCollection &rhs);

  /**
   *  Add variables to DataCollection object
   *  @param name name given to data element
   *  @param value value of data element
   */
  void addValue(const char *name, const int value);
  void addValue(const char *name, const long value);
  void addValue(const char *name, const bool value);
  void addValue(const char *name, const char *value);
  void addValue(const char *name, const float value);
  void addValue(const char *name, const double value);
  void addValue(const char *name, const gridpack::ComplexType value);

  /**
   *  Add variables to DataCollection object with an additional index to keep
   *  track of items that can appear more than once. Item appears in
   *  DataCollection with the tag "name:idx"
   *  @param name name given to data element
   *  @param value value of data element
   *  @param idx index of value
   */
  void addValue(const char *name, const int value, const int idx);
  void addValue(const char *name, const long value, const int idx);
  void addValue(const char *name, const bool value, const int idx);
  void addValue(const char *name, const char *value, const int idx);
  void addValue(const char *name, const float value, const int idx);
  void addValue(const char *name, const double value, const int idx);
  void addValue(const char *name, const gridpack::ComplexType value, const int idx);

  /**
   *  Modify current value of existing data element in
   *  DataCollection object
   *  @param name name of data element
   *  @param value new value of data element
   *  @return false if no element of the correct name and type exists in
   *  DataCollection object
   */
  bool setValue(const char *name, const int value);
  bool setValue(const char *name, const long value);
  bool setValue(const char *name, const bool value);
  bool setValue(const char *name, const char *value);
  bool setValue(const char *name, const float value);
  bool setValue(const char *name, const double value);
  bool setValue(const char *name, const gridpack::ComplexType value);

  /**
   *  Modify current value of existing data element in
   *  DataCollection object. Assume that name appears in DataCollection with an
   *  additional index in the form "name:idx"
   *  @param name name of data element
   *  @param value new value of data element
   *  @param idx index of value
   *  @return false if no element of the correct name and type exists in
   *  DataCollection object
   */
  bool setValue(const char *name, const int value, const int idx);
  bool setValue(const char *name, const long value, const int idx);
  bool setValue(const char *name, const bool value, const int idx);
  bool setValue(const char *name, const char *value, const int idx);
  bool setValue(const char *name, const float value, const int idx);
  bool setValue(const char *name, const double value, const int idx);
  bool setValue(const char *name, const gridpack::ComplexType value, const int idx);

  /**
   *  Retrieve current value of existing data element in
   *  DataCollection object
   *  @param name name of data element
   *  @param value current value of data element
   *  @return false if no element of the correct name and type exists in
   *  DataCollection object
   */
  bool getValue(const char *name, int *value);
  bool getValue(const char *name, long *value);
  bool getValue(const char *name, bool *value);
  bool getValue(const char *name, std::string *value);
  bool getValue(const char *name, float *value);
  bool getValue(const char *name, double *value);
  bool getValue(const char *name, gridpack::ComplexType *value);

  /**
   *  Retrieve current value of existing data element in
   *  DataCollection object. Assume that item appears in DataCollection with the
   *  tag "name:idx"
   *  @param name name of data element
   *  @param value current value of data element
   *  @param idx index of value
   *  @return false if no element of the correct name and type exists in
   *  DataCollection object
   */
  bool getValue(const char *name, int *value, const int idx);
  bool getValue(const char *name, long *value, const int idx);
  bool getValue(const char *name, bool *value, const int idx);
  bool getValue(const char *name, std::string *value, const int idx);
  bool getValue(const char *name, float *value, const int idx);
  bool getValue(const char *name, double *value, const int idx);
  bool getValue(const char *name, gridpack::ComplexType *value, const int idx);

  /**
   * Dump contents of data collection to standard out
   */
  void dump(void);
private:
#ifdef OLD_MAP
  std::map<std::string, int> p_ints; 
  std::map<std::string, long> p_longs; 
  std::map<std::string, bool> p_bools; 
  std::map<std::string, std::string> p_strings; 
  std::map<std::string, float> p_floats; 
  std::map<std::string, double> p_doubles; 
  std::map<std::string, gridpack::ComplexType> p_complexType; 
#else
  boost::unordered_map<std::string, int> p_ints; 
  boost::unordered_map<std::string, long> p_longs; 
  boost::unordered_map<std::string, bool> p_bools; 
  boost::unordered_map<std::string, std::string> p_strings; 
  boost::unordered_map<std::string, float> p_floats; 
  boost::unordered_map<std::string, double> p_doubles; 
  boost::unordered_map<std::string, gridpack::ComplexType> p_complexType; 
#endif

private:
  friend class boost::serialization::access;

  /// Serialization method
  template<class Archive> void serialize(Archive &ar, const unsigned int)
  {
    ar & p_ints
      & p_longs
      & p_bools
      & p_strings
      & p_floats
      & p_doubles
      & p_complexType;
  }

};


}    // component
}    // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::component::DataCollection);

#endif // _data_collection_h
