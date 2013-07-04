#ifndef _data_collection_h
#define _data_collection_h

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
   *  @param name: name given to data element
   *  @param value: value of data element
   */
  void addValue(char *name, int value);
  void addValue(char *name, long value);
  void addValue(char *name, bool value);
  void addValue(char *name, char *value);
  void addValue(char *name, float value);
  void addValue(char *name, double value);
  void addValue(char *name, gridpack::ComplexType value);

  /**
   *  Modify current value of existing data element in
   *  DataCollection object
   *  @param name: name of data element
   *  @param value: new value of data element
   *  @return: false if no element of the correct name and type exists in
   *  DataCollection object
   */
  bool setValue(char *name, int value);
  bool setValue(char *name, long value);
  bool setValue(char *name, bool value);
  bool setValue(char *name, char *value);
  bool setValue(char *name, float value);
  bool setValue(char *name, double value);
  bool setValue(char *name, gridpack::ComplexType value);

  /**
   *  Retrieve current value of existing data element in
   *  DataCollection object
   *  @param name: name of data element
   *  @param value: current value of data element
   *  @return: false if no element of the correct name and type exists in
   *  DataCollection object
   */
  bool getValue(char *name, int *value);
  bool getValue(char *name, long *value);
  bool getValue(char *name, bool *value);
  bool getValue(char *name, std::string *value);
  bool getValue(char *name, float *value);
  bool getValue(char *name, double *value);
  bool getValue(char *name, gridpack::ComplexType *value);
private:
  std::map<std::string, int> p_ints; 
  std::map<std::string, long> p_longs; 
  std::map<std::string, bool> p_bools; 
  std::map<std::string, std::string> p_strings; 
  std::map<std::string, float> p_floats; 
  std::map<std::string, double> p_doubles; 
  std::map<std::string, gridpack::ComplexType> p_complexType; 

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
