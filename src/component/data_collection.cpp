#include "gridpack/component/data_collection.hpp"

/**
 * Simple constructor
 */
gridpack::component::DataCollection::DataCollection(void)
{
}

/**
 * Simple destructor
 */
gridpack::component::DataCollection::~DataCollection(void)
{
}

/**
 * Assignment operator
 */
gridpack::component::DataCollection & gridpack::component::DataCollection::operator=
  (const gridpack::component::DataCollection &rhs)
{
  if (this == &rhs) return *this;
  p_ints = rhs.p_ints;
  p_longs = rhs.p_longs;
  p_bools = rhs.p_bools;
  p_strings = rhs.p_strings;
  p_floats = rhs.p_floats;
  p_doubles = rhs.p_doubles;
  p_complexType = rhs.p_complexType;
}

/**
 *  Add variables to DataCollection object
 *  @param name: name given to data element
 *  @param value: value of data element
 */
void gridpack::component::DataCollection::addValue(char *name, int value)
{
  std::string str = name;
  p_ints.insert(std::pair<std::string, int>(str,value));
}

void gridpack::component::DataCollection::addValue(char *name, long value)
{
  std::string str = name;
  p_longs.insert(std::pair<std::string, long>(str,value));
}

void gridpack::component::DataCollection::addValue(char *name, bool value)
{
  std::string str = name;
  p_bools.insert(std::pair<std::string, bool>(str,value));
}

void gridpack::component::DataCollection::addValue(char *name, char *value)
{
  std::string str = name;
  std::string val = value;
  p_strings.insert(std::pair<std::string, std::string>(str,value));
}

void gridpack::component::DataCollection::addValue(char *name, float value)
{
  std::string str = name;
  p_floats.insert(std::pair<std::string, float>(str,value));
}

void gridpack::component::DataCollection::addValue(char *name, double value)
{
  std::string str = name;
  p_doubles.insert(std::pair<std::string, double>(str,value));
}

#if 1
void gridpack::component::DataCollection::addValue(char *name, gridpack::ComplexType value)
{
  std::string str = name;
  p_complexType.insert(std::pair<std::string, gridpack::ComplexType>(str,value));
}
#endif

/**
 *  Modify current value of existing data element in
 *  DataCollection object
 *  @param name: name of data element
 *  @param value: new value of data element
 *  @return: false if no element of the correct name and type exists in
 *  DataCollection object
 */
bool gridpack::component::DataCollection::setValue(char *name, int value)
{
  std::string str = name;
  std::map<std::string, int>::iterator it;
  it = p_ints.find(str);
  if (it != p_ints.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::setValue(char *name, long value)
{
  std::string str = name;
  std::map<std::string, long>::iterator it;
  it = p_longs.find(str);
  if (it != p_longs.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::setValue(char *name, bool value)
{
  std::string str = name;
  std::map<std::string, bool>::iterator it;
  it = p_bools.find(str);
  if (it != p_bools.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::setValue(char *name, char *value)
{
  std::string str = name;
  std::string val = value;
  std::map<std::string, std::string>::iterator it;
  it = p_strings.find(str);
  if (it != p_strings.end()) {
    it->second = val;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::setValue(char *name, float value)
{
  std::string str = name;
  std::map<std::string, float>::iterator it;
  it = p_floats.find(str);
  if (it != p_floats.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::setValue(char *name, double value)
{
  std::string str = name;
  std::map<std::string, double>::iterator it;
  it = p_doubles.find(str);
  if (it != p_doubles.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

#if 1
bool gridpack::component::DataCollection::setValue(char *name, gridpack::ComplexType value)
{
  std::string str = name;
  std::map<std::string, gridpack::ComplexType>::iterator it;
  it = p_complexType.find(str);
  if (it != p_complexType.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}
#endif

/**
 *  Retrieve current value of existing data element in
 *  DataCollection object
 *  @param name: name of data element
 *  @param value: current value of data element
 *  @return: false if no element of the correct name and type exists in
 *  DataCollection object
 */
bool gridpack::component::DataCollection::getValue(char *name, int *value)
{
  std::string str = name;
  std::map<std::string, int>::iterator it;
  it = p_ints.find(str);
  if (it != p_ints.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::getValue(char *name, long *value)
{
  std::string str = name;
  std::map<std::string, long>::iterator it;
  it = p_longs.find(str);
  if (it != p_longs.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::getValue(char *name, bool *value)
{
  std::string str = name;
  std::map<std::string, bool>::iterator it;
  it = p_bools.find(str);
  if (it != p_bools.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::getValue(char *name, std::string *value)
{
  std::string str = name;
  std::map<std::string, std::string>::iterator it;
  it = p_strings.find(str);
  if (it != p_strings.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::getValue(char *name, float *value)
{
  std::string str = name;
  std::map<std::string, float>::iterator it;
  it = p_floats.find(str);
  if (it != p_floats.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::getValue(char *name, double *value)
{
  std::string str = name;
  std::map<std::string, double>::iterator it;
  it = p_doubles.find(str);
  if (it != p_doubles.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

#if 1
bool gridpack::component::DataCollection::getValue(char *name, gridpack::ComplexType *value)
{
  std::string str = name;
  std::map<std::string, gridpack::ComplexType>::iterator it;
  it = p_complexType.find(str);
  if (it != p_complexType.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}
#endif
