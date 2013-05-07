#ifndef _data_collection_h
#define _data_collection_h

#include <map>
/**
 * Simple constructor
 */
DataCollection::DataCollection(void)
{
}

/**
 * Simple destructor
 */
DataCollection::~DataCollection(void)
{
}

/**
 *  Add variables to DataCollection object
 *  @param name: name given to data element
 *  @param value: value of data element
 */
void DataCollection::addValue(char *name, int value)
{
  string str = name;
  p_ints.insert(std::pair<std::string, int>(str,value));
}

void DataCollection::addValue(char *name, long value)
{
  string str = name;
  p_longs.insert(std::pair<std::string, long>(str,value));
}

void DataCollection::addValue(char *name, bool value)
{
  string str = name;
  p_bools.insert(std::pair<std::string, bool>(str,value));
}

void DataCollection::addValue(char *name, char *value)
{
  string str = name;
  string val = value;
  p_strings.insert(std::pair<std::string, std::string>(str,value));
}

void DataCollection::addValue(char *name, float value)
{
  string str = name;
  p_floats.insert(std::pair<std::string, float>(str,value));
}

void DataCollection::addValue(char *name, double value)
{
  string str = name;
  p_doubles.insert(std::pair<std::string, double>(str,value));
}

void DataCollection::addValue(char *name, SingleComplex value)
{
  string str = name;
  p_singleComplex.insert(std::pair<std::string, singleComplex>(str,value));
}

void DataCollection::addValue(char *name, DoubleComplex value)
{
  string str = name;
  p_doubleComplex.insert(std::pair<std::string, doubleComplex>(str,value));
}

/**
 *  Modify current value of existing data element in
 *  DataCollection object
 *  @param name: name of data element
 *  @param value: new value of data element
 *  @return: false if no element of the correct name and type exists in
 *  DataCollection object
 */
bool DataCollection::setValue(char *name, int value)
{
  string str = name;
  std::map<std::string, int>::iterator it;
  it = p_ints.find(str);
  if (it != p_ints.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::setValue(char *name, long value)
{
  string str = name;
  std::map<std::string, long>::iterator it;
  it = p_longs.find(str);
  if (it != p_longs.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::setValue(char *name, bool value)
{
  string str = name;
  std::map<std::string, bool>::iterator it;
  it = p_bools.find(str);
  if (it != p_bools.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::setValue(char *name, char *value)
{
  string str = name;
  string val = value;
  std::map<std::string, std::string>::iterator it;
  it = p_strings.find(str);
  if (it != p_strings.end()) {
    it->second = val;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::setValue(char *name, float value)
{
  string str = name;
  std::map<std::string, float>::iterator it;
  it = p_floats.find(str);
  if (it != p_floats.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::setValue(char *name, double value)
{
  string str = name;
  std::map<std::string, double>::iterator it;
  it = p_doubles.find(str);
  if (it != p_doubles.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::setValue(char *name, SingleComplex value)
{
  string str = name;
  std::map<std::string, singleComplex>::iterator it;
  it = p_singleComplex.find(str);
  if (it != p_singleComplex.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::setValue(char *name, DoubleComplex value)
{
  string str = name;
  std::map<std::string, doubleComplex>::iterator it;
  it = p_doubleComplex.find(str);
  if (it != p_doubleComplex.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

/**
 *  Retrieve current value of existing data element in
 *  DataCollection object
 *  @param name: name of data element
 *  @param value: current value of data element
 *  @return: false if no element of the correct name and type exists in
 *  DataCollection object
 */
bool DataCollection::getValue(char *name, int *value)
{
  string str = name;
  std::map<std::string, int>::iterator it;
  it = p_ints.find(str);
  if (it != p_ints.end()) {
    value = it->second;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::getValue(char *name, long *value)
{
  string str = name;
  std::map<std::string, long>::iterator it;
  it = p_longs.find(str);
  if (it != p_longs.end()) {
    value = it->second;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::getValue(char *name, bool *value)
{
  string str = name;
  std::map<std::string, bool>::iterator it;
  it = p_bools.find(str);
  if (it != p_bools.end()) {
    value = it->second;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::getValue(char *name, std::string *value)
{
  string str = name;
  std::map<std::string, std::string>::iterator it;
  it = p_strings.find(str);
  if (it != p_bools.end()) {
    value = it->second;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::getValue(char *name, float *value)
{
  string str = name;
  std::map<std::string, float>::iterator it;
  it = p_floats.find(str);
  if (it != p_floats.end()) {
    value = it->second;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::getValue(char *name, double *value)
{
  string str = name;
  std::map<std::string, double>::iterator it;
  it = p_doubles.find(str);
  if (it != p_doubles.end()) {
    value = it->second;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::getValue(char *name, SingleComplex *value)
{
  string str = name;
  std::map<std::string, singleComplex>::iterator it;
  it = p_singleComplex.find(str);
  if (it != p_singleComplex.end()) {
    value = it->second;
    return true;
  } else {
    return false;
  }
}

bool DataCollection::getValue(char *name, DoubleComplex *value)
{
  string str = name;
  std::map<std::string, doubleComplex>::iterator it;
  it = p_doubleComplex.find(str);
  if (it != p_doubleComplex.end()) {
    value = it->second;
    return true;
  } else {
    return false;
  }
}
