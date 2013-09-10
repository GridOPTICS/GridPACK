#include "gridpack/component/data_collection.hpp"
#include <iostream.h>

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
 *  @param name name given to data element
 *  @param value value of data element
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

void gridpack::component::DataCollection::addValue(char *name, gridpack::ComplexType value)
{
  std::string str = name;
  p_complexType.insert(std::pair<std::string, gridpack::ComplexType>(str,value));
}

/**
 *  Add variables to DataCollection object with an additional index to keep
 *  track of items that can appear more than once. Item appears in
 *  DataCollection with the tag "name:idx"
 *  @param name name given to data element
 *  @param value value of data element
 *  @param idx index of value
 */
void gridpack::component::DataCollection::addValue(char *name, int value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  p_ints.insert(std::pair<std::string, int>(str,value));
}

void gridpack::component::DataCollection::addValue(char *name, long value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  p_longs.insert(std::pair<std::string, long>(str,value));
}

void gridpack::component::DataCollection::addValue(char *name, bool value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  p_bools.insert(std::pair<std::string, bool>(str,value));
}

void gridpack::component::DataCollection::addValue(char *name, char *value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::string val = value;
  p_strings.insert(std::pair<std::string, std::string>(str,value));
}

void gridpack::component::DataCollection::addValue(char *name, float value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  p_floats.insert(std::pair<std::string, float>(str,value));
}

void gridpack::component::DataCollection::addValue(char *name, double value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  p_doubles.insert(std::pair<std::string, double>(str,value));
}

void gridpack::component::DataCollection::addValue(char *name,
    gridpack::ComplexType value, int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  p_complexType.insert(std::pair<std::string, gridpack::ComplexType>(str,value));
}

/**
 *  Modify current value of existing data element in
 *  DataCollection object
 *  @param name name of data element
 *  @param value new value of data element
 *  @return false if no element of the correct name and type exists in
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
bool gridpack::component::DataCollection::setValue(char *name, int value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, int>::iterator it;
  it = p_ints.find(str);
  if (it != p_ints.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::setValue(char *name, long value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, long>::iterator it;
  it = p_longs.find(str);
  if (it != p_longs.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::setValue(char *name, bool value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, bool>::iterator it;
  it = p_bools.find(str);
  if (it != p_bools.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::setValue(char *name, char *value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
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

bool gridpack::component::DataCollection::setValue(char *name, float value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, float>::iterator it;
  it = p_floats.find(str);
  if (it != p_floats.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::setValue(char *name, double value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, double>::iterator it;
  it = p_doubles.find(str);
  if (it != p_doubles.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::setValue(char *name,
    gridpack::ComplexType value, int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, gridpack::ComplexType>::iterator it;
  it = p_complexType.find(str);
  if (it != p_complexType.end()) {
    it->second = value;
    return true;
  } else {
    return false;
  }
}

/**
 *  Retrieve current value of existing data element in
 *  DataCollection object
 *  @param name name of data element
 *  @param value current value of data element
 *  @return false if no element of the correct name and type exists in
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
bool gridpack::component::DataCollection::getValue(char *name, int *value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, int>::iterator it;
  it = p_ints.find(str);
  if (it != p_ints.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::getValue(char *name, long *value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, long>::iterator it;
  it = p_longs.find(str);
  if (it != p_longs.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::getValue(char *name, bool *value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, bool>::iterator it;
  it = p_bools.find(str);
  if (it != p_bools.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::getValue(char *name, std::string *value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, std::string>::iterator it;
  it = p_strings.find(str);
  if (it != p_strings.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::getValue(char *name, float *value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, float>::iterator it;
  it = p_floats.find(str);
  if (it != p_floats.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::getValue(char *name, double *value,
    int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, double>::iterator it;
  it = p_doubles.find(str);
  if (it != p_doubles.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

bool gridpack::component::DataCollection::getValue(char *name,
    gridpack::ComplexType *value, int idx)
{
  std::string str = name;
  str.append(":");
  char buf[8];
  sprintf(buf,"%d",idx);
  str.append(buf);
  std::map<std::string, gridpack::ComplexType>::iterator it;
  it = p_complexType.find(str);
  if (it != p_complexType.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

/**
 * Dump contents of data collection to standard out
 */
void gridpack::component::DataCollection::dump(void)
{
  // print out integers
  std::map<std::string, int>::iterator int_it;
  int_it = p_ints.begin();
  int i = 0;
  while (int_it != p_ints.end()) {
    std::string key = int_it->first;
    int ival = int_it->second;
    cout << "  (INTEGER) key: "<<key<<" value: "<<ival<<endl;
    int_it++;
  }
  // print out longs
  std::map<std::string, long>::iterator long_it;
  long_it = p_longs.begin();
  while (long_it != p_longs.end()) {
    std::string key = long_it->first;
    long lval = long_it->second;
    cout << "  (LONG) key: "<<key<<" value: "<<lval<<endl;
    long_it++;
  }
  // print out bools
  std::map<std::string, bool>::iterator bool_it;
  bool_it = p_bools.begin();
  while (bool_it != p_bools.end()) {
    std::string key = bool_it->first;
    bool bval = bool_it->second;
    cout << "  (BOOL) key: "<<key<<" value: "<<bval<<endl;
    bool_it++;
  }
  // print out strings
  std::map<std::string, std::string>::iterator str_it;
  str_it = p_strings.begin();
  while (str_it != p_strings.end()) {
    std::string key = str_it->first;
    std::string sval = str_it->second;
    cout << "  (STRING) key: "<<key<<" value: "<<sval<<endl;
    str_it++;
  }
  // print out floats
  std::map<std::string, float>::iterator flt_it;
  flt_it = p_floats.begin();
  while (flt_it != p_floats.end()) {
    std::string key = flt_it->first;
    float fval = flt_it->second;
    cout << "  (FLOAT) key: "<<key<<" value: "<<fval<<endl;
    flt_it++;
  }
  // print out doubles
  std::map<std::string, double>::iterator dbl_it;
  dbl_it = p_doubles.begin();
  while (dbl_it != p_doubles.end()) {
    std::string key = dbl_it->first;
    double dval = dbl_it->second;
    cout << "  (DOUBLE) key: "<<key<<" value: "<<dval<<endl;
    dbl_it++;
  }
  // print out complex
  std::map<std::string, ComplexType>::iterator cmplx_it;
  cmplx_it = p_complexType.begin();
  while (cmplx_it != p_complexType.end()) {
    std::string key = cmplx_it->first;
    ComplexType cval = cmplx_it->second;
    cout << "  (COMPLEX) key: "<<key<<" value: "<<cval<<endl;
    cmplx_it++;
  }
}
