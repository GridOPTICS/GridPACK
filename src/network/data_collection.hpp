#ifndef _data_collection_h
#define _data_collection_h

#include <map>
// Simple outline of data collection object
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
   *  Add variables to DataCollection object
   *  @param name: name given to data element
   *  @param value: value of data element
   */
  void AddValue(char *name, int value);
  void AddValue(char *name, long value);
  void AddValue(char *name, bool value);
  void AddValue(char *name, char *value);
  void AddValue(char *name, float value);
  void AddValue(char *name, double value);
  void AddValue(char *name, SingleComplex value);
  void AddValue(char *name, DoubleComplex value);

  /**
   *  Modify current value of existing data element in
   *  DataCollection object
   *  @param name: name of data element
   *  @param value: new value of data element
   *  @return: false if no element of the correct name and type exists in
   *  DataCollection object
   */
  bool SetValue(char *name, int value);
  bool SetValue(char *name, long value);
  bool SetValue(char *name, bool value);
  bool SetValue(char *name, char *value);
  bool SetValue(char *name, float value);
  bool SetValue(char *name, double value);
  bool SetValue(char *name, SingleComplex value);
  bool SetValue(char *name, DoubleComplex value);

  /**
   *  Retrieve current value of existing data element in
   *  DataCollection object
   *  @param name: name of data element
   *  @param value: current value of data element
   *  @return: false if no element of the correct name and type exists in
   *  DataCollection object
   */
  bool GetValue(char *name, int *value);
  bool GetValue(char *name, long *value);
  bool GetValue(char *name, bool *value);
  bool GetValue(char *name, char *std::string);
  bool GetValue(char *name, float *value);
  bool GetValue(char *name, double *value);
  bool GetValue(char *name, SingleComplex *value);
  bool GetValue(char *name, DoubleComplex *value);
private:
  std::map<std::string, int> p_ints; 
  std::map<std::string, long> p_longs; 
  std::map<std::string, bool> p_bools; 
  std::map<std::string, std::string> p_strings; 
  std::map<std::string, float> p_floats; 
  std::map<std::string, double> p_doubles; 
  std::map<std::string, SingleComplex> p_singleComplex; 
  std::map<std::string, DoubleComplex> p_doubleComplex; 
}

#endif // _data_collection_h
