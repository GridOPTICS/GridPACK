/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#ifndef _base_analytics_intfc_h
#define _base_analytics_intfc_h

#include <map>
#include <string>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#include "boost/smart_ptr/shared_ptr.hpp"
//#include <boost/serialization/export.hpp>

#include "gridpack/component/data_collection.hpp"

// Simple outline of data collection object

namespace gridpack{
namespace component{

typedef boost::shared_ptr<DataCollection> DataPtr;

class BaseBusAnalyticsInterface {
public:
  /**
   * Simple constructor
   */
  BaseBusAnalyticsInterface(void);

  /**
   * Simple destructor
   */
  virtual ~BaseBusAnalyticsInterface(void);

  /**
   * Set data collection object inside interface
   * @param data_collection pointer to data collection object
   */
  virtual void setData(DataPtr &data);

  /**
   * return number of generators on bus
   */
  virtual int numGenerators();
  
  /**
   * return number of storage units on bus
   */
  virtual int numStorage();
  
  /**
   * return number of loads on bus
   */
  virtual int numLoads();
  
private:

  DataCollection *p_data;

  friend class boost::serialization::access;

  /// Serialization method
  template<class Archive> void serialize(Archive &ar, const unsigned int version)
  {
    ar & p_data;
  }

};

//BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseBusAnalyticsInterface)

class BaseBranchAnalyticsInterface {
public:
  /**
   * Simple constructor
   */
  BaseBranchAnalyticsInterface(void);

  /**
   * Simple destructor
   */
  ~BaseBranchAnalyticsInterface(void);

  /**
   * Set data collection object inside interface
   * data_collection pointer to data collection object
   */
  virtual void setData(DataPtr &data);

  /**
   * return number of lines on branch
   */
  virtual int numLines();
  
private:

  DataCollection *p_data;

  friend class boost::serialization::access;

  /// Serialization method
  template<class Archive> void serialize(Archive &ar, const unsigned int version)
  {
    ar & p_data;
  }

};

//BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseBranchAnalyticsInterface)

}    // component
}    // gridpack


#endif // _base_analytics_intfc_h
