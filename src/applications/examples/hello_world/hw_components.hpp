/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hw_components.hpp
 * @author Bruce Palmer
 * @date   2013-10-24 14:30:43 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _hw_components_h_
#define _hw_components_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/network/base_network.hpp"

namespace gridpack {
namespace hello_world {

class HWBus
  : public gridpack::component::BaseBusComponent {
  public:
    /**
     *  Simple constructor
     */
    HWBus(void);

    /**
     *  Simple destructor
     */
    ~HWBus(void);

    /**
     * Load values stored in DataCollection object into HWBus object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       bus that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);

    /**
     * Write output from buses to standard out
     * @param string (output) string with information to be printed out
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if bus is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const char *signal = NULL);

  private:
    int p_original_idx;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<gridpack::component::BaseBusComponent>(*this)
      & p_original_idx;
  }  

};

class HWBranch
  : public gridpack::component::BaseBranchComponent {
  public:
    /**
     *  Simple constructor
     */
    HWBranch(void);

    /**
     *  Simple destructor
     */
    ~HWBranch(void);

    /**
     * Load values stored in DataCollection object into HWBranch object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       branch that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);

    /**
     * Write output from branches to standard out
     * @param string (output) string with information to be printed out
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if branch is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const char *signal = NULL);

  private:
    int p_original_idx1;
    int p_original_idx2;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<gridpack::component::BaseBranchComponent>(*this)
      & p_original_idx1
      & p_original_idx2;
  }  

};


/// The type of network used in the examples application
typedef network::BaseNetwork<HWBus, HWBranch > HWNetwork;


}     // hello_world
}     // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::hello_world::HWBus);
BOOST_CLASS_EXPORT_KEY(gridpack::hello_world::HWBranch);


#endif
