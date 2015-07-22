/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   uc_components.hpp
 * @author 
 * @date   
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _uc_components_h_
#define _uc_components_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/network/base_network.hpp"

namespace gridpack {
namespace unit_commitment {

class UCBus
  : public gridpack::component::BaseBusComponent
 {
  public:
    /**
     *  Simple constructor
     */
    UCBus(void);

    /**
     *  Simple destructor
     */
    ~UCBus(void);

    /**
     * Load values stored in DataCollection object into UCBus object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       bus that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);
    /**
     * get objective function
    */
    double objectiveFunction(void);
    /** 
      * get solution
      */
    bool solution(void);

    int numGen;
    std::vector<double> p_iniLevel;
    std::vector<double> p_minUpTime;
    std::vector<double> p_minDownTime;
    std::vector<double> p_minPower;
    std::vector<double> p_maxPower;
    std::vector<double> p_costConst;
    std::vector<double> p_costLinear;
    std::vector<double> p_costQuad;
    std::vector<double> p_rampUp;
    std::vector<double> p_rampDown;
    std::vector<double> p_startUp;
    std::vector<double> p_initPeriod;
    std::vector<double> p_startCap;
    std::vector<double> p_shutCap;
    std::vector<double> p_opMaxGen;
    std::vector<double> p_powerProduced;
    std::vector<double> p_powerReserved;
    std::vector<int> p_gen;
  private:
    int p_num_generator;
    int p_bus_id;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<gridpack::component::BaseBusComponent>(*this)
      & p_num_generator;
  }  

};

class UCBranch
  : public gridpack::component::BaseBranchComponent {
  public:
    /**
     *  Simple constructor
     */
    UCBranch(void);

    /**
     *  Simple destructor
     */
    ~UCBranch(void);

    /**
     * Load values stored in DataCollection object into UCBranch object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       branch that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);

     

  private:


};


/// The type of network used in the examples application
typedef network::BaseNetwork<UCBus, UCBranch > UCNetwork;


}     // unit_commitment
}     // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::unit_commitment::UCBus);
BOOST_CLASS_EXPORT_KEY(gridpack::unit_commitment::UCBranch);


#endif
