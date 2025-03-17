/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_factory_module.hpp
 * @author Yousu Chen, Bruce Palmer 
 * @date   1/23/2015
 * 
 * @brief  
 * @update Yousu Chen
 *         Adding functions of applying virtual measurements
 * @date   2025-03-05
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _se_factory_module_h_
#define _se_factory_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/applications/components/se_matrix/se_components.hpp"

namespace gridpack {
namespace state_estimation {

class SEFactoryModule
  : public gridpack::factory::BaseFactory<SENetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    SEFactoryModule(NetworkPtr network);

    /**
     * Basic destructor
     */
    ~SEFactoryModule();

    /**
     * Create the admittance (Y-Bus) matrix
     */
    void setYBus(void);

    /**
     * Disribute measurements
     * @param measurements a vector containing all measurements
     */
    void setMeasurements(std::vector<Measurement> measurements);

    /**
     * Initialize some parameters in state estimation components
     */
    void configureSE(void);

    void identifyBadData(int idx);

    /**
     * Set the mode for all network components 
     * @param mode enumerated constant for different modes
     */
    void setMode(int mode);

    /**
     * Report Measurement details 
     * @Para: int idx for the measurement index 
     * std::string& details to return the measurement details.
     */
    //void reportMeasurement(int idx);
    bool reportMeasurement(int idx, std::string& details);
    //bool reportMeasurement(int busID, const std::string& type, std::string& details, int& globalIdx);

    /**
     * Set voltage limits for all buses
     * @param v_min minimum voltage limit
     * @param v_max maximum voltage limit
     * @param enforce whether to enforce limits
     */
    void setVoltageLimits(double v_min, double v_max, bool enforce);

  private:

    NetworkPtr p_network;
};

} // state_estimation
} // gridpack
#endif
