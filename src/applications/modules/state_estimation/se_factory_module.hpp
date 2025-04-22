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
 *         Adding more functions to handle bad data detection and more comprehensive outputs
 * @date   2025-04-20
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _se_factory_module_h_
#define _se_factory_module_h_

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/applications/components/se_matrix/se_components.hpp"

namespace gridpack {
namespace state_estimation {

// State estimation factory class
class SEFactoryModule
  : public gridpack::factory::BaseFactory<SENetwork> {
  public:
    typedef boost::shared_ptr<SENetwork> NetworkPtr;

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
     * Distribute measurements
     * @param measurements a vector containing all measurements
     */
    void setMeasurements(
        std::vector<gridpack::state_estimation::Measurement> measurements);

    /**
     * Initialize some parameters in state estimation components
     */
    void configureSE(void);

    /**
     * Identify bad data based on residual index
     * @param idx: index of maximum normalized residual
     */
    void identifyBadData(int idx);

    /**
     * Set mode for all components in the network
     * @param mode: operation mode for components
     */
    void setMode(int mode);

    /**
     * Report measurement details based on index
     * @param idx: index of measurement
     * @param details: string to store measurement details
     * @return true if measurement found, false otherwise
     */
    bool reportMeasurement(int idx, std::string& details);

    /**
     * Set voltage limits for all buses
     * @param v_min: minimum voltage magnitude
     * @param v_max: maximum voltage magnitude
     * @param enforce: whether to enforce limits
     */
    void setVoltageLimits(double v_min, double v_max, bool enforce);

    /**
     * Get mode name as string for reporting
     * @param mode: mode enum value
     * @return string representation of mode
     */
    const char* getModeNameString(int mode) const;

  private:
    NetworkPtr p_network;
};

} // state_estimation
} // gridpack
#endif // _se_factory_module_h_
