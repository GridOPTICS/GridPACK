/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_factory_module.hpp
 * @author Bruce Palmer
 * @date   2014-01-28 11:33:42 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _pf_factory_module_h_
#define _pf_factory_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/applications/components/pf_matrix/pf_components.hpp"

namespace gridpack {
namespace powerflow {

/// The type of network used in the powerflow application
typedef gridpack::network::BaseNetwork<PFBus, PFBranch > PFNetwork;

class PFFactoryModule
  : public gridpack::factory::BaseFactory<PFNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    PFFactoryModule(NetworkPtr network);

    /**
     * Basic destructor
     */
    ~PFFactoryModule();

    /**
     * Create the admittance (Y-Bus) matrix
     */
    void setYBus(void);

    /**
     * Make SBus vector 
     */
    void setSBus(void);

    /**
      * Update pg of specified bus element based on their genID
      * @param name 
      * @param busID
      * @param genID
      * @param value
      */
//    void updatePg(std::string &name, int busID, std::string genID, double value);
    void updatePg(int busID, std::string genID, double value);
    void updateQg(int busID, std::string genID, double value);

    /**
     * Create the PQ 
     */
    void setPQ(void);

    /**
     * Check for lone buses in the system. Do this by looking for buses that
     * have no branches attached to them or for whom all the branches attached
     * to the bus have all transmission elements with status false (the element
     * is off). Set status of bus to isolated so that it does not contribute to
     * powerflow matrix
     * @param stream optional stream pointer that can be used to print out IDs
     * of isolated buses
     * @return false if there is an isolated bus in the network
     */
    bool checkLoneBus(std::ofstream *stream = NULL);

    /**
     * Set lone buses back to their original status.
     */
    void clearLoneBus();

    /**
     * Check to see if there are any voltage violations in the network
     * @param minV maximum voltage limit
     * @param maxV maximum voltage limit
     * @param area only check for voltage violations in this area
     * @return true if no violations found
     */
    bool checkVoltageViolations(double Vmin, double Vmax);
    bool checkVoltageViolations(int area, double Vmin, double Vmax);

    /**
     * Set "ignore" parameter on all buses with violations so that subsequent
     * checks are not counted as violations
     * @param minV maximum voltage limit
     * @param maxV maximum voltage limit
     */
    void ignoreVoltageViolations(double Vmin, double Vmax);

    /**
     * Clear "ignore" parameter on all buses
     */
    void clearVoltageViolations();

    /**
     * Check to see if there are any line overload violations in the network
     * @param area only check for voltage violations in this area
     * @return true if no violations found
     */
    bool checkLineOverloadViolations();
    bool checkLineOverloadViolations(int area);

    /**
     * Set "ignore" paramter on all lines with violations so that subsequent
     * checks are not counted as violations
     */
    void ignoreLineOverloadViolations();

    /**
     * Clear "ignore" parameter on all lines
     */
    void clearLineOverloadViolations();

    /**
     * Reinitialize voltages
     */
    void resetVoltages();
  private:

    NetworkPtr p_network;
    std::vector<bool> p_saveIsolatedStatus;
};

} // powerflow
} // gridpack
#endif
