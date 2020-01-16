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
     * Struct for storing information on violations in contingency calculations
     */
    struct Violation {
      bool line_violation; /* branch overload violation */
      bool bus_violation;  /* bus voltage violation */
      int bus1;  /* bus ID for voltage violations or from bus for branch violations */
      int bus2;  /* to bus ID for branch violations */
      char tag[3]; /* 2 character identifier */
    };
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
     * Set voltage limits on all buses
     * @param Vmin lower bound on voltages
     * @param Vmax upper bound on voltages
     */
    void setVoltageLimits(double Vmin, double Vmax);

    /**
     * Check to see if there are any voltage violations in the network
     * @param area only check for voltage violations in this area
     * @return true if no violations found
     */
    bool checkVoltageViolations();
    bool checkVoltageViolations(int area);

    /**
     * Check to see if there are any Q limit violations in the network
     * @param area only check for Q limit violations in this area
     * @return true if no violations found
     */
    bool checkQlimViolations();
    bool checkQlimViolations(int area);

    /**
     * Clear changes that were made for Q limit violations and reset
     * system to its original state
     */
    void clearQlimViolations();

    /**
     * Set "ignore" parameter on all buses with violations so that subsequent
     * checks are not counted as violations
     */
    void ignoreVoltageViolations();

    /**
     * Clear "ignore" parameter on all buses
     */
    void clearVoltageViolations();

    /**
     * Check to see if there are any line overload violations in
     * the network. The last call checks for overloads on specific lines.
     * @param area only check for voltage violations in this area
     * @param bus1 original index of "from" bus for branch
     * @param bus2 original index of "to" bus for branch
     * @param tags line IDs for individual lines
     * @param violations true if violation detected on branch, false otherwise
     * @return true if no violations found
     */
    bool checkLineOverloadViolations();
    bool checkLineOverloadViolations(int area);
    bool checkLineOverloadViolations(std::vector<int> &bus1, std::vector<int> &bus2,
        std::vector<std::string> &tags, std::vector<bool> &violations);

    /**
     * Set "ignore" parameter on all lines with violations so that subsequent
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

    /**
     * Scale generator real power. If zone less than 1 then scale all
     * generators in the area
     * @param scale factor to scale real power generation
     * @param area index of area for scaling generation
     * @param zone index of zone for scaling generation
     */
    void scaleGeneratorRealPower(double scale, int area, int zone);

    /**
     * Scale load real power. If zone less than 1 then scale all
     * loads in the area
     * @param scale factor to scale load real power
     * @param area index of area for scaling load
     * @param zone index of zone for scaling load
     * @return false if there is not enough capacity to change generation
     *         by requested amount
     */
    void scaleLoadPower(double scale, int area, int zone);

    /**
     * Return the total real power load for all loads in the zone. If zone
     * less than 1, then return the total load for the area
     * @param area index of area
     * @param zone index of zone
     * @return total load
     */
    double getTotalLoadRealPower(int area, int zone);

    /**
     * Return the current real power generation and the maximum and minimum total
     * power generation for all generators in the zone. If zone is less than 1
     * then return values for all generators in the area
     * @param area index of area
     * @param zone index of zone
     * @param total total real power generation
     * @param pmin minimum allowable real power generation
     * @param pmax maximum available real power generation
     */
    void getGeneratorMargins(int area, int zone, double *total, double *pmin,
        double *pmax);

    /**
     * Reset power of loads and generators to original values
     */
    void resetPower();

    /**
     * Set parameters for real time path rating diagnostics
     * @param src_area generation area
     * @param src_zone generation zone
     * @param load_area load area
     * @param load_zone load zone
     * @param gen_scale scale factor for generation
     * @param load_scale scale factor for loads
     */
    void setRTPRParams(int src_area, int src_zone, int load_area,
            int load_zone, double gen_scale, double load_scale);

    /**
     * Return vector describing all violations
     * @return violation vector
     */
    std::vector<Violation> getViolations();

    /**
     * Clear violation vector
     */
    void clearViolations();

    /**
     * User rate B parameter for line overload violations
     * @param flag if true, use RATEB parameter
     */
    void useRateB(bool flag);

  private:

    NetworkPtr p_network;
    std::vector<bool> p_saveIsolatedStatus;

    std::vector<Violation> p_violations;

    bool p_rateB;
};

} // powerflow
} // gridpack
#endif
