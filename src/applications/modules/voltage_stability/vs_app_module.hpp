/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   vs_app_module.hpp
 * @author Bruce Palmer
 * @date   2014-01-28 11:30:49 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _vs_app_module_h_
#define _vs_app_module_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/serial_io/serial_io.hpp"
#include "gridpack/configuration/configuration.hpp"
#include "vs_factory_module.hpp"
#include "gridpack/applications/modules/voltage_stability/vs_app_module.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"

namespace gridpack {
namespace voltage_stability {

// Structs that are used for some applications

struct pathVSBranch{
  int fromBus;
  int toBus;
  std::string branchID;
};

struct stressVSBus{
  int busID;
  std::string genID;
  double participation;
};

struct pathVSStress{
  std::string name;
  std::vector<pathVSBranch> path;
  std::vector<stressVSBus> sourceArea;
  std::vector<stressVSBus> sinkArea;
};

// Contingency types
enum VSContingencyType{VS_Generator, VS_Branch};

// Struct that is used to define a collection of contingencies

struct VSContingency
{
  int p_type;
  std::string p_name;
  // Line contingencies
  std::vector<int> p_from;
  std::vector<int> p_to;
  std::vector<std::string> p_ckt;
  // Status of line before contingency
  std::vector<bool> p_saveLineStatus;
  // Generator contingencies
  std::vector<int> p_busid;
  std::vector<std::string> p_genid;
  // Status of generator before contingency
  std::vector<bool> p_saveGenStatus;
};

// Calling program for powerflow application

class VSAppModule
{
  public:

    /**
     * Basic constructor
     */
    VSAppModule(void);

    /**
     * Basic destructor
     */
    ~VSAppModule(void);

    void IncrementGeneratorRealPower(double inc, int area, int zone, double gtotal);
    
    /**
     * Increment load power based off specified value. 
     * Increment loads in specified area.
     * @param transfer value to increment load real power
     * @param area index of area for incrementing load
     * @param zone index of zone for incrementing load
     * @param total active power demand of the area
     */
    void IncrementLoadPower(double inc, int area, int zone, double ltotal);
    
    /**
     * Return the total real power load for all loads in the zone. If zone
     * less than 1, then return the total load real power for the area
     * @param area index of area
     * @param zone index of zone
     * @return total load
     */
    double getTotalLoadRealPower(int area, int zone);
        
    /**
     * Check to see if PV Analysis is complete
     * @return return true if parameter is not found
     */
    bool isPVAnlyDone(double increment, double max_increment);
    
    /**
     * Set up PV Curve internal parameters and initialize
     */
    void InitializePVCurve(std::string filename);
    
    /**
     * Execute one transfer increment
     */
    void IncrementPVCurveStep();

  private:
    

    // pointer to network
    boost::shared_ptr<gridpack::powerflow::PFNetwork> p_network;

    // communicator for network
    gridpack::parallel::Communicator p_comm;

    // pointer to factory
    boost::shared_ptr<gridpack::powerflow::PFFactoryModule> p_factory;

    // pointer to factory
    boost::shared_ptr<gridpack::voltage_stability::VSFactoryModule> v_factory;
    
    // pointer to bus IO module
    boost::shared_ptr<gridpack::serial_io::SerialBusIO<gridpack::powerflow::PFNetwork> > p_busIO;

    // pointer to branch IO module
    boost::shared_ptr<gridpack::serial_io::SerialBranchIO<gridpack::powerflow::PFNetwork> > p_branchIO;

    // pointer to configuration module
    gridpack::utility::Configuration *p_config;
    // Variables for PV Analysis
    bool p_bPVAnlyDone,PV_header;
    double max_increment, increment, gen_scale, load_scale, zone, lt, gt;
    double current_increment;
    int sink_area, src_area;

#ifdef USE_GOSS
    gridpack::goss::GOSSClient p_goss_client;

    std::string p_simID;
#endif
};

} // voltage_stability
} // gridpack
#endif
