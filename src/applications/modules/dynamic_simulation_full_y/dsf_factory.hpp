/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsf_factory.hpp
 * @author Shuangshuang Jin 
 * @date   Feb 04, 2015
 * @last modified date   May 13, 2015
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _dsf_factory_h_
#define _dsf_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "dsf_components.hpp"
#include <vector>

namespace gridpack {
namespace dynamic_simulation {

class DSFullFactory
  : public gridpack::factory::BaseFactory<DSFullNetwork> {
  public:

    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    DSFullFactory(NetworkPtr network);

    /**
     * Basic destructor
     */
    ~DSFullFactory();

    /**
     * Create the admittance (Y-Bus) matrix
     */
    void setYBus(void);

    /**
     * Create the admittance (Y-Bus) matrix for Branch only
     */
    void setYBranch(void);

    /**
     * Get the updating factor for posfy11 stage ybus
     */
    gridpack::ComplexType setFactor(int sw2_2, int sw3_2);

    /**
     * Apply an event to all branches in the system
     * @param event a struct describing a fault
     */
    void setEvent(const Event &event);

    /**
     * Check network to see if there is a process with no generators
     * @return true if all processors have at least on generator
     */
    bool checkGen(void);

    /**
     * Initialize init vectors for integration
     * @param ts time step
     */
    void initDSVect(double ts);

    /**
     * Update vectors in each integration time step (Predictor)
     */
    void predictor_currentInjection(bool flag);

    /**
     * Update vectors in each integration time step (Corrector)
     */
    void corrector_currentInjection(bool flag);

    /**
     * Update vectors in each integration time step (Predictor)
     */
    void predictor(double t_inc, bool flag);

    /**
     * Update vectors in each integration time step (Corrector)
     */
    void corrector(double t_inc, bool flag);
	
	/**
     * Update dynamic load internal relays action
     */
	void dynamicload_post_process(double t_inc, bool flag);
	
	/**
     * load parameters for the extended buses from composite load model
     */
	void LoadExtendedCmplBus( );
	
	/**
     * set voltage for the extended buses from composite load model
     */
	void setExtendedCmplBusVoltage( );

    /**
     * Set volt from volt_full
     */
    void setVolt(bool flag);
	
	/**
    * update bus frequecy
    */
    void updateBusFreq(double delta_t);
	
	std::vector<double> grabWideAreaFreq() ;  //renke hard coded
	void setWideAreaFreqforPSS(double freq); //renke hard coded
	
   /**
     * update bus relay status
     */
    bool updateBusRelay(bool flag,double delta_t);
	
	/**
     * update branch relay status
     */
    bool updateBranchRelay(bool flag,double delta_t);
	
	/**
     * update old bus voltage
     */
	void updateoldbusvoltage();
	void printallbusvoltage();

    /**
     * Add constant impedance load admittance to diagonal elements of
     * Y-matrix
     */
    void addLoadAdmittance();

    bool securityCheck();

    /**
     * Scale generator real power. If zone less than 1 then scale all
     * generators in the area
     * @param scale factor to scale real power generation
     * @param area index of area for scaling generation
     * @param zone index of zone for scaling generation
     */
    void scaleGeneratorRealPower(double scale, int area, int zone);

    /**
     * Scale load power. If zone less than 1 then scale all
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

#ifdef USE_FNCS
    /**
     * Scatter load from FNCS framework to buses
     */
    void scatterLoad();

    /**
     * Gather voltage magnitude and phase angle to root node and transmit to
     * FNCS framework
     */
    void gatherVoltage();
#endif

  private:

    NetworkPtr p_network;
#ifdef USE_FNCS
    gridpack::hash_distr::HashDistribution<DSFullNetwork,
           gridpack::ComplexType,gridpack::ComplexType> *p_hash;

    gridpack::serial_io::SerialBusIO<DSFullNetwork> *p_busIO;
#endif

    int p_numBus;

    DSFullBus **p_buses;

    int p_numBranch;

    DSFullBranch **p_branches;
};

} // dynamic_simulation
} // gridpack
#endif
