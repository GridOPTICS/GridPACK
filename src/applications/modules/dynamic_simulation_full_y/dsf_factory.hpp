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
    void setEvent(const DSFullBranch::Event &event);

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
