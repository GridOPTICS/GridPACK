/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   kds_factory_module.hpp
 * @author Da Meng and Yousu Chen 
 * @date   1/06/2015
 * 
 * @brief  
 * 
 * Modified by Xinya Li, July 2015
 */
// -------------------------------------------------------------

#ifndef _kalman_factory_h_
#define _kalman_factory_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/applications/components/kds_matrix/kds_components.hpp"


namespace gridpack {
namespace kalman_filter {

typedef gridpack::network::BaseNetwork<gridpack::kalman_filter::KalmanBus,
        gridpack::kalman_filter::KalmanBranch > KalmanNetwork;

class KalmanFactory
  : public gridpack::factory::BaseFactory<KalmanNetwork> {
  public:
    /**
     * Basic constructor
     * @param network: network associated with factory
     */
    KalmanFactory(NetworkPtr network);

    /**
     * Basic destructor
     */
    ~KalmanFactory();

    /**
     * Create the admittance (Y-Bus) matrix
     */
    void setYBus(void);

    /**
     * Get the updating factor for posfy11 stage ybus
     */
    gridpack::ComplexType setFactor(int sw2_2, int sw3_2);

    /**
     * Apply an event to all branches in the system
     * @param event a struct describing a fault
     */
    void setEvent(const KalmanBranch::Event &event);

    /**
     * Set the number of ensembles that will be used
     * @param nsize number of ensembles
     */
    void setEnsembleSize(int nsize);

    /**
     * Set the distribution width
     * @param sigma width of guassian distribution
     * @param noise width of guassian distribution
     */
    void setGaussianWidth(double sigma, double noise);

    /**
     * Set the current time step
     * @param istep current time step
     */
    void setCurrentTimeStep(int istep);

    /**
     * Create ensemble of initial conditions
     */
    void createEnsemble(void);

    /**
     * Evaluate X2 ensemble
     */
    void evaluateX2(void);

    /**
     * Evaluate X3 ensemble
     */
    void evaluateX3(void);

    /**
     * Test if each processor has at least one generator
     */
    bool checkGenerators();

  private:

    NetworkPtr p_network;
};

} // kalman_filter
} // gridpack
#endif
