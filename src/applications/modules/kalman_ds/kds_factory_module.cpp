/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   kds_factory_module.cpp
 * @author Da Meng and Yousu Chen
 * @date   1/06/2015
 * 
 * @brief  
 * 
 * Modified by Xinya Li, July 2015 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "kds_factory_module.hpp"

namespace gridpack {
namespace kalman_filter {

// State estimation factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
KalmanFactory::KalmanFactory(KalmanFactory::NetworkPtr network)
  : gridpack::factory::BaseFactory<KalmanNetwork>(network)
{
  p_network = network;
}

/**
 * Basic destructor
 */
gridpack::kalman_filter::KalmanFactory::~KalmanFactory()
{
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::kalman_filter::KalmanFactory::setYBus(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // Invoke setYBus method on all bus objects
  for (i=0; i<numBus; i++) {
    (dynamic_cast<KalmanBus*>(p_network->getBus(i).get()))->setYBus();
  }

  // Invoke setYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
    (dynamic_cast<KalmanBranch*>(p_network->getBranch(i).get()))->setYBus();
  }
}

/**
 * Get the updating factor for posfy11 stage ybus
 */
gridpack::ComplexType
gridpack::kalman_filter::KalmanFactory::setFactor(int sw2_2, int sw3_2)
{
  gridpack::ComplexType dummy(-999.0, -999.0);

  int numBranch = p_network->numBranches();
  int i;

  // Invoke getPosfy11YbusUpdateFactor method on all branch objects
  for (i=0; i<numBranch; i++) {
    gridpack::ComplexType ret = (dynamic_cast<KalmanBranch*>(p_network->getBranch(i).get()))
      ->getPosfy11YbusUpdateFactor(sw2_2, sw3_2);
    if (ret != dummy) {
      return ret;
    }
  }
  return dummy;
}


/**
 * Apply an event to all branches in the system
 * @param event a struct describing a fault
 */
void gridpack::kalman_filter::KalmanFactory::setEvent(const
    gridpack::kalman_filter::KalmanBranch::Event &event)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;
  for (i=0; i<numBus; i++) {
    dynamic_cast<KalmanBus*>(p_network->getBus(i).get())->clearEvent();
  }
  for (i=0; i<numBranch; i++) {
    dynamic_cast<KalmanBranch*>(p_network->getBranch(i).get())->setEvent(event);
  }
}


/**
 * Set the number of ensembles that will be used
 * @param nsize number of ensembles
 */
void gridpack::kalman_filter::KalmanFactory::setEnsembleSize(int nsize)
{
  int numBus = p_network->numBuses();
  int i;

  // Invoke setEnsembleSize method on all bus objects
  for (i=0; i<numBus; i++) {
    (dynamic_cast<KalmanBus*>(p_network->getBus(i).get()))->
      setEnsembleSize(nsize);
  }
}

/**
 * Set the distribution width
 * @param sigma width of guassian distribution
 * @param noise width of guassian distribution
 */
void gridpack::kalman_filter::KalmanFactory::setGaussianWidth(double sigma,
    double noise)
{
  int numBus = p_network->numBuses();
  int i;

  // Invoke setGaussianWidth method on all bus objects
  for (i=0; i<numBus; i++) {
    (dynamic_cast<KalmanBus*>(p_network->getBus(i).get()))->
      setGaussianWidth(sigma,noise);
  }
}

/**
 * Set the current time step
 * @param istep current time step
 */
void gridpack::kalman_filter::KalmanFactory::setCurrentTimeStep(int istep)
{
  int numBus = p_network->numBuses();
  int i;

  // Invoke setCurrentTimeStep on all bus objects
  for (i=0; i<numBus; i++) {
    (dynamic_cast<KalmanBus*>(p_network->getBus(i).get()))->
      setCurrentTimeStep(istep);
  }
}

/**
 * Create ensemble of initial conditions
 */
void gridpack::kalman_filter::KalmanFactory::createEnsemble(void)
{
  int numBus = p_network->numBuses();
  int i;

  // Invoke createEnsemble on all bus objects
  for (i=0; i<numBus; i++) {
    (dynamic_cast<KalmanBus*>(p_network->getBus(i).get()))->
      createEnsemble();
  }
}


/**
 * Create evaluate X2 distribution
 */
void gridpack::kalman_filter::KalmanFactory::evaluateX2(void)
{
  int numBus = p_network->numBuses();
  int i;

  // Invoke evaluateX2 on all bus objects
  for (i=0; i<numBus; i++) {
    (dynamic_cast<KalmanBus*>(p_network->getBus(i).get()))->
      evaluateX2();
  }
}

/**
 * Create evaluate X3 distribution
 */
void gridpack::kalman_filter::KalmanFactory::evaluateX3(void)
{
  int numBus = p_network->numBuses();
  int i;

  // Invoke evaluateX3 on all bus objects
  for (i=0; i<numBus; i++) {
    (dynamic_cast<KalmanBus*>(p_network->getBus(i).get()))->
      evaluateX3();
  }
}

/**
 * Test if each processor has at least one generator
 */
bool gridpack::kalman_filter::KalmanFactory::checkGenerators()
{
  int size = p_network->communicator().size();
  int rank = p_network->communicator().rank();
  std::vector<int> ngen(size);
  int numBus = p_network->numBuses();
  int i;
  for (i=0; i<size; i++) ngen[i] = 0;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      ngen[rank] += (dynamic_cast<KalmanBus*>(p_network->getBus(i).get()))->
        numGenerators();
    }
  }
  p_network->communicator().sum(&ngen[0],size);
  bool ok = true;
  int tgen = 0;
  for (i=0; i<size; i++) {
    if (ngen[i] == 0) ok = false;
    tgen += ngen[i];
  }
  if (rank == 0) {
    printf ("Total number of generators: %d\n",tgen);
    for (i=0; i<size; i++) {
      if (i<10) {
        printf("  Process[%1d]:   %6d Generators\n",i,ngen[i]);
      } else if (i<100) {
        printf("  Process[%2d]:  %6d Generators\n",i,ngen[i]);
      } else {
        printf("  Process[%2d]: %6d Generators\n",i,ngen[i]);
      }
    }
  }
  return ok;
}

} // namespace kalman_filter
} // namespace gridpack
