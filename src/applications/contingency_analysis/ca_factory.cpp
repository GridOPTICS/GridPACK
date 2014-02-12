/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_factory.cpp
 * @author Yousu Chen 
 * @date   Feb 11, 2014 
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "ca_components.hpp"
#include "ca_factory.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"


namespace gridpack {
namespace contingency_analysis {

// Powerflow factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
CAFactory::CAFactory(CAFactory::NetworkPtr network)
  : gridpack::factory::BaseFactory<CANetwork>(network)
{
  p_network = network;
}

/**
 * Basic destructor
 */
gridpack::contingency_analysis::CAFactory::~CAFactory()
{
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::contingency_analysis::CAFactory::setYBus(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // Invoke setYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
    dynamic_cast<CABranch*>(p_network->getBranch(i).get())->setYBus();
  }

  // Invoke setYBus method on all bus objects
  for (i=0; i<numBus; i++) {
    dynamic_cast<CABus*>(p_network->getBus(i).get())->setYBus();
  }

}

/**
 * Find GBus vector 
 */
void gridpack::contingency_analysis::CAFactory::setGBus(void)
{
  int numBus = p_network->numBuses();
  int i;

  // Invoke setGBus method on all bus objects
  for (i=0; i<numBus; i++) {
    dynamic_cast<CABus*>(p_network->getBus(i).get())->setGBus();
  }
}

/**
  * Make SBus vector 
  */
void gridpack::contingency_analysis::CAFactory::setSBus(void)
{
  int numBus = p_network->numBuses();
  int i;

  // Invoke setSBus method on all bus objects
  for (i=0; i<numBus; i++) {
    dynamic_cast<CABus*>(p_network->getBus(i).get())->setSBus();
  }
}

/**
  * Create the PQ 
  */
void gridpack::contingency_analysis::CAFactory::setPQ(void)
{
  int numBus = p_network->numBuses();
  int i;
  ComplexType values[2];

  for (i=0; i<numBus; i++) {
    dynamic_cast<CABus*>(p_network->getBus(i).get())->vectorValues(values);
  }
}

/**
 * Create the Jacobian matrix
 */
void gridpack::contingency_analysis::CAFactory::setJacobian(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  /*for (i=0; i<numBus; i++) {
    dynamic_cast<CABus*>(p_network->getBus(i).get())->setYBus();
  }

  for (i=0; i<numBranch; i++) {
    dynamic_cast<CABranch*>(p_network->getBranch(i).get())->getJacobian();
  }*/
}

/**
 * Get values from the admittance (Y-Bus) matrix
 */
#if 0
gridpack::ComplexType gridpack::contingency_analysis::CAFactory::calMis(
    gridpack::math::Vector V,
    gridpack::math::Vector SBUS)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;
  gridpack::ComplexType ibus;
  gridpack::ComplexType mis;

  // MIS = V * conj (YBus * V) - SBUS

  // Invoke getYBus method on all bus objects
  /* for (i=0; i<numBus; i++) {
     ibus =
     (dynamic_cast<gridpack::contingency_analysis::CABus*>(p_network->getBus(i).get()))->getYBus()
   * V(i) ;
   mis(i) = V * conj(ibus
   }

  // Invoke getYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
  (dynamic_cast<gridpack::contingency_analysis::CABranch*>(p_network->getBranch(i).get()))->getYBus();
  }*/
}
#endif


} // namespace contingency_analysis 
} // namespace gridpack
