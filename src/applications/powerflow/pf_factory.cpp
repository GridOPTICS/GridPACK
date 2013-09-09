// -------------------------------------------------------------
/**
 * @file   pf_factory.cpp
 * @author Bruce Palmer
 * @date   2013-09-09 15:36:40 d3g096
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
#include "gridpack/applications/powerflow/pf_components.hpp"
#include "gridpack/applications/powerflow/pf_factory.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"


namespace gridpack {
namespace powerflow {

// Powerflow factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
PFFactory::PFFactory(PFFactory::NetworkPtr network)
  : gridpack::factory::BaseFactory<PFNetwork>(network)
{
  p_network = network;
}

/**
 * Basic destructor
 */
gridpack::powerflow::PFFactory::~PFFactory()
{
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::powerflow::PFFactory::setYBus(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // Invoke setYBus method on all bus objects
  for (i=0; i<numBus; i++) {
    dynamic_cast<PFBus*>(p_network->getBus(i).get())->setYBus();
  }

  // Invoke setYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
    dynamic_cast<PFBranch*>(p_network->getBranch(i).get())->setYBus();
  }
}

/**
 * Find GBus vector 
 */
void gridpack::powerflow::PFFactory::setGBus(void)
{
  int numBus = p_network->numBuses();
  int i;

  // Invoke setGBus method on all bus objects
  for (i=0; i<numBus; i++) {
    dynamic_cast<PFBus*>(p_network->getBus(i).get())->setGBus();
  }
}

/**
  * Make SBus vector 
  */
void gridpack::powerflow::PFFactory::setSBus(void)
{
  int numBus = p_network->numBuses();
  int i;

  // Invoke setSBus method on all bus objects
  for (i=0; i<numBus; i++) {
    dynamic_cast<PFBus*>(p_network->getBus(i).get())->setSBus();
  }
}

/**
  * Create the PQ 
  */
void gridpack::powerflow::PFFactory::setPQ(void)
{
  int numBus = p_network->numBuses();
  int i;
  ComplexType values[2];

  for (i=0; i<numBus; i++) {
    dynamic_cast<PFBus*>(p_network->getBus(i).get())->vectorValues(values);
  }
}

/**
 * Create the Jacobian matrix
 */
void gridpack::powerflow::PFFactory::setJacobian(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  /*for (i=0; i<numBus; i++) {
    dynamic_cast<PFBus*>(p_network->getBus(i).get())->setYBus();
  }

  for (i=0; i<numBranch; i++) {
    dynamic_cast<PFBranch*>(p_network->getBranch(i).get())->getJacobian();
  }*/
}

/**
 * Get values from the admittance (Y-Bus) matrix
 */
#if 0
gridpack::ComplexType gridpack::powerflow::PFFactory::calMis(
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
     (dynamic_cast<gridpack::powerflow::PFBus*>(p_network->getBus(i).get()))->getYBus()
   * V(i) ;
   mis(i) = V * conj(ibus
   }

  // Invoke getYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
  (dynamic_cast<gridpack::powerflow::PFBranch*>(p_network->getBranch(i).get()))->getYBus();
  }*/
}
#endif

// -------------------------------------------------------------
// PFFactory::operator()
// -------------------------------------------------------------
void
PFFactory::operator() (const math::Vector& X, math::Matrix& J)
{
  // In both the Netwon-Raphson and PETSc nonlinear solver (some
  // methods) implementations, the RHS function builder before this,
  // so we may be able to count on the current solution being pushed
  // back on the netork in the RHS function builder method

  // Push current values in X vector back into network components
  // Need to implement setValues method in PFBus class in order for this to
  // work
  // gridpack::mapper::BusVectorMap<PFNetwork> vMap(*p_network);
  // vMap.mapToBus(X);

  // Exchange data between ghost buses (I don't think we need to exchange data
  // between branches)
  // FIXME: called twice per iteration -- may be able to skip it here
  // network->updateBuses();
  
  // Set to build Jacobian
  setJacobian();

  // build the Jacobian
  gridpack::mapper::FullMatrixMap<PFNetwork> jMap(p_network);
  jMap.mapToMatrix(J);
}

void
PFFactory::operator() (const math::Vector& X, math::Vector& PQ)
{
  // In both the Netwon-Raphson and PETSc nonlinear solver
  // implementations, this is called before the Jacobian builder, so
  // we may only need to map the solution back on to the network here.

  // Push current values in X vector back into network components
  // Need to implement setValues method in PFBus class in order for this to
  // work

  // FIXME: set proper mode for state (voltage, phase)
  
  gridpack::mapper::BusVectorMap<PFNetwork> vMap(p_network);
  vMap.mapToBus(X);

  // Exchange data between ghost buses (I don't think we need to exchange data
  // between branches)
  p_network->updateBuses();
  
  // set to build RHS vector
  setPQ();

  // build the RHS vector
  vMap.mapToVector(PQ);
}


} // namespace powerflow
} // namespace gridpack
