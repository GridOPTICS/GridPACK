// -------------------------------------------------------------
/**
 * @file   base_component.cpp
 * @author Bruce Palmer
 * @date   May 13, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/parallel/distributed.hpp"
#include "gridpack/component/base_component.hpp"

// Base implementation of the MatVecInterface. These functions should be
// overwritten in actual components

/**
 * Constructor
 */
gridpack::component::MatVecInterface::MatVecInterface(void)
{
}

/**
 * Constructor
 */
gridpack::component::MatVecInterface::~MatVecInterface(void)
{
}

/**
 * Return size of matrix block on the diagonal contributed by component
 * @param isize, jsize: number of rows and columns of matrix
 *        block
 * @return: false if network component does not contribute
 *        matrix element
 */
bool gridpack::component::MatVecInterface::matrixDiagSize(int *isize, int *jsize) const
{
  *isize = 0;
  *jsize = 0;
  return false;
}

/**
 * Return the values of for a diagonal matrix block. The values are
 * returned in row-major order.
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute
 *        matrix element
 */
bool gridpack::component::MatVecInterface::matrixDiagValues(void *values)
{
  return false;
}

/**
 * Return size of off-diagonal matrix block contributed by component. The values
 * are for the forward direction
 * @param isize, jsize: number of rows and columns of matrix
 *        block
 * @return: false if network component does not contribute
 *        matrix element
 */
bool gridpack::component::MatVecInterface::matrixForwardSize(int *isize, int *jsize) const
{
  *isize = 0;
  *jsize = 0;
  return false;
}

/**
 * Return the values of for an off-diagonl matrix block. The values are for the
 * forward direction and are returned in row-major order.
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute
 *        matrix element
 */
bool gridpack::component::MatVecInterface::matrixForwardValues(void *values)
{
  return false;
}

/**
 * Return size of off-diagonal matrix block contributed by component. The values
 * are for the reverse direction
 * @param isize, jsize: number of rows and columns of matrix
 *        block
 * @return: false if network component does not contribute
 *        matrix element
 */
bool gridpack::component::MatVecInterface::matrixReverseSize(int *isize, int *jsize) const
{
  *isize = 0;
  *jsize = 0;
  return false;
}

/**
 * Return the values of for an off-diagonl matrix block. The values are for the
 * reverse direction and are returned in row-major order.
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute
 *        matrix element
 */
bool gridpack::component::MatVecInterface::matrixReverseValues(void *values)
{
  return false;
}

/**
 * Return size of vector block contributed by component
 * @param isize: number of vector elements
 * @return: false if network component does not contribute
 *        vector element
 */
bool gridpack::component::MatVecInterface::vectorSize(int *isize) const
{
  *isize = 0;
  return false;
}

/**
 * Return the values of the vector block.
 * @param values: pointer to vector values
 * @return: false if network component does not contribute
 *        vector element
 */
bool gridpack::component::MatVecInterface::vectorValues(void *values)
{
  return false;
}

// The base implementation for bus and branch components.

/**
 * Simple constructor
 */
gridpack::component::BaseComponent::BaseComponent(void)
{
}

/**
 * Destructor
 */
gridpack::component::BaseComponent::~BaseComponent(void)
{
}

/**
 * Return the size of the component for use in packing and
 * unpacking routines. This might not be needed, but throw
 * it in for now.
 * @return: size of network component
 */
int gridpack::component::BaseComponent::size(void) const
{
  return sizeof(*this);
}

// Base implementation for a bus object. Provides a mechanism for the bus to
// provide a list of the branches that are directly connected to it as well as a
// mechanism for returning a list of the buses that are connected to it via a
// single branch

/**
 * Simple constructor
 */
gridpack::component::BaseBusComponent::BaseBusComponent(void)
{
}

/**
 * Simple destructor
 */
gridpack::component::BaseBusComponent::~BaseBusComponent(void)
{
}

/**
 * Add a pointer to the list of branches that a bus is connected to
 * @param branch: pointer to a branch that is connected to bus
 */
void
gridpack::component::BaseBusComponent::addBranch(boost::shared_ptr<gridpack::component::BaseComponent> branch)
{
  p_branches.push_back(branch);
}

/**
 * Add a pointer to the list of buses that a bus is connected to via
 * a branch
 * @param bus: pointer to a branch that is connected to bus
 */
void
gridpack::component::BaseBusComponent::addBus(boost::shared_ptr<gridpack::component::BaseComponent> bus)
{
  p_buses.push_back(bus);
}

/**
 * Get pointers to branches that are connected to bus
 * @return: list of pointers to neighboring branches
 */
std::vector<boost::shared_ptr<gridpack::component::BaseComponent> >
gridpack::component::BaseBusComponent::getNeighborBranches(void) const
{
  return p_branches;
}

/**
 * Get pointers to buses that are connected to calling bus via a branch
 * @return: list of pointers to neighboring buses
 */
std::vector<boost::shared_ptr<gridpack::component::BaseComponent> >
gridpack::component::BaseBusComponent::getNeighborBuses(void) const
{
  return p_buses;
}

/**
 * Clear all pointers to neighboring branches
 */
void gridpack::component::BaseBusComponent::clearBranches(void)
{
  p_branches.clear();
}

/**
 * Clear all pointers to neighboring buses
 */
void gridpack::component::BaseBusComponent::clearBuses(void)
{
  p_buses.clear();
}

// Base implementation for a branch object. Provides a mechanism for the branch to
// provide the buses at either end of the branch

/**
 * Simple constructor
 */
gridpack::component::BaseBranchComponent::BaseBranchComponent(void)
{
}

/**
 * Simple destructor
 */
gridpack::component::BaseBranchComponent::~BaseBranchComponent(void)
{
}

/**
 * Set pointer to bus at one end of branch
 * @param: pointer to bus
 */
void
gridpack::component::BaseBranchComponent::setBus1(boost::shared_ptr<gridpack::component::BaseComponent> bus)
{
  p_bus2 = bus;
}

/**
 * Set pointer to bus at other end of branch
 * @param: pointer to bus
 */
void
gridpack::component::BaseBranchComponent::setBus2(boost::shared_ptr<gridpack::component::BaseComponent> bus)
{
  p_bus2 = bus;
}

/**
 * Get pointer to bus at one end of branch
 */
boost::shared_ptr<gridpack::component::BaseComponent>
  gridpack::component::BaseBranchComponent::getBus1(void) const
{
  return p_bus1;
}

/**
 * Get pointer to bus at other end of branch
 */
boost::shared_ptr<gridpack::component::BaseComponent>
  gridpack::component::BaseBranchComponent::getBus2(void) const
{
  return p_bus2;
}

/**
 * Clear bus pointers
 */
void gridpack::component::BaseBranchComponent::clearBuses(void)
{
  p_bus1.reset();
  p_bus2.reset();
}
