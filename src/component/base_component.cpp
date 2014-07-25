/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_component.cpp
 * @author Bruce Palmer
 * @date   2013-07-11 12:25:44 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/component/base_component.hpp"

// Base implementation of the MatVecInterface. These functions should be
// overwritten in actual components

namespace gridpack {
namespace component {
/**
 * Constructor
 */
MatVecInterface::MatVecInterface(void)
{
  p_ival = 0;
  p_idx = 0;
  p_jdx = 0;
}

/**
 * Constructor
 */
MatVecInterface::~MatVecInterface(void)
{
}

/**
 * Return size of matrix block on the diagonal contributed by component
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool MatVecInterface::matrixDiagSize( int *isize,
             int *jsize) const
{
  *isize = 0;
  *jsize = 0;
  return false;
}

/**
 * Return the values of for a diagonal matrix block. The values are
 * returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool MatVecInterface::matrixDiagValues(ComplexType *values)
{
  return false;
}

/**
 * Return size of off-diagonal matrix block contributed by component. The
 * values are for the forward direction.
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool MatVecInterface::matrixForwardSize(int *isize,
       int *jsize) const
{
  *isize = 0;
  *jsize = 0;
  return false;
}

/**
 * Return the values of for an off-diagonl matrix block. The values are
 * for the forward direction and are returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool MatVecInterface::matrixForwardValues(ComplexType *values)
{
  return false;
}

/**
 * Return size of off-diagonal matrix block contributed by component. The
 * values are for the reverse direction.
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool MatVecInterface::matrixReverseSize(int *isize,
       int *jsize) const
{
  *isize = 0;
  *jsize = 0;
  return false;
}

/**
 * Return the values of for an off-diagonl matrix block. The values are
 * for the reverse direction and are returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool MatVecInterface::matrixReverseValues(ComplexType *values)
{
  return false;
}

/**
 * Return size of vector block contributed by component
 * @param isize number of vector elements
 * @return false if network component does not contribute
 *        vector element
 */
bool MatVecInterface::vectorSize(int *isize) const
{
  *isize = 0;
  return false;
}

/**
 * Return the values of the vector block
 * @param values pointer to vector values
 * @return false if network component does not contribute
 *        vector element
 */
bool MatVecInterface::vectorValues(ComplexType *values)
{
  return false;
}

/**
 * Set values in the bus or branch component based on values in a vector or
 * matrix
 * @param values values in vector or matrix
 */
void MatVecInterface::setValues(ComplexType *values)
{
}

/**
 * Set the matrix index for diagonal matrix components or vector component,
 * based on location of component in network
 * @param idx value of index
 */
void MatVecInterface::setMatVecIndex(int idx)
{
  p_ival = idx;
}

/**
 * Get the matrix index for diagonal matrix components or vector component,
 * based on location of component in network
 * @return value of index
 */
void MatVecInterface::getMatVecIndex(int *idx) const
{
  *idx = p_ival;
}

/**
 * Set the matrix indices for matrix components, based on location of
 * component
 * in network
 * @param idx,jdx value of indices
 */
void MatVecInterface::setMatVecIndices(int idx, int jdx)
{
  p_idx = idx;
  p_jdx = jdx;
}

/**
 * Get the matrix indices for matrix components, * based on location of
 * component
 * in network
 * @param idx,jdx value of indices
 */
void MatVecInterface::getMatVecIndices(int *idx, int *jdx) const
{
  *idx = p_idx;
  *jdx = p_jdx;
}

// base implementation for the generalized matrix-vector interface

/**
 * Constructor
 */
GenMatVecInterface::GenMatVecInterface(void)
{
}

/**
 * Destructor
 */
GenMatVecInterface::~GenMatVecInterface(void)
{
}

/**
 * Return number of rows in matrix from component
 * @return number of rows from component
 */
int GenMatVecInterface::matrixNumRows() const
{
  return 0;
}

/**
 * Return number of columns in matrix from component
 * @return number of columnsows from component
 */
int GenMatVecInterface::matrixNumCols() const
{
  return 0;
}

/**
 * Set row indices corresponding to the rows contributed by this
 * component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @param idx matrix index of row irow
 */
void GenMatVecInterface::matrixSetRowIndex(int irow, int idx)
{
}

/**
 * Set column indices corresponding to the columns contributed by this
 * component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @param idx matrix index of column icol
 */
void GenMatVecInterface::matrixSetColIndex(int icol, int idx)
{
}

/**
 * Get the row index corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @return matrix index of row irow
 */
int GenMatVecInterface::matrixGetRowIndex(int idx)
{
  return -1;
}

/**
 * Get the column index corresponding to the columns contributed by this component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @return matrix index of column icol
 */
int GenMatVecInterface::matrixGetColIndex(int idx)
{
  return -1;
}

/**
 * Return the number of matrix values contributed by this component
 * @return number of matrix values
 */
int GenMatVecInterface::matrixNumValues() const
{
  return 0;
}

/**
 * Get a list of matrix values contributed by this component and their
 * matrix indices
 * @param values list of matrix element values
 * @param rows row indices for the matrix elements
 * @param cols column indices for the matrix elements
 */
void GenMatVecInterface::matrixGetValues(ComplexType *values, int *rows, int*cols)
{
}

/**
 * Return number of elements in vector from component
 * @return number of elements contributed from component
 */
int GenMatVecInterface::vectorNumElements() const
{
  return 0;
}

/**
 * Set indices corresponding to the elements contributed by this
 * component
 * @param ielem index of element contributed by this component
 * (e.g. if component contributes 3 elements then ielem is between
 * 0 and 2)
 * @param idx vector index of element ielem
 */
void GenMatVecInterface::vectorSetElementIndex(int ielem, int idx)
{
}

/**
 * Get list of element indices from component
 * @param idx list of indices that component maps onto
 */
void GenMatVecInterface::vectorGetElementIndices(int *idx)
{
}

/**
 * Get a list of vector values contributed by this component and their
 * indices
 * @param values list of vector element values
 * @param idx indices for the vector elements
 */
void GenMatVecInterface::vectorGetElementValues(ComplexType *values, int *idx)
{
}

/**
 * Transfer vector values to component
 * @param values list of vector element values
 */
void GenMatVecInterface::vectorSetElementValues(ComplexType *values)
{
}
 
// The base implementation for bus and branch components.

/**
 * Simple constructor
 */
BaseComponent::BaseComponent(void)
  : p_XCBuf(NULL), p_XCBufSize(0), p_mode(0)
{
}

/**
 * Destructor
 */
BaseComponent::~BaseComponent(void)
{
}

/**
 * Load data from DataCollection object into corresponding
 * component. This needs to be implemented by every component
 * @param data data collection associated with component
 */
void BaseComponent::load(
  const boost::shared_ptr<DataCollection> &data)
{
  // This implementation is a no-op and is included in BaseComponent so that
  // a generic load method can be defined in the base factory class.
}

/**
 * Return the size of the buffer needed for data exchanges. Note that this
 * must be the same size for all bus and all branch objects (branch buffers
 * do not need to be the same size as bus buffers), even if all objects
 * do not require the same parameters. Thus, the buffer must be big enough
 * to exchange all variables that an object might need, even if individual
 * objects don't need all the variables
 */
int BaseComponent::getXCBufSize(void)
{
  return p_XCBufSize;
}
/**
 * Assign the location of the data exchange buffer. These buffers are
 * allocated and deallocated by the network
 * @param buf void pointer to exchange buffer
 */
void BaseComponent::setXCBuf(void *buf)
{
  // FIXME: ?
  if (buf == NULL) {
    BaseComponent::p_XCBuf = NULL;
  } else {
    BaseComponent::p_XCBuf = buf;
  }
}

/**
 * Set an internal variable that can be used to control the behavior of the
 * component. This function doesn't need to be implemented, but if needed,
 * it can be used to change the behavior of the network in different phases
 * of the calculation. For example, if a different matrix needs to be
 * generated at different times, the mode of the calculation can changed to
 * get different values from the MatVecInterface functions
 * @param mode integer indicating which mode should be used
 */
void BaseComponent::setMode(int mode)
{
  p_mode = mode;
}

/**
 * Copy a string for output into buffer. The behavior of this method can be
 * altered by inputting different values for the signal string
 * @param string buffer containing string to be written to output
 * @param bufsize size of string buffer in bytes
 * @param signal string to control behavior of routine (e.g. what
 * properties to write
 * @return true if component is writing a contribution, false otherwise
 */
bool BaseComponent::serialWrite(char *string, const int bufsize, const char *signal)
{
  return false;
  // This is defined so that generic operations for writing strings from buses
  // and branches can be built
}

// Base implementation for a bus object. Provides a mechanism for the bus to
// provide a list of the branches that are directly connected to it as well as a
// mechanism for returning a list of the buses that are connected to it via a
// single branch

/**
 * Simple constructor
 */
BaseBusComponent::BaseBusComponent(void)
  : p_refBus(false)
{
  
}

/**
 * Simple destructor
 */
BaseBusComponent::~BaseBusComponent(void)
{
}

/**
 * Add a pointer to the list of branches that a bus is connected to
 * @param branch pointer to a branch that is connected to bus
 */
void
BaseBusComponent::addBranch(const
  boost::shared_ptr<BaseComponent> & branch)
{
  boost::weak_ptr<BaseComponent> tbranch(branch);
  p_branches.push_back(tbranch);
}

/**
 * Add a pointer to the list of buses that a bus is connected to via
 * a branch
 * @param bus pointer to a branch that is connected to bus
 */
void
BaseBusComponent::addBus(const
  boost::shared_ptr<BaseComponent> & bus)
{
  boost::weak_ptr<BaseComponent> tbus(bus);
  p_buses.push_back(tbus);
}

/**
 * Get pointers to branches that are connected to bus
 * @param nghbrs list of pointers to neighboring branches
 */
void BaseBusComponent::getNeighborBranches(
  std::vector<boost::shared_ptr<BaseComponent> > &nghbrs) const
{
  nghbrs.clear();
  int i;
  int size = p_branches.size();
  for (i=0; i<size; i++) {
    boost::shared_ptr<BaseComponent> branch = p_branches[i].lock();
    nghbrs.push_back(branch);
  }
}

/**
 * Get pointers to buses that are connected to calling bus via a branch
 * @param nghbrs list of pointers to neighboring buses
 */
void BaseBusComponent::getNeighborBuses(
  std::vector<boost::shared_ptr<BaseComponent> > &nghbrs) const
{
  nghbrs.clear();
  int i;
  int size = p_buses.size();
  for (i=0; i<size; i++) {
    boost::shared_ptr<BaseComponent> bus = p_buses[i].lock();
    nghbrs.push_back(bus);
  }
}

/**
 * Clear all pointers to neighboring branches
 */
void BaseBusComponent::clearBranches(void)
{
  p_branches.clear();
}

/**
 * Clear all pointers to neighboring buses
 */
void BaseBusComponent::clearBuses(void)
{
  p_buses.clear();
}

/**
 * Set reference bus status
 * @param status reference bus status
 */
void BaseBusComponent::setReferenceBus(bool status)
{
  p_refBus = status;
}

/**
 * Get reference bus status
 * @return reference bus status
 */
bool BaseBusComponent::getReferenceBus(void) const
{
  return p_refBus;
}

/**
 * Set original index (from input file)
 * @param idx original index from network
 */
void BaseBusComponent::setOriginalIndex(int idx)
{
  p_originalIndex = idx;
}

/**
 * Get original index
 * @return original index from network
 */
int BaseBusComponent::getOriginalIndex(void) const
{
  return p_originalIndex;
}

/**
 * Set global index
 * @param idx global index from network
 */
void BaseBusComponent::setGlobalIndex(int idx)
{
  p_globalIndex = idx;
}

/**
 * Get global index
 * @return global index from network
 */
int BaseBusComponent::getGlobalIndex(void) const
{
  return p_globalIndex;
}

// Base implementation for a branch object. Provides a mechanism for the branch to
// provide the buses at either end of the branch

/**
 * Simple constructor
 */
BaseBranchComponent::BaseBranchComponent(void)
{
}

/**
 * Simple destructor
 */
BaseBranchComponent::~BaseBranchComponent(void)
{
}

/**
 * Set pointer to bus at one end of branch
 * @param bus pointer to bus
 */
void BaseBranchComponent::setBus1(const
  boost::shared_ptr<BaseComponent> & bus)
{
  p_bus1 = boost::weak_ptr<BaseComponent>(bus);
}

/**
 * Set pointer to bus at other end of branch
 * @param bus pointer to bus
 */
void BaseBranchComponent::setBus2(const
  boost::shared_ptr<BaseComponent> & bus)
{
  p_bus2 = boost::weak_ptr<BaseComponent>(bus);
}

/**
 * Get pointer to bus at one end of branch
 * @return pointer to bus 1
 */
boost::shared_ptr<BaseComponent>
  BaseBranchComponent::getBus1(void) const
{
  boost::shared_ptr<BaseComponent> ret(p_bus1);
  return ret;
}

/**
 * Get pointer to bus at other end of branch
 * @return pointer to bus 2
 */
boost::shared_ptr<BaseComponent>
  BaseBranchComponent::getBus2(void) const
{
  boost::shared_ptr<BaseComponent> ret(p_bus2);
  return ret;
}

/**
 * Clear bus pointers
 */
void BaseBranchComponent::clearBuses(void)
{
  p_bus1.reset();
  p_bus2.reset();
}

/**
 * Set original index for bus 1
 * @param idx original index for bus 1 (assigned from input * file)
 */
void BaseBranchComponent::setBus1OriginalIndex(int idx)
{
  p_originalBus1Index = idx;
}

/**
 * Set original index for bus 2
 * @param idx original index for bus 2 (assigned from input * file)
 */
void BaseBranchComponent::setBus2OriginalIndex(int idx)
{
  p_originalBus2Index = idx;
}

/**
 * Set global index (from network) for bus 1
 * @param idx global index for bus 1
 */
void BaseBranchComponent::setBus1GlobalIndex(int idx)
{
  p_globalBus1Index = idx;
}

/**
 * Set global index (from network) for bus 2
 * @param idx global index for bus 2
 */
void BaseBranchComponent::setBus2GlobalIndex(int idx)
{
  p_globalBus2Index = idx;
}

/**
 * Get original index for bus 1
 * @return original index for bus 1
 */
int BaseBranchComponent::getBus1OriginalIndex(void) const
{
  return p_originalBus1Index;
}

/**
 * Get original index for bus 2
 * @return original index for bus 2
 */
int BaseBranchComponent::getBus2OriginalIndex(void) const
{
  return p_originalBus2Index;
}

/**
 * Get global index for bus 1
 * @return global index for bus 1
 */
int BaseBranchComponent::getBus1GlobalIndex(void) const
{
  return p_globalBus1Index;
}

/**
 * Get global index for bus 2
 * @return global index for bus 2
 */
int BaseBranchComponent::getBus2GlobalIndex(void) const
{
  return p_globalBus2Index;
}

}  // component
}  // gridpack
