// -------------------------------------------------------------
/**
 * @file   base_component.hpp
 * @author Bruce Palmer
 * @date   April 8, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _base_component_h_
#define _base_component_h_

#include "boost/smart_ptr/shared_ptr.hpp"

// TODO: Might want to put MatrixIndices and VectorIndex operations into a
//       separate class since these can probably be implemented once for all
//       network components
namespace gridpack{
namespace component{

class MatVecInterface {
  public:
    /**
     * Constructor
     */
    MatVecInterface(void);

    /**
     * Destructor
     */
    virtual ~MatVecInterface(void);

    /**
     * Return size of matrix block contributed by component
     * @param isize, jsize: number of rows and columns of matrix
     *        block
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixSize(int *isize, int *jsize) const;

    /**
     * Return the values of the matrix block. The values are
     * returned in row-major order.
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixValues(void *values);

    /**
     * Return size of vector block contributed by component
     * @param isize: number of vector elements
     * @return: false if network component does not contribute
     *        vector element
     */
    virtual bool vectorSize(int *isize) const;

    /**
     * Return the values of the vector block.
     * @param values: pointer to vector values
     * @return: false if network component does not contribute
     *        vector element
     */
    virtual bool vectorValues(void *values);

    //TODO: May need to include routines that support moving values from vectors
    //      back into network components.

  private:

};

// -------------------------------------------------------------
//  class BaseField:
//  This class implements some basic functions that can be
//  expected from any field on the network.
// -------------------------------------------------------------
class BaseComponent
  : public MatVecInterface {
  public:
    /**
     * Simple constructor
     */
    BaseComponent(void);

    /**
     * Destructor
     */
    virtual ~BaseComponent(void);

    /**
     * Return the size of the component for use in packing and
     * unpacking routines. This might not be needed, but throw
     * it in for now.
     * @return: size of network component
     */
    virtual int size(void) const;

  private:

};

class BaseBusComponent
  : public BaseComponent {
  public:
    /**
     * Simple constructor
     */
    BaseBusComponent(void);

    /**
     * Simple destructor
     */
    virtual ~BaseBusComponent(void);

    /**
     * Add a pointer to the list of branches that a bus is connected to
     * @param branch: pointer to a branch that is connected to bus
     */
    void addBranch(boost::shared_ptr<BaseComponent> branch);

    /**
     * Add a pointer to the list of buses that a bus is connected to via
     * a branch
     * @param bus: pointer to a branch that is connected to bus
     */
    void addBus(boost::shared_ptr<BaseComponent> bus);

    /**
     * Get pointers to branches that are connected to bus
     * @return: list of pointers to neighboring branches
     */
    std::vector<boost::shared_ptr<BaseComponent> > getNeighborBranches(void) const;

    /**
     * Get pointers to buses that are connected to calling bus via a branch
     * @return: list of pointers to neighboring buses
     */
    std::vector<boost::shared_ptr<BaseComponent> > getNeighborBuses(void) const;

    /**
     * Clear all pointers to neighboring branches
     */
    void clearBranches();

    /**
     * Clear all pointers to neighboring buses
     */
    void clearBuses();

  private:
    /**
     * Branches that are connect to bus
     */
    std::vector<boost::shared_ptr<BaseComponent> > p_branches;
    /**
     * Buses that are connect to bus via a branch
     */
    std::vector<boost::shared_ptr<BaseComponent> > p_buses;
};

class BaseBranchComponent
  : public BaseComponent {
  public:
    /**
     * Simple constructor
     */
    BaseBranchComponent(void);

    /**
     * Simple destructor
     */
    virtual ~BaseBranchComponent(void);

    /**
     * Set pointer to bus at one end of branch
     * @param bus: pointer to bus
     */
    void setBus1(boost::shared_ptr<BaseComponent> bus);

    /**
     * Set pointer to bus at other end of branch
     * @param bus: pointer to bus
     */
    void setBus2(boost::shared_ptr<BaseComponent> bus);

    /**
     * Get pointer to bus at one end of branch
     */
    boost::shared_ptr<BaseComponent> getBus1(void) const;

    /**
     * Get pointer to bus at other end of branch
     */
    boost::shared_ptr<BaseComponent> getBus2(void) const;

    /**
     * Clear bus pointers
     */
    void clearBuses(void);

  private:
    /**
     *  Pointers to buses at either end of branch
     */
    boost::shared_ptr<BaseComponent> p_bus1;
    boost::shared_ptr<BaseComponent> p_bus2;

};

}    // component
}    // gridpack
#endif
