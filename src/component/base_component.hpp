// -------------------------------------------------------------
/**
 * @file   base_component.hpp
 * @author Bruce Palmer
 * @date   2013-07-11 12:29:46 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _base_component_h_
#define _base_component_h_

#include <vector>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "boost/smart_ptr/weak_ptr.hpp"
#include "gridpack/component/data_collection.hpp"

#include <boost/serialization/export.hpp>

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
     * Return size of matrix block on the diagonal contributed by component and
     * the global index of this component
     * @param idx: global index of this component
     * @param isize, jsize: number of rows and columns of matrix
     *        block
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixDiagSize(int *idx, int *isize, int *jsize) const;

    /**
     * Return the values of for a diagonal matrix block. The values are
     * returned in row-major order. Also return the global index of component
     * @param idx: global index of this component
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixDiagValues(int *idx, void *values);

    /**
     * Return size of off-diagonal matrix block contributed by component. The
     * values are for the forward direction. Also return matrix location using
     * global indices
     * @param idx, jdx: matrix location using global indices
     * @param isize, jsize: number of rows and columns of matrix
     *        block
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixForwardSize(int *idx, int *jdx, int *isize, int *jsize) const;

    /**
     * Return the values of for an off-diagonl matrix block. The values are
     * for the forward direction and are returned in row-major order. Also
     * return matrix location using global indices
     * @param idx, jdx: matrix location using global indices
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixForwardValues(int *idx, int *jdx, void *values);

    /**
     * Return size of off-diagonal matrix block contributed by component. The
     * values are for the reverse direction. Also return matrix location using
     * global indices
     * @param idx, jdx: matrix location using global indices
     * @param isize, jsize: number of rows and columns of matrix
     *        block
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixReverseSize(int *idx, int *jdx, int *isize, int *jsize) const;

    /**
     * Return the values of for an off-diagonl matrix block. The values are
     * for the reverse direction and are returned in row-major order. Also
     * return matrix location using global indices
     * @param idx, jdx: matrix location using global indices
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixReverseValues(int *idx, int *jdx, void *values);

    /**
     * Return size of vector block contributed by component and location using
     * global indices
     * @param idx: vector location using global indices
     * @param isize: number of vector elements
     * @return: false if network component does not contribute
     *        vector element
     */
    virtual bool vectorSize(int *idx, int *isize) const;

    /**
     * Return the values of the vector block and location using global indices
     * @param idx: vector location using global indices
     * @param values: pointer to vector values
     * @return: false if network component does not contribute
     *        vector element
     */
    virtual bool vectorValues(int *idx, void *values);

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

    /**
     * Load data from DataCollection object into corresponding
     * component. This needs to be implemented by every component
     * @param data: data collection associated with component
     */
    virtual void load(
        const boost::shared_ptr<gridpack::component::DataCollection> &data);

    /**
     * Return the size of the buffer needed for data exchanges. Note that this
     * must be the same size for all bus and all branch objects (branch buffers
     * do not need to be the same size as bus buffers), even if all objects
     * do not require the same parameters. Thus, the buffer must be big enough
     * to exchange all variables that an object might need, even if individual
     * objects don't need all the variables
     */
    virtual int getXCBufSize(void);

    /**
     * Assign the location of the data exchange buffer. These buffers are
     * allocated and deallocated by the network
     * @param buf: void pointer to exchange buffer
     */
    virtual void setXCBuf(void *buf);

    /**
     * Set an internal variable that can be used to control the behavior of the
     * component. This function doesn't need to be implemented, but if needed,
     * it can be used to change the behavior of the network in different phases
     * of the calculation. For example, if a different matrix needs to be
     * generated at different times, the mode of the calculation can changed to
     * get different values from the MatVecInterface functions
     * @param mode: integer indicating which mode should be used
     */
    virtual void setMode(int mode);

  protected:
    /**
     * A buffer that can be used for exchanging component data. This is
     * allocated by the network based on an inquiry to the getXCBusSize method
     */
     void *p_XCBuf;

    /**
     * Size (in bytes) of buffer p_XCBuf
     */
     int p_XCBufSize;

     /**
      * Current mode
      */
     int p_mode;

  private:

  friend class boost::serialization::access;

  /// Serialization "pack" routine
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const
  {
    // p_XCBuf and p_XCBufSize are managed somewhere else; they will
    // have to be initialized
    ar << p_XCBufSize
       << p_mode;
  }

  /// Serialization "unpack" routine
  template<class Archive>
  void load(Archive & ar, const unsigned int version)
  {
    ar >> p_XCBufSize
       >> p_mode;
  }

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    boost::serialization::split_member(ar, *this, version);
  }

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
    void addBranch(const boost::shared_ptr<BaseComponent> & branch);

    /**
     * Add a pointer to the list of buses that a bus is connected to via
     * a branch
     * @param bus: pointer to a branch that is connected to bus
     */
    void addBus(const boost::shared_ptr<BaseComponent> & bus);

    /**
     * Get pointers to branches that are connected to bus
     * @param nghbrs: list of pointers to neighboring branches
     */
    void getNeighborBranches(std::vector<boost::shared_ptr<BaseComponent> > &nghbrs) const;

    /**
     * Get pointers to buses that are connected to calling bus via a branch
     * @param nghbrs: list of pointers to neighboring buses
     */
    void getNeighborBuses(std::vector<boost::shared_ptr<BaseComponent> > &nghbrs) const;

    /**
     * Clear all pointers to neighboring branches
     */
    void clearBranches(void);

    /**
     * Clear all pointers to neighboring buses
     */
    void clearBuses(void);

    /**
     * Set reference bus status
     * @param status: reference bus status
     */
    void setReferenceBus(bool status);

    /**
     * Get reference bus status
     * @return: reference bus status
     */
    bool getReferenceBus(void) const;

    /**
     * Set global index of bus
     * @param idx: global index of bus
     */
    void setGlobalIndex(int idx);

    /**
     * Get global index of bus
     * @param idx: global index of bus
     */
    void getGlobalIndex(int *idx) const;

  private:
    /**
     * Branches that are connect to bus
     */
    std::vector<boost::weak_ptr<BaseComponent> > p_branches;
    /**
     * Buses that are connect to bus via a branch
     */
    std::vector<boost::weak_ptr<BaseComponent> > p_buses;

    /**
     * Is this a reference bus?
     */
    bool p_refBus;

    /**
     * Global index of bus
     */
    int p_idx;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
     ar & boost::serialization::base_object<BaseComponent>(*this)
       & p_refBus
       & p_idx;

     // p_branches and p_buses are ignored, but that may change
  }
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
    void setBus1(const boost::shared_ptr<BaseComponent> & bus);

    /**
     * Set pointer to bus at other end of branch
     * @param bus: pointer to bus
     */
    void setBus2(const boost::shared_ptr<BaseComponent> & bus);

    /**
     * Get pointer to bus at one end of branch
     * @return: pointer to bus 1
     */
    boost::shared_ptr<BaseComponent> getBus1(void) const;

    /**
     * Get pointer to bus at other end of branch
     * @return: pointer to bus 2
     */
    boost::shared_ptr<BaseComponent> getBus2(void) const;

    /**
     * Clear bus pointers
     */
    void clearBuses(void);

    /**
     * Set the global index of the branch and set the index values for the two
     * buses at either end of the branch.
     * @param branch_idx: global index of branch
     * @param bus1_idx: global index of "from" bus
     * @param bus2_idx: global index of "to" bus
     */
    void setGlobalIndices(int branch_idx, int bus1_idx, int bus2_idx);

    /**
     * Get the global index of the branch and get the index values for the two
     * buses at either end of the branch.
     * @param branch_idx: global index of branch
     * @param bus1_idx: global index of "from" bus
     * @param bus2_idx: global index of "to" bus
     */
    void getGlobalIndices(int *branch_idx, int *bus1_idx, int *bus2_idx) const;

  private:
    /**
     *  Pointers to buses at either end of branch
     */
    boost::weak_ptr<BaseComponent> p_bus1;
    boost::weak_ptr<BaseComponent> p_bus2;

    /**
     *  Global index of branch and global indices of buses at each end of branch
     */
    int p_branch_idx;
    int p_bus1_idx, p_bus2_idx;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
     ar & boost::serialization::base_object<BaseComponent>(*this)
       & p_branch_idx
       & p_bus1_idx
       & p_bus2_idx;
     // p_bus1 and p_bus2 are ignored, but that may change
  }

};

}    // component
}    // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseComponent);
BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseBusComponent);
BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseBranchComponent);


#endif
