// -------------------------------------------------------------
/**
 * @file   base_component.hpp
 * @author Bruce Palmer
 * @date   2013-07-17 10:05:04 d3g096
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
     * Return size of matrix block on the diagonal contributed by component
     * @param isize, jsize: number of rows and columns of matrix
     *        block
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixDiagSize(int *isize, int *jsize) const;

    /**
     * Return the values of for a diagonal matrix block. The values are
     * returned in row-major order.
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixDiagValues(void *values);

    /**
     * Return size of off-diagonal matrix block contributed by component. The
     * values are for the forward direction.
     * @param isize, jsize: number of rows and columns of matrix
     *        block
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixForwardSize(int *isize, int *jsize) const;

    /**
     * Return the values of for an off-diagonl matrix block. The values are
     * for the forward direction and are returned in row-major order.
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixForwardValues(void *values);

    /**
     * Return size of off-diagonal matrix block contributed by component. The
     * values are for the reverse direction.
     * @param isize, jsize: number of rows and columns of matrix
     *        block
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixReverseSize(int *isize, int *jsize) const;

    /**
     * Return the values of for an off-diagonl matrix block. The values are
     * for the reverse direction and are returned in row-major order.
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixReverseValues(void *values);

    /**
     * Return size of vector block contributed by component
     * @param isize: number of vector elements
     * @return: false if network component does not contribute
     *        vector element
     */
    virtual bool vectorSize(int *isize) const;

    /**
     * Return the values of the vector block
     * @param values: pointer to vector values
     * @return: false if network component does not contribute
     *        vector element
     */
    virtual bool vectorValues(void *values);

    /**
     * Set values in the bus or branch component based on values in a vector or
     * matrix
     * @param values: values in vector or matrix
     */
    virtual void setValues(void *values);

    /**
     * Set the matrix index for diagonal matrix components or vector component,
     * based on location of component in network
     * @param idx: value of index
     */
    void setMatVecIndex(int idx);

    /**
     * Get the matrix index for diagonal matrix components or vector component,
     * based on location of component in network
     * @return: value of index
     */
    void getMatVecIndex(int *idx) const;

    /**
     * Set the matrix indices for matrix components, based on location of component
     * in network
     * @param idx, jdx: value of indices
     */
    void setMatVecIndices(int idx, int jdx);

    /**
     * Get the matrix indices for matrix components, * based on location of component
     * in network
     * @param idx, jdx: value of indices
     */
    void getMatVecIndices(int *idx, int *jdx) const;

    //TODO: May need to include routines that support moving values from vectors
    //      back into network components.

  private:
    int p_ival, p_idx, p_jdx;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & p_ival
       & p_idx
       & p_jdx;
  }


};

// -------------------------------------------------------------
//  class BaseComponent:
//  This class implements some basic functions that can be
//  expected from any component on the network.
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

    /**
     * Copy a string for output into buffer. The behavior of this method can be
     * altered by inputting different values for the signal string
     * @param string: buffer containing string to be written to output
     * @param signal: string to control behavior of routine (e.g. what
     * properties to write
     * @return: true if component is writing a contribution, false otherwise
     */
    virtual bool serialWrite(char *string, char *signal = NULL);

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
    ar << boost::serialization::base_object<MatVecInterface>(*this)
       << p_XCBufSize
       << p_mode;
  }

  /// Serialization "unpack" routine
  template<class Archive>
  void load(Archive & ar, const unsigned int version)
  {
    ar >> boost::serialization::base_object<MatVecInterface>(*this)
       >> p_XCBufSize
       >> p_mode;
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER();

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
     * Set original index (from input file)
     * @param idx: original index from network
     */
    void setOriginalIndex(int idx);

    /**
     * Get original index
     * @return: original index from network
     */
    int getOriginalIndex(void) const;

    /**
     * Set global index
     * @param idx: global index from network
     */
    void setGlobalIndex(int idx);

    /**
     * Get global index
     * @return: global index from network
     */
    int getGlobalIndex(void) const;

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
     * Original and global indices, derived from network class
     */
    int p_originalIndex;
    int p_globalIndex;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
     ar & boost::serialization::base_object<BaseComponent>(*this)
       & p_refBus;

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
     * Set original index for bus 1
     * @param idx: original index for bus 1 (assigned from input file)
     */
    void setBus1OriginalIndex(int idx);

    /**
     * Set original index for bus 2
     * @param idx: original index for bus 2 (assigned from input file)
     */
    void setBus2OriginalIndex(int idx);

    /**
     * Set global index (from network) for bus 1
     * @param idx: global index for bus 1
     */
    void setBus1GlobalIndex(int idx);

    /**
     * Set global index (from network) for bus 2
     * @param idx: global index for bus 2
     */
    void setBus2GlobalIndex(int idx);

    /**
     * Get original index for bus 1
     * @return: original index for bus 1
     */
    int getBus1OriginalIndex(void) const;

    /**
     * Get original index for bus 2
     * @return: original index for bus 2
     */
    int getBus2OriginalIndex(void) const;

    /**
     * Get global index for bus 1
     * @return: global index for bus 1
     */
    int getBus1GlobalIndex(void) const;

    /**
     * Get global index for bus 2
     * @return: global index for bus 2
     */
    int getBus2GlobalIndex(void) const;

  private:
    /**
     *  Pointers to buses at either end of branch
     */
    boost::weak_ptr<BaseComponent> p_bus1;
    boost::weak_ptr<BaseComponent> p_bus2;

    /**
     *  Original indices for bus 1 and bus 2 (assigned from input file)
     */
    int p_originalBus1Index;
    int p_originalBus2Index;

    /**
     *  Global indices for bus 1 and bus 2
     */
    int p_globalBus1Index;
    int p_globalBus2Index;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
     ar & boost::serialization::base_object<BaseComponent>(*this);
  }

};

}    // component
}    // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::component::MatVecInterface);
BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseComponent);
BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseBusComponent);
BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseBranchComponent);


#endif
