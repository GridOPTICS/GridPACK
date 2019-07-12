/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_component.hpp
 * @author Bruce Palmer
 * @date   2019-07-12 11:16:04 d3g096
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
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/component/base_matrix_ifc.hpp"

#include <boost/serialization/export.hpp>

// TODO: Might want to put MatrixIndices and VectorIndex operations into a
//       separate class since these can probably be implemented once for all
//       network components
namespace gridpack{
namespace component{
// -------------------------------------------------------------
//  class BaseComponent:
//  This class implements some basic functions that can be
//  expected from any component on the network.
// -------------------------------------------------------------
class BaseComponent : public BaseMatrixInterface, public GenMatVecInterface {
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
     * @param data data collection associated with component
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
     * @return size of buffer
     */
    virtual int getXCBufSize(void);

    /**
     * Assign the location of the data exchange buffer. These buffers are
     * allocated and deallocated by the network
     * @param buf void pointer to exchange buffer
     */
    virtual void setXCBuf(void *buf);

    /**
     * Return the location of the data exchange buffer
     * @param buf void pointer to exchange buffer
     */
    virtual void getXCBuf(void **bus);

    /**
     * Set an internal variable that can be used to control the behavior of the
     * component. This function doesn't need to be implemented, but if needed,
     * it can be used to change the behavior of the network in different phases
     * of the calculation. For example, if a different matrix needs to be
     * generated at different times, the mode of the calculation can changed to
     * get different values from the MatVecInterface functions
     * @param mode integer indicating which mode should be used
     */
    virtual void setMode(int mode);

    /**
     * Copy a string for output into buffer. The behavior of this method can be
     * altered by inputting different values for the signal string
     * @param string buffer containing string to be written to output
     * @param bufsize size of string buffer in bytes
     * @param signal string to control behavior of routine (e.g. what
     * properties to write)
     * @return true if component is writing a contribution, false otherwise
     */
    virtual bool serialWrite(char *string, const int bufsize, const char *signal = NULL);

    /**
     * Save state variables inside the component to a DataCollection object.
     * This can be used as a way of moving data in a way that is useful for
     * creating output or for copying state data from one network to another.
     * @param data data collection object into which new values are inserted
     */
    virtual void saveData(boost::shared_ptr<gridpack::component::DataCollection>
        data);

    /**
     * Retrieve an opaque data item from component. Different items may be
     * returned based on the value of signal.
     * @param data item to retrieve from component
     * @param signal string to control behavior of routine (e.g. what
     * data item to return)
     * @return true if component is returning data element, false otherwise
     */
    virtual bool getDataItem(void *data, const char *signal = NULL);

    /**
     * Set rank holding the component
     * @param rank processor rank holding the component
     */
    void setRank(int rank);

    /**
     * Get rank holding the component
     * @return processor rank holding the component
     */
    int getRank(void) const;

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

    /**
     * Rank holding the component. Useful for debugging
     */
    int p_rank;

  private:

  friend class boost::serialization::access;

  /// Serialization "pack" routine
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const
  {
    // p_XCBuf and p_XCBufSize are managed somewhere else; they will
    // have to be initialized
    ar << boost::serialization::base_object<BaseMatrixInterface>(*this)
       << boost::serialization::base_object<GenMatVecInterface>(*this)
       << p_XCBufSize
       << p_mode;
  }

  /// Serialization "unpack" routine
  template<class Archive>
  void load(Archive & ar, const unsigned int version)
  {
    ar >> boost::serialization::base_object<BaseMatrixInterface>(*this)
       >> boost::serialization::base_object<GenMatVecInterface>(*this)
       >> p_XCBufSize
       >> p_mode;
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()

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
     * @param branch pointer to a branch that is connected to bus
     */
    void addBranch(const boost::shared_ptr<BaseComponent> & branch);

    /**
     * Add a pointer to the list of buses that a bus is connected to via
     * a branch
     * @param bus pointer to a bus that is connected to bus via a branch
     */
    void addBus(const boost::shared_ptr<BaseComponent> & bus);

    /**
     * Get pointers to branches that are connected to bus
     * @param nghbrs list of pointers to neighboring branches
     */
    void getNeighborBranches(std::vector<boost::shared_ptr<BaseComponent> > &nghbrs) const;

    /**
     * Get pointers to buses that are connected to calling bus via a branch
     * @param nghbrs list of pointers to neighboring buses
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
     * @param status reference bus status
     */
    void setReferenceBus(bool status);

    /**
     * Get reference bus status
     * @return true if bus is reference bus, false otherwise
     */
    bool getReferenceBus(void) const;

    /**
     * Set original index (from input file)
     * @param idx original index from network
     */
    void setOriginalIndex(int idx);

    /**
     * Get original index
     * @return original index from network
     */
    int getOriginalIndex(void) const;

    /**
     * Set global index
     * @param idx global index from network
     */
    void setGlobalIndex(int idx);

    /**
     * Get global index
     * @return global index from network
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
       & p_refBus
       & p_originalIndex
       & p_globalIndex;

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
     * @param bus pointer to bus
     */
    void setBus1(const boost::shared_ptr<BaseComponent> & bus);

    /**
     * Set pointer to bus at other end of branch
     * @param bus pointer to bus
     */
    void setBus2(const boost::shared_ptr<BaseComponent> & bus);

    /**
     * Get pointer to bus at one end of branch
     * @return pointer to bus 1
     */
    boost::shared_ptr<BaseComponent> getBus1(void) const;

    /**
     * Get pointer to bus at other end of branch
     * @return pointer to bus 2
     */
    boost::shared_ptr<BaseComponent> getBus2(void) const;

    /**
     * Clear bus pointers
     */
    void clearBuses(void);

    /**
     * Set global index for branch
     * @param idx global index of branch
     */
    void setGlobalIndex(int idx);

    /**
     * Set original index for bus 1
     * @param idx original index for bus 1 (assigned from input file)
     */
    void setBus1OriginalIndex(int idx);

    /**
     * Set original index for bus 2
     * @param idx original index for bus 2 (assigned from input file)
     */
    void setBus2OriginalIndex(int idx);

    /**
     * Set global index (from network) for bus 1
     * @param idx global index for bus 1
     */
    void setBus1GlobalIndex(int idx);

    /**
     * Set global index (from network) for bus 2
     * @param idx global index for bus 2
     */
    void setBus2GlobalIndex(int idx);

    /**
     * Get original index for bus 1
     * @return original index for bus 1
     */
    int getBus1OriginalIndex(void) const;

    /**
     * Get original index for bus 2
     * @return original index for bus 2
     */
    int getBus2OriginalIndex(void) const;

    /**
     * Get global index for bus 1
     * @return global index for bus 1
     */
    int getBus1GlobalIndex(void) const;

    /**
     * Get global index for bus 2
     * @return global index for bus 2
     */
    int getBus2GlobalIndex(void) const;

    /**
     * Get global index for branch
     */
    int getGlobalIndex(void) const;

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
    int p_globalIndex;
    int p_globalBus1Index;
    int p_globalBus2Index;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
     ar & boost::serialization::base_object<BaseComponent>(*this)
       & p_originalBus1Index
       & p_originalBus2Index
       & p_globalBus1Index
       & p_globalBus2Index;
  }

};

}    // component
}    // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseComponent)
BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseBusComponent)
BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseBranchComponent)


#endif
