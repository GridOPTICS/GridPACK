/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_component.hpp
 * @author Bruce Palmer
 * @date   2016-07-14 13:50:54 d3g096
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
     * @param isize,jsize number of rows and columns of matrix
     *        block
     * @return false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixDiagSize(int *isize, int *jsize) const;

    /**
     * Return the values for a diagonal matrix block. The values are
     * returned in row-major order.
     * @param values pointer to matrix block values
     * @return false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixDiagValues(ComplexType *values);
    virtual bool matrixDiagValues(RealType *values);

    /**
     * Return size of off-diagonal matrix block contributed by component. The
     * values are for the forward direction.
     * @param isize,jsize number of rows and columns of matrix
     *        block
     * @return false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixForwardSize(int *isize, int *jsize) const;

    /**
     * Return the values for an off-diagonl matrix block. The values are
     * for the forward direction and are returned in row-major order.
     * @param values pointer to matrix block values
     * @return false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixForwardValues(ComplexType *values);
    virtual bool matrixForwardValues(RealType *values);

    /**
     * Return size of off-diagonal matrix block contributed by component. The
     * values are for the reverse direction.
     * @param isize,jsize number of rows and columns of matrix
     *        block
     * @return false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixReverseSize(int *isize, int *jsize) const;

    /**
     * Return the values for an off-diagonl matrix block. The values are
     * for the reverse direction and are returned in row-major order.
     * @param values pointer to matrix block values
     * @return false if network component does not contribute
     *        matrix element
     */
    virtual bool matrixReverseValues(ComplexType *values);
    virtual bool matrixReverseValues(RealType *values);

    /**
     * Return size of vector block contributed by component
     * @param isize number of vector elements
     * @return false if network component does not contribute
     *        vector element
     */
    virtual bool vectorSize(int *isize) const;

    /**
     * Return the values of the vector block
     * @param values pointer to vector values
     * @return false if network component does not contribute
     *        vector element
     */
    virtual bool vectorValues(ComplexType *values);
    virtual bool vectorValues(RealType *values);

    /**
     * Set values in the bus or branch component based on values in a vector or
     * matrix
     * @param values values in vector or matrix
     */
    virtual void setValues(ComplexType *values);
    virtual void setValues(RealType *values);

    /**
     * Set the matrix index for diagonal matrix components or vector component,
     * based on location of component in network
     * @param idx value of index
     */
    void setMatVecIndex(int idx);

    /**
     * Get the matrix index for diagonal matrix components or vector component,
     * based on location of component in network
     * @param value of index
     */
    void getMatVecIndex(int *idx) const;

    /**
     * Set the matrix indices for matrix components, based on location of component
     * in network
     * @param idx,jdx value of indices
     */
    void setMatVecIndices(int idx, int jdx);

    /**
     * Get the matrix indices for matrix components, * based on location of component
     * in network
     * @param idx,jdx value of indices
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
//  class GenMatVecInterface:
//    A general interface for defining matrices generated from
//    the network that do not follow the mapping of buses to
//    diagonal blocks and branches to off-diagonal blocks
// -------------------------------------------------------------

class GenMatVecInterface {
  public:
    /**
     * Constructor
     */
    GenMatVecInterface(void);

    /**
     * Destructor
     */
    virtual ~GenMatVecInterface(void);

    /**
     * Return number of rows in matrix from component
     * @return number of rows from component
     */
    virtual int matrixNumRows() const;

    /**
     * Return number of columns in matrix from component
     * @return number of columnsows from component
     */
    virtual int matrixNumCols() const;

    /**
     * Set row indices corresponding to the rows contributed by this
     * component
     * @param irow index of row contributed by this component (e.g. if component
     * contributes 3 rows then irow is between 0 and 2)
     * @param idx matrix index of row irow
     */
    virtual void matrixSetRowIndex(int irow, int idx);

    /**
     * Set column indices corresponding to the columns contributed by this
     * component
     * @param icol index of column contributed by this component (e.g. if component
     * contributes 3 columns then icol is between 0 and 2)
     * @param idx matrix index of column icol
     */
    virtual void matrixSetColIndex(int icol, int idx);

    /**
     * Get the row indices corresponding to the rows contributed by this component
     * @param irow index of row contributed by this component (e.g. if component
     * contributes 3 rows then irow is between 0 and 2)
     * @return matrix index of row irow
     */
    virtual int matrixGetRowIndex(int irow);

    /**
     * Get the column indices corresponding to the columns contributed by this component
     * @param icol index of column contributed by this component (e.g. if component
     * contributes 3 columns then icol is between 0 and 2)
     * @return matrix index of column icol
     */
    virtual int matrixGetColIndex(int icol);

    /**
     * Return the number of matrix values contributed by this component
     * @return number of matrix values
     */
    virtual int matrixNumValues() const;

  /**
   * Return values from a matrix block
   * @param nvals: number of values to be inserted
   * @param values: pointer to matrix block values
   * @param rows: pointer to matrix block rows
   * @param cols: pointer to matrix block cols
   */
  virtual void matrixGetValues(int *nvals,ComplexType *values, int *rows, int*cols);
  virtual void matrixGetValues(int *nvals,RealType *values, int *rows, int*cols);

    /**
     * Return number of elements in vector from component
     * @return number of elements contributed from component
     */
    virtual int vectorNumElements() const;

    /**
     * Set indices corresponding to the elements contributed by this
     * component
     * @param ielem index of element contributed by this component (e.g. if component
     * contributes 3 elements then ielem is between 0 and 2)
     * @param idx vector index of element ielem
     */
    virtual void vectorSetElementIndex(int ielem, int idx);

    /**
     * Get list of element indices from component
     * @param idx list of indices that component maps onto
     */
    virtual void vectorGetElementIndices(int *idx);

    /**
     * Get a list of vector values contributed by this component and their
     * indices
     * @param values list of vector element values
     * @param idx indices for the vector elements
     */
    virtual void vectorGetElementValues(ComplexType *values, int *idx);
    virtual void vectorGetElementValues(RealType *values, int *idx);

    /**
     * Transfer vector values to component
     * @param values list of vector element values
     */
    virtual void vectorSetElementValues(ComplexType *values);
    virtual void vectorSetElementValues(RealType *values);

    /**
     * Return number of rows and columns in matrix from component
     * Number of columns must be the same for all components
     * @return size of block contributed by component
     */
    virtual void slabSize(int *rows, int *cols) const;

    /**
     * Set indices corresponding to the rows contributed by this
     * component
     * @param irow index of row contributed by this component (e.g. if component
     * contributes 3 rows then irow is between 0 and 2)
     * @param idx row index of row irow
     */
    virtual void slabSetRowIndex(int irow, int idx);

    /**
     * Get list of row indices from component
     * @param idx list of row indices that component maps onto
     */
    virtual void slabGetRowIndices(int *idx);

    /**
     * Get a list of row values contributed by this component and their
     * indices
     * @param values list of values for rows
     * @param idx indices for the matrix rows
     */
    virtual void slabGetValues(std::vector<ComplexType*> &values, int *idx);
    virtual void slabGetValues(std::vector<RealType*> &values, int *idx);

    /**
     * Transfer slab values to component
     * @param values list of slab values
     */
    virtual void slabSetValues(ComplexType **values);
    virtual void slabSetValues(RealType **values);

  private:

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    // yar;
  }


};

// -------------------------------------------------------------
//  class BaseComponent:
//  This class implements some basic functions that can be
//  expected from any component on the network.
// -------------------------------------------------------------
class BaseComponent
  : public MatVecInterface, public GenMatVecInterface {
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
     * Save state variables inside the component to a DataCollection object.
     * This can be used as a way of moving data in a way that is useful for
     * creating output or for copying state data from one network to another.
     * @param data data collection object into which new values are inserted
	 * added by Renke, also modify the original bus mag, ang, 
	 * and the original generator PG QG in the datacollection
     */
    virtual void saveDataAlsotoOrg(boost::shared_ptr<gridpack::component::DataCollection>
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
    ar << boost::serialization::base_object<MatVecInterface>(*this)
       << boost::serialization::base_object<GenMatVecInterface>(*this)
       << p_XCBufSize
       << p_mode;
  }

  /// Serialization "unpack" routine
  template<class Archive>
  void load(Archive & ar, const unsigned int version)
  {
    ar >> boost::serialization::base_object<MatVecInterface>(*this)
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

BOOST_CLASS_EXPORT_KEY(gridpack::component::MatVecInterface)
BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseComponent)
BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseBusComponent)
BOOST_CLASS_EXPORT_KEY(gridpack::component::BaseBranchComponent)


#endif
