/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   fortran_component.hpp
 * @author Bruce Palmer
 * @date   2013-10-07 11:05:08 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _fortran_component_h_
#define _fortran_component_h_

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"

namespace gridpack{
namespace fortran_component{

class FortranBusComponent
  : public gridpack::component::BaseBusComponent {
  public:
    /**
     * Simple constructor
     */
    FortranBusComponent(void);

    /**
     * Simple destructor
     */
    virtual ~FortranBusComponent(void);

    /**
     * Return size of matrix block on the diagonal contributed by component
     * @param isize,jsize number of rows and columns of matrix
     *        block
     * @return false if network component does not contribute
     *        matrix element
     */
    bool matrixDiagSize(int *isize, int *jsize) const;

    /**
     * Return the values of for a diagonal matrix block. The values are
     * returned in row-major order.
     * @param values pointer to matrix block values
     * @return false if network component does not contribute
     *        matrix element
     */
    bool matrixDiagValues(ComplexType *values);

    /**
     * Return size of off-diagonal matrix block contributed by component. The
     * values are for the forward direction.
     * @param isize,jsize number of rows and columns of matrix
     *        block
     * @return false if network component does not contribute
     *        matrix element
     */
    bool matrixForwardSize(int *isize, int *jsize) const;

    /**
     * Return the values of for an off-diagonl matrix block. The values are
     * for the forward direction and are returned in row-major order.
     * @param values pointer to matrix block values
     * @return false if network component does not contribute
     *        matrix element
     */
    bool matrixForwardValues(ComplexType *values);

    /**
     * Return size of off-diagonal matrix block contributed by component. The
     * values are for the reverse direction.
     * @param isize,jsize number of rows and columns of matrix
     *        block
     * @return false if network component does not contribute
     *        matrix element
     */
    bool matrixReverseSize(int *isize, int *jsize) const;

    /**
     * Return the values of for an off-diagonl matrix block. The values are
     * for the reverse direction and are returned in row-major order.
     * @param values pointer to matrix block values
     * @return false if network component does not contribute
     *        matrix element
     */
    bool matrixReverseValues(ComplexType *values);

    /** 
     * Return size of vector block contributed by component 
     * @param isize number of vector elements 
     * @return false if network component does not contribute 
     *        vector element 
     */ 
    bool vectorSize(int *isize) const;

    /**
     * Return the values of the vector block
     * @param values pointer to vector values
     * @return false if network component does not contribute
     *        vector element
     */
    bool vectorValues(ComplexType *values);

    /**
     * Set values in the bus or branch component based on values in a
     * vector or matrix
     * @param values values in vector or matrix
     */
    void setValues(ComplexType *values);

    /**
     * Load data from DataCollection object into corresponding
     * component. This needs to be implemented by every component
     * @param n_handle index of network object
     * @param idx index of data collection object
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection> data);

    /**
     * Return the size of the buffer needed for data exchanges. Note that this
     * must be the same size for all bus and all branch objects (branch buffers
     * do not need to be the same size as bus buffers), even if all objects
     * do not require the same parameters. Thus, the buffer must be big enough
     * to exchange all variables that an object might need, even if individual
     * objects don't need all the variables
     * @return size of buffer
     */
    int getXCBufSize(void);

    /**
     * Return the location of the data exchange buffer
     * @param buf void pointer to exchange buffer
     */
    void getXCBuf(void **bus);

    /**
     * Set an internal variable that can be used to control the behavior of the
     * component. This function doesn't need to be implemented, but if needed,
     * it can be used to change the behavior of the network in different phases
     * of the calculation. For example, if a different matrix needs to be
     * generated at different times, the mode of the calculation can changed to
     * get different values from the MatVecInterface functions
     * @param mode integer indicating which mode should be used
     */
    void setMode(int mode);

    /**
     * Copy a string for output into buffer. The behavior of this method can be
     * altered by inputting different values for the signal string
     * @param string buffer containing string to be written to output
     * @param bufsize size of string buffer in bytes
     * @param signal string to control behavior of routine (e.g. what
     * properties to write
     * @return true if component is writing a contribution, false otherwise
     */
    bool serialWrite(char *string, const int bufsize, const char
        *signal = NULL);

    /**
     * Return the number of neighbors that are attached to bus
     * @return number of branches/buses that are attached to calling bus
     */
    int getNumNeighbors();

    /**
     * Set local index
     * @param idx local index of bus
     */
    void setLocalIndex(int idx);

    /**
     * Get local index
     * @return local index of bus
     */
    int getLocalIndex(void) const;

    /**
     * Return pointer to bus to calling program
     * @param idx neighbor index (runs between 0 and number of neighbors - 1)
     * @return pointer to bus object
     */
    void* getNeighborBus(int idx) const;

    /**
     * Return pointer to branch to calling program
     * @param idx neighbor index (runs between 0 and number of neighbors - 1)
     * @return pointer to branch object
     */
    void* getNeighborBranch(int idx) const;

    /**
     * Return pointer to imbedded Fortran object
     * @return pointer to Fortran wrapper
     */
    void* getFortranPointer() const;
  private:

    int p_local_index;
    void* p_fortran_bus_ptr;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
     ar & boost::serialization::base_object<BaseBusComponent>(*this)
        & p_local_index;
  }
};

class FortranBranchComponent
  : public gridpack::component::BaseBranchComponent {
  public:
    /**
     * Simple constructor
     */
    FortranBranchComponent(void);

    /**
     * Simple destructor
     */
    virtual ~FortranBranchComponent(void);

    /**
     * Return size of matrix block on the diagonal contributed by component
     * @param isize,jsize number of rows and columns of matrix
     *        block
     * @return false if network component does not contribute
     *        matrix element
     */
    bool matrixDiagSize(int *isize, int *jsize) const;

    /**
     * Return the values of for a diagonal matrix block. The values are
     * returned in row-major order.
     * @param values pointer to matrix block values
     * @return false if network component does not contribute
     *        matrix element
     */
    bool matrixDiagValues(ComplexType *values);

    /**
     * Return size of off-diagonal matrix block contributed by component. The
     * values are for the forward direction.
     * @param isize,jsize number of rows and columns of matrix
     *        block
     * @return false if network component does not contribute
     *        matrix element
     */
    bool matrixForwardSize(int *isize, int *jsize) const;

    /**
     * Return the values of for an off-diagonl matrix block. The values are
     * for the forward direction and are returned in row-major order.
     * @param values pointer to matrix block values
     * @return false if network component does not contribute
     *        matrix element
     */
    bool matrixForwardValues(ComplexType *values);

    /**
     * Return size of off-diagonal matrix block contributed by component. The
     * values are for the reverse direction.
     * @param isize,jsize number of rows and columns of matrix
     *        block
     * @return false if network component does not contribute
     *        matrix element
     */
    bool matrixReverseSize(int *isize, int *jsize) const;

    /**
     * Return the values of for an off-diagonl matrix block. The values are
     * for the reverse direction and are returned in row-major order.
     * @param values pointer to matrix block values
     * @return false if network component does not contribute
     *        matrix element
     */
    bool matrixReverseValues(ComplexType *values);

    /** 
     * Return size of vector block contributed by component 
     * @param isize number of vector elements 
     * @return false if network component does not contribute 
     *        vector element 
     */ 
    bool vectorSize(int *isize) const;

    /**
     * Return the values of the vector block
     * @param values pointer to vector values
     * @return false if network component does not contribute
     *        vector element
     */
    bool vectorValues(ComplexType *values);

    /**
     * Set values in the bus or branch component based on values in a
     * vector or matrix
     * @param values values in vector or matrix
     */
    void setValues(ComplexType *values);

    /**
     * Load data from DataCollection object into corresponding
     * component. This needs to be implemented by every component
     * @param data data collection associated with component
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection> data);

    /**
     * Return the size of the buffer needed for data exchanges. Note that this
     * must be the same size for all bus and all branch objects (branch buffers
     * do not need to be the same size as bus buffers), even if all objects
     * do not require the same parameters. Thus, the buffer must be big enough
     * to exchange all variables that an object might need, even if individual
     * objects don't need all the variables
     * @return size of buffer
     */
    int getXCBufSize(void);

    /**
     * Return the location of the data exchange buffer
     * @param buf void pointer to exchange buffer
     */
    void getXCBuf(void **bus);

    /**
     * Set an internal variable that can be used to control the behavior of the
     * component. This function doesn't need to be implemented, but if needed,
     * it can be used to change the behavior of the network in different phases
     * of the calculation. For example, if a different matrix needs to be
     * generated at different times, the mode of the calculation can changed to
     * get different values from the MatVecInterface functions
     * @param mode integer indicating which mode should be used
     */
    void setMode(int mode);

    /**
     * Copy a string for output into buffer. The behavior of this method can be
     * altered by inputting different values for the signal string
     * @param string buffer containing string to be written to output
     * @param bufsize size of string buffer in bytes
     * @param signal string to control behavior of routine (e.g. what
     * properties to write
     * @return true if component is writing a contribution, false otherwise
     */
    bool serialWrite(char *string, const int bufsize, const char
        *signal = NULL);

    /**
     * Get local index to bus at one end of branch
     * @return local index to bus 1
     */
    int getBus1Index(void) const;

    /**
     * Get local index to bus at other end of branch
     * @return local index to bus 2
     */
    int getBus2Index(void) const;

    /**
     * Set local index
     * @param network handle for network
     * @param idx local index of branch
     */
    void setLocalIndex(int idx);

    /**
     * Get local index
     * @return local index of branch
     */
    int getLocalIndex(void) const;

    /**
     * Return pointer to imbedded Fortran object
     * @return pointer to Fortran wrapper
     */
    void* getFortranPointer() const;
  private:

    int p_local_index;
    void* p_fortran_branch_ptr;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
     ar & boost::serialization::base_object<BaseComponent>(*this)
        & p_local_index;
  }

};

}    // fortran_component
}    // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::fortran_component::FortranBusComponent);
BOOST_CLASS_EXPORT_KEY(gridpack::fortran_component::FortranBranchComponent);
#endif
