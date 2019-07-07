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

#ifndef _base_matrix_ifc_h_
#define _base_matrix_ifc_h_

#include <vector>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "boost/smart_ptr/weak_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/component/matvec_ifc.hpp"

#include <boost/serialization/export.hpp>

// TODO: Might want to put MatrixIndices and VectorIndex operations into a
//       separate class since these can probably be implemented once for all
//       network components
namespace gridpack{
namespace component{


class BaseMatrixInterface : public MatVecInterface {
  public:

    enum MapMode{STANDARD};

    /**
     * Constructor
     */
    BaseMatrixInterface(void);

    /**
     * Destructor
     */
    virtual ~BaseMatrixInterface(void);

    /**
     * Return size of matrix block on the diagonal contributed by component
     * @param isize,jsize number of rows and columns of matrix
     *        block
     * @return false if network component does not contribute
     *        matrix element
     */
    virtual bool baseMatrixDiagSize(int *isize, int *jsize) const;

    /**
     * Return the values for a diagonal matrix block. The values are
     * returned in row-major order.
     * @param values pointer to matrix block values
     * @return false if network component does not contribute
     *        matrix element
     */
    virtual bool baseMatrixDiagValues(ComplexType *values);
    virtual bool baseMatrixDiagValues(RealType *values);

    /**
     * Return size of off-diagonal matrix block contributed by component. The
     * values are for the forward direction.
     * @param isize,jsize number of rows and columns of matrix
     *        block
     * @return false if network component does not contribute
     *        matrix element
     */
    virtual bool baseMatrixForwardSize(int *isize, int *jsize) const;

    /**
     * Return the values for an off-diagonl matrix block. The values are
     * for the forward direction and are returned in row-major order.
     * @param values pointer to matrix block values
     * @return false if network component does not contribute
     *        matrix element
     */
    virtual bool baseMatrixForwardValues(ComplexType *values);
    virtual bool baseMatrixForwardValues(RealType *values);

    /**
     * Return size of off-diagonal matrix block contributed by component. The
     * values are for the reverse direction.
     * @param isize,jsize number of rows and columns of matrix
     *        block
     * @return false if network component does not contribute
     *        matrix element
     */
    virtual bool baseMatrixReverseSize(int *isize, int *jsize) const;

    /**
     * Return the values for an off-diagonl matrix block. The values are
     * for the reverse direction and are returned in row-major order.
     * @param values pointer to matrix block values
     * @return false if network component does not contribute
     *        matrix element
     */
    virtual bool baseMatrixReverseValues(ComplexType *values);
    virtual bool baseMatrixReverseValues(RealType *values);

    /**
     * Return size of vector block contributed by component
     * @param isize number of vector elements
     * @return false if network component does not contribute
     *        vector element
     */
    virtual bool baseVectorSize(int *isize) const;

    /**
     * Return the values of the vector block
     * @param values pointer to vector values
     * @return false if network component does not contribute
     *        vector element
     */
    virtual bool baseVectorValues(ComplexType *values);
    virtual bool baseVectorValues(RealType *values);

    /**
     * Set values in the bus or branch component based on values in a vector or
     * matrix
     * @param values values in vector or matrix
     */
    virtual void baseSetValues(ComplexType *values);
    virtual void baseSetValues(RealType *values);

    /**
     * Set the matrix index for diagonal matrix components or vector component,
     * based on location of component in network
     * @param idx value of index
     */
    void baseSetMatVecIndex(int idx);

    /**
     * Get the matrix index for diagonal matrix components or vector component,
     * based on location of component in network
     * @param value of index
     */
    void baseGetMatVecIndex(int *idx) const;

    /**
     * Set the matrix indices for matrix components, based on location of component
     * in network
     * @param idx,jdx value of indices
     */
    void baseSetMatVecIndices(int idx, int jdx);

    /**
     * Get the matrix indices for matrix components, * based on location of component
     * in network
     * @param idx,jdx value of indices
     */
    void baseGetMatVecIndices(int *idx, int *jdx) const;

    /**
     * Set the internal mode variable
     * @param mode current value of mode
     */
    void baseSetMode(int mode);

    //TODO: May need to include routines that support moving values from vectors
    //      back into network components.

  private:
    int p_ival, p_idx, p_jdx;

    int p_mode;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & p_ival
       & p_idx
       & p_jdx
       & p_mode;
  }


};
}
}
#endif
