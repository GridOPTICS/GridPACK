/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_component.hpp
 * @author Bruce Palmer
 * @date   2019-07-12 11:24:20 d3g096
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
    ar & boost::serialization::base_object<MatVecInterface>(*this)& p_ival
       & p_idx
       & p_jdx
       & p_mode;
  }


};

// -------------------------------------------------------------
//  class BaseGenMatVecInterface:
//    A general interface for defining matrices generated from
//    the network that do not follow the mapping of buses to
//    diagonal blocks and branches to off-diagonal blocks
// -------------------------------------------------------------

class BaseGenMatVecInterface : public GenMatVecInterface {
  public:

    enum MapMode{STANDARD};

    /**
     * Constructor
     */
    BaseGenMatVecInterface(void);

    /**
     * Destructor
     */
    virtual ~BaseGenMatVecInterface(void);

    /**
     * Return number of rows in matrix from component
     * @return number of rows from component
     */
    virtual int baseMatrixNumRows() const;

    /**
     * Return number of columns in matrix from component
     * @return number of columnsows from component
     */
    virtual int baseMatrixNumCols() const;

    /**
     * Set row indices corresponding to the rows contributed by this
     * component
     * @param irow index of row contributed by this component (e.g. if component
     * contributes 3 rows then irow is between 0 and 2)
     * @param idx matrix index of row irow
     */
    virtual void baseMatrixSetRowIndex(int irow, int idx);

    /**
     * Set column indices corresponding to the columns contributed by this
     * component
     * @param icol index of column contributed by this component (e.g. if component
     * contributes 3 columns then icol is between 0 and 2)
     * @param idx matrix index of column icol
     */
    virtual void baseMatrixSetColIndex(int icol, int idx);

    /**
     * Get the row indices corresponding to the rows contributed by this component
     * @param irow index of row contributed by this component (e.g. if component
     * contributes 3 rows then irow is between 0 and 2)
     * @return matrix index of row irow
     */
    virtual int baseMatrixGetRowIndex(int irow);

    /**
     * Get the column indices corresponding to the columns contributed by this component
     * @param icol index of column contributed by this component (e.g. if component
     * contributes 3 columns then icol is between 0 and 2)
     * @return matrix index of column icol
     */
    virtual int baseMatrixGetColIndex(int icol);

    /**
     * Return the number of matrix values contributed by this component
     * @return number of matrix values
     */
    virtual int baseMatrixNumValues() const;

    /**
     * Get a list of matrix values contributed by this component and their
     * matrix indices
     * @param values list of matrix element values
     * @param rows row indices for the matrix elements
     * @param cols column indices for the matrix elements
     */
    virtual void baseMatrixGetValues(ComplexType *values, int *rows, int*cols);
    virtual void baseMatrixGetValues(RealType *values, int *rows, int*cols);

    /**
     * Return number of elements in vector from component
     * @return number of elements contributed from component
     */
    virtual int baseVectorNumElements() const;

    /**
     * Set indices corresponding to the elements contributed by this
     * component
     * @param ielem index of element contributed by this component (e.g. if component
     * contributes 3 elements then ielem is between 0 and 2)
     * @param idx vector index of element ielem
     */
    virtual void baseVectorSetElementIndex(int ielem, int idx);

    /**
     * Get list of element indices from component
     * @param idx list of indices that component maps onto
     */
    virtual void baseVectorGetElementIndices(int *idx);

    /**
     * Get a list of vector values contributed by this component and their
     * indices
     * @param values list of vector element values
     * @param idx indices for the vector elements
     */
    virtual void baseVectorGetElementValues(ComplexType *values, int *idx);
    virtual void baseVectorGetElementValues(RealType *values, int *idx);

    /**
     * Transfer vector values to component
     * @param values list of vector element values
     */
    virtual void baseVectorSetElementValues(ComplexType *values);
    virtual void baseVectorSetElementValues(RealType *values);

    /**
     * Return number of rows and columns in matrix from component
     * Number of columns must be the same for all components
     * @return size of block contributed by component
     */
    virtual void baseSlabSize(int *rows, int *cols) const;

    /**
     * Set indices corresponding to the rows contributed by this
     * component
     * @param irow index of row contributed by this component (e.g. if component
     * contributes 3 rows then irow is between 0 and 2)
     * @param idx row index of row irow
     */
    virtual void baseSlabSetRowIndex(int irow, int idx);

    /**
     * Get list of row indices from component
     * @param idx list of row indices that component maps onto
     */
    virtual void baseSlabGetRowIndices(int *idx);

    /**
     * Get a list of row values contributed by this component and their
     * indices
     * @param values list of values for rows
     * @param idx indices for the matrix rows
     */
    virtual void baseSlabGetValues(std::vector<ComplexType*> &values, int *idx);
    virtual void baseSlabGetValues(std::vector<RealType*> &values, int *idx);

    /**
     * Transfer slab values to component
     * @param values list of slab values
     */
    virtual void baseSlabSetValues(ComplexType **values);
    virtual void baseSlabSetValues(RealType **values);

  private:

    int p_mode;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<GenMatVecInterface>(*this)& p_mode;
    }


};

}    // component
}    // gridpack

#endif
