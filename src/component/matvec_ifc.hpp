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

#ifndef _matvec_ifc_h_
#define _matvec_ifc_h_

#include <vector>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "boost/smart_ptr/weak_ptr.hpp"
#include "gridpack/utilities/complex.hpp"

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

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    // ar;
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
     * Get a list of matrix values contributed by this component and their
     * matrix indices
     * @param values list of matrix element values
     * @param rows row indices for the matrix elements
     * @param cols column indices for the matrix elements
     */
    virtual void matrixGetValues(ComplexType *values, int *rows, int*cols);
    virtual void matrixGetValues(RealType *values, int *rows, int*cols);

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

}    // component
}    // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::component::MatVecInterface)
#endif
