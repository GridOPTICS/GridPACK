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

#include <stdio.h>
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/base_matrix_ifc.hpp"

// Base implementation of the BaseMatrixInterface. These functions should be
// overwritten in actual components

namespace gridpack {
namespace component {

/**
 * Constructor
 */
BaseMatrixInterface::BaseMatrixInterface(void)
{
  p_mode = STANDARD;
  p_ival = 0;
  p_idx = 0;
  p_jdx = 0;
}

/**
 * Constructor
 */
BaseMatrixInterface::~BaseMatrixInterface(void)
{
}

/**
 * Return size of matrix block on the diagonal contributed by component
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool BaseMatrixInterface::baseMatrixDiagSize( int *isize,
             int *jsize) const
{
  bool ret = false;
  switch (p_mode) {
    case STANDARD:
      ret = matrixDiagSize(isize,jsize);
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Return the values for a diagonal matrix block. The values are
 * returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool BaseMatrixInterface::baseMatrixDiagValues(ComplexType *values)
{
  bool ret = false;
  switch (p_mode) {
    case STANDARD:
      ret = matrixDiagValues(values);
      break;
    default:
      break;
  }
  return ret;
}
bool BaseMatrixInterface::baseMatrixDiagValues(RealType *values)
{
  bool ret = false;
  switch (p_mode) {
    case STANDARD:
      ret = matrixDiagValues(values);
      break;
    default:
      break;
  }
  return ret;
}


/**
 * Return size of off-diagonal matrix block contributed by component. The
 * values are for the forward direction.
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool BaseMatrixInterface::baseMatrixForwardSize(int *isize,
       int *jsize) const
{
  bool ret = false;
  switch (p_mode) {
    case STANDARD:
      ret = matrixForwardSize(isize,jsize);
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Return the values for an off-diagonl matrix block. The values are
 * for the forward direction and are returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool BaseMatrixInterface::baseMatrixForwardValues(ComplexType *values)
{
  bool ret = false;
  switch (p_mode) {
    case STANDARD:
      ret = matrixForwardValues(values);
      break;
    default:
      break;
  }
  return ret;
}
bool BaseMatrixInterface::baseMatrixForwardValues(RealType *values)
{
  bool ret = false;
  switch (p_mode) {
    case STANDARD:
      ret = matrixForwardValues(values);
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Return size of off-diagonal matrix block contributed by component. The
 * values are for the reverse direction.
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool BaseMatrixInterface::baseMatrixReverseSize(int *isize,
       int *jsize) const
{
  bool ret = false;
  switch (p_mode) {
    case STANDARD:
      ret = matrixReverseSize(isize,jsize);
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Return the values for an off-diagonl matrix block. The values are
 * for the reverse direction and are returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool BaseMatrixInterface::baseMatrixReverseValues(ComplexType *values)
{
  bool ret = false;
  switch (p_mode) {
    case STANDARD:
      ret = matrixReverseValues(values);
      break;
    default:
      break;
  }
  return ret;
}
bool BaseMatrixInterface::baseMatrixReverseValues(RealType *values)
{
  bool ret = false;
  switch (p_mode) {
    case STANDARD:
      ret = matrixReverseValues(values);
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Return size of vector block contributed by component
 * @param isize number of vector elements
 * @return false if network component does not contribute
 *        vector element
 */
bool BaseMatrixInterface::baseVectorSize(int *isize) const
{
  *isize = 0;
  bool ret = false;
  switch (p_mode) {
    case STANDARD:
      ret = vectorSize(isize);
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Return the values of the vector block
 * @param values pointer to vector values
 * @return false if network component does not contribute
 *        vector element
 */
bool BaseMatrixInterface::baseVectorValues(ComplexType *values)
{
  bool ret = false;
  switch (p_mode) {
    case STANDARD:
      ret = vectorValues(values);
      break;
    default:
      break;
  }
  return ret;
}
bool BaseMatrixInterface::baseVectorValues(RealType *values)
{
  bool ret = false;
  switch (p_mode) {
    case STANDARD:
      ret = vectorValues(values);
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Set values in the bus or branch component based on values in a vector or
 * matrix
 * @param values values in vector or matrix
 */
void BaseMatrixInterface::baseSetValues(ComplexType *values)
{
  switch (p_mode) {
    case STANDARD:
      setValues(values);
      break;
    default:
      break;
  }
}
void BaseMatrixInterface::baseSetValues(RealType *values)
{
  switch (p_mode) {
    case STANDARD:
      setValues(values);
      break;
    default:
      break;
  }
}

/**
 * Set the matrix index for diagonal matrix components or vector component,
 * based on location of component in network
 * @param idx value of index
 */
void BaseMatrixInterface::baseSetMatVecIndex(int idx)
{
  p_ival = idx;
}

/**
 * Get the matrix index for diagonal matrix components or vector component,
 * based on location of component in network
 * @return value of index
 */
void BaseMatrixInterface::baseGetMatVecIndex(int *idx) const
{
  *idx = p_ival;
}

/**
 * Set the matrix indices for matrix components, based on location of
 * component
 * in network
 * @param idx,jdx value of indices
 */
void BaseMatrixInterface::baseSetMatVecIndices(int idx, int jdx)
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
void BaseMatrixInterface::baseGetMatVecIndices(int *idx, int *jdx) const
{
  *idx = p_idx;
  *jdx = p_jdx;
}

/**
 * Set the internal mode variable
 * @param mode current value of mode
 */
void BaseMatrixInterface::baseSetMode(int mode)
{
  p_mode = mode;
}

// base implementation for the generalized matrix-vector interface

/**
 * Constructor
 */
BaseGenMatVecInterface::BaseGenMatVecInterface(void)
{
  p_mode = STANDARD;
}

/**
 * Destructor
 */
BaseGenMatVecInterface::~BaseGenMatVecInterface(void)
{
}

/**
 * Return number of rows in matrix from component
 * @return number of rows from component
 */
int BaseGenMatVecInterface::baseMatrixNumRows() const
{
  int ret = 0;
  switch (p_mode) {
    case STANDARD:
      ret = matrixNumRows();
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Return number of columns in matrix from component
 * @return number of columnsows from component
 */
int BaseGenMatVecInterface::baseMatrixNumCols() const
{
  int ret = 0;
  switch (p_mode) {
    case STANDARD:
      ret = matrixNumCols();
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Set row indices corresponding to the rows contributed by this
 * component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @param idx matrix index of row irow
 */
void BaseGenMatVecInterface::baseMatrixSetRowIndex(int irow, int idx)
{
  switch (p_mode) {
    case STANDARD:
      matrixSetRowIndex(irow,idx);
      break;
    default:
      break;
  }
}

/**
 * Set column indices corresponding to the columns contributed by this
 * component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @param idx matrix index of column icol
 */
void BaseGenMatVecInterface::baseMatrixSetColIndex(int icol, int idx)
{
  switch (p_mode) {
    case STANDARD:
      matrixSetColIndex(icol,idx);
      break;
    default:
      break;
  }
}

/**
 * Get the row indices corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @return matrix index of row irow
 */
int BaseGenMatVecInterface::baseMatrixGetRowIndex(int irow)
{
  int ret = -1;
  switch (p_mode) {
    case STANDARD:
      ret = matrixGetRowIndex(irow);
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Get the column indices corresponding to the columns contributed by this component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @return matrix index of column icol
 */
int BaseGenMatVecInterface::baseMatrixGetColIndex(int icol)
{
  int ret = -1;
  switch (p_mode) {
    case STANDARD:
      ret = matrixGetColIndex(icol);
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Return the number of matrix values contributed by this component
 * @return number of matrix values
 */
int BaseGenMatVecInterface::baseMatrixNumValues() const
{
  int ret = 0;
  switch (p_mode) {
    case STANDARD:
      ret = matrixNumValues();
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Get a list of matrix values contributed by this component and their
 * matrix indices
 * @param values list of matrix element values
 * @param rows row indices for the matrix elements
 * @param cols column indices for the matrix elements
 */
void BaseGenMatVecInterface::baseMatrixGetValues(ComplexType *values, int *rows, int*cols)
{
  switch (p_mode) {
    case STANDARD:
      matrixGetValues(values, rows, cols);
      break;
    default:
      break;
  }
}
void BaseGenMatVecInterface::baseMatrixGetValues(RealType *values, int *rows, int*cols)
{
  switch (p_mode) {
    case STANDARD:
      matrixGetValues(values, rows, cols);
      break;
    default:
      break;
  }
}

/**
 * Return number of elements in vector from component
 * @return number of elements contributed from component
 */
int BaseGenMatVecInterface::baseVectorNumElements() const
{
  int ret = 0;
  switch (p_mode) {
    case STANDARD:
      vectorNumElements();
      break;
    default:
      break;
  }
  return ret;
}

/**
 * Set indices corresponding to the elements contributed by this
 * component
 * @param ielem index of element contributed by this component
 * (e.g. if component contributes 3 elements then ielem is between
 * 0 and 2)
 * @param idx vector index of element ielem
 */
void BaseGenMatVecInterface::baseVectorSetElementIndex(int ielem, int idx)
{
  switch (p_mode) {
    case STANDARD:
      vectorSetElementIndex(ielem, idx);
      break;
    default:
      break;
  }
}

/**
 * Get list of element indices from component
 * @param idx list of indices that component maps onto
 */
void BaseGenMatVecInterface::baseVectorGetElementIndices(int *idx)
{
  switch (p_mode) {
    case STANDARD:
      vectorGetElementIndices(idx);
      break;
    default:
      break;
  }
}

/**
 * Get a list of vector values contributed by this component and their
 * indices
 * @param values list of vector element values
 * @param idx indices for the vector elements
 */
void BaseGenMatVecInterface::baseVectorGetElementValues(ComplexType *values, int *idx)
{
  switch (p_mode) {
    case STANDARD:
      vectorGetElementValues(values,idx);
      break;
    default:
      break;
  }
}
void BaseGenMatVecInterface::baseVectorGetElementValues(RealType *values, int *idx)
{
  switch (p_mode) {
    case STANDARD:
      vectorGetElementValues(values,idx);
      break;
    default:
      break;
  }
}

/**
 * Transfer vector values to component
 * @param values list of vector element values
 */
void BaseGenMatVecInterface::baseVectorSetElementValues(ComplexType *values)
{
  switch (p_mode) {
    case STANDARD:
      vectorSetElementValues(values);
      break;
    default:
      break;
  }
}
void BaseGenMatVecInterface::baseVectorSetElementValues(RealType *values)
{
  switch (p_mode) {
    case STANDARD:
      vectorSetElementValues(values);
      break;
    default:
      break;
  }
}

/**
 * Return number of rows and columns in matrix from component
 * Number of columns must be the same for all components
 * @return size of block contributed by component
 */
void BaseGenMatVecInterface::baseSlabSize(int *rows, int *cols) const
{
  switch (p_mode) {
    case STANDARD:
      slabSize(rows,cols);
      break;
    default:
      break;
  }
}

/**
 * Set indices corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @param idx row index of row irow
 */
void BaseGenMatVecInterface::baseSlabSetRowIndex(int irow, int idx)
{
  switch (p_mode) {
    case STANDARD:
      slabSetRowIndex(irow,idx);
      break;
    default:
      break;
  }
}

/**
 * Get list of row indices from component
 * @param idx list of row indices that component maps onto
 */
void BaseGenMatVecInterface::baseSlabGetRowIndices(int *idx)
{
  switch (p_mode) {
    case STANDARD:
      slabGetRowIndices(idx);
      break;
    default:
      break;
  }
}

/**
 * Get a list of row values contributed by this component and their
 * indices
 * @param values list of values for rows
 * @param idx indices for the matrix rows
 */
void BaseGenMatVecInterface::baseSlabGetValues(std::vector<ComplexType*> &values,
    int *idx)
{
  switch (p_mode) {
    case STANDARD:
      slabGetValues(values, idx);
      break;
    default:
      break;
  }
}
void BaseGenMatVecInterface::baseSlabGetValues(std::vector<RealType*> &values,
    int *idx)
{
  switch (p_mode) {
    case STANDARD:
      slabGetValues(values, idx);
      break;
    default:
      break;
  }
}

/**
 * Transfer slab values to component
 * @param values list of slab values
 */
void BaseGenMatVecInterface::baseSlabSetValues(ComplexType **values)
{
  switch (p_mode) {
    case STANDARD:
      slabSetValues(values);
      break;
    default:
      break;
  }
}
void BaseGenMatVecInterface::baseSlabSetValues(RealType **values)
{
  switch (p_mode) {
    case STANDARD:
      slabSetValues(values);
      break;
    default:
      break;
  }
}

}   // component
}   // gridpack
