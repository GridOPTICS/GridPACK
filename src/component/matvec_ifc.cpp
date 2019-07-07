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

#include "gridpack/component/base_component.hpp"
#include "gridpack/component/base_matrix_ifc.hpp"

// Base implementation of the MatVecInterface. These functions should be
// overwritten in actual components

namespace gridpack {
namespace component {
/**
 * Constructor
 */
MatVecInterface::MatVecInterface(void)
{
}

/**
 * Constructor
 */
MatVecInterface::~MatVecInterface(void)
{
}

/**
 * Return size of matrix block on the diagonal contributed by component
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool MatVecInterface::matrixDiagSize( int *isize,
             int *jsize) const
{
  *isize = 0;
  *jsize = 0;
  return false;
}

/**
 * Return the values for a diagonal matrix block. The values are
 * returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool MatVecInterface::matrixDiagValues(ComplexType *values)
{
  return false;
}
bool MatVecInterface::matrixDiagValues(RealType *values)
{
  return false;
}


/**
 * Return size of off-diagonal matrix block contributed by component. The
 * values are for the forward direction.
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool MatVecInterface::matrixForwardSize(int *isize,
       int *jsize) const
{
  *isize = 0;
  *jsize = 0;
  return false;
}

/**
 * Return the values for an off-diagonl matrix block. The values are
 * for the forward direction and are returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool MatVecInterface::matrixForwardValues(ComplexType *values)
{
  return false;
}
bool MatVecInterface::matrixForwardValues(RealType *values)
{
  return false;
}

/**
 * Return size of off-diagonal matrix block contributed by component. The
 * values are for the reverse direction.
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool MatVecInterface::matrixReverseSize(int *isize,
       int *jsize) const
{
  *isize = 0;
  *jsize = 0;
  return false;
}

/**
 * Return the values for an off-diagonl matrix block. The values are
 * for the reverse direction and are returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool MatVecInterface::matrixReverseValues(ComplexType *values)
{
  return false;
}
bool MatVecInterface::matrixReverseValues(RealType *values)
{
  return false;
}

/**
 * Return size of vector block contributed by component
 * @param isize number of vector elements
 * @return false if network component does not contribute
 *        vector element
 */
bool MatVecInterface::vectorSize(int *isize) const
{
  *isize = 0;
  return false;
}

/**
 * Return the values of the vector block
 * @param values pointer to vector values
 * @return false if network component does not contribute
 *        vector element
 */
bool MatVecInterface::vectorValues(ComplexType *values)
{
  return false;
}
bool MatVecInterface::vectorValues(RealType *values)
{
  return false;
}

/**
 * Set values in the bus or branch component based on values in a vector or
 * matrix
 * @param values values in vector or matrix
 */
void MatVecInterface::setValues(ComplexType *values)
{
}
void MatVecInterface::setValues(RealType *values)
{
}

// base implementation for the generalized matrix-vector interface

/**
 * Constructor
 */
GenMatVecInterface::GenMatVecInterface(void)
{
}

/**
 * Destructor
 */
GenMatVecInterface::~GenMatVecInterface(void)
{
}

/**
 * Return number of rows in matrix from component
 * @return number of rows from component
 */
int GenMatVecInterface::matrixNumRows() const
{
  return 0;
}

/**
 * Return number of columns in matrix from component
 * @return number of columnsows from component
 */
int GenMatVecInterface::matrixNumCols() const
{
  return 0;
}

/**
 * Set row indices corresponding to the rows contributed by this
 * component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @param idx matrix index of row irow
 */
void GenMatVecInterface::matrixSetRowIndex(int irow, int idx)
{
}

/**
 * Set column indices corresponding to the columns contributed by this
 * component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @param idx matrix index of column icol
 */
void GenMatVecInterface::matrixSetColIndex(int icol, int idx)
{
}

/**
 * Get the row indices corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @return matrix index of row irow
 */
int GenMatVecInterface::matrixGetRowIndex(int irow)
{
  return -1;
}

/**
 * Get the column indices corresponding to the columns contributed by this component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @return matrix index of column icol
 */
int GenMatVecInterface::matrixGetColIndex(int icol)
{
  return -1;
}

/**
 * Return the number of matrix values contributed by this component
 * @return number of matrix values
 */
int GenMatVecInterface::matrixNumValues() const
{
  return 0;
}

/**
 * Get a list of matrix values contributed by this component and their
 * matrix indices
 * @param values list of matrix element values
 * @param rows row indices for the matrix elements
 * @param cols column indices for the matrix elements
 */
void GenMatVecInterface::matrixGetValues(ComplexType *values, int *rows, int*cols)
{
}
void GenMatVecInterface::matrixGetValues(RealType *values, int *rows, int*cols)
{
}

/**
 * Return number of elements in vector from component
 * @return number of elements contributed from component
 */
int GenMatVecInterface::vectorNumElements() const
{
  return 0;
}

/**
 * Set indices corresponding to the elements contributed by this
 * component
 * @param ielem index of element contributed by this component
 * (e.g. if component contributes 3 elements then ielem is between
 * 0 and 2)
 * @param idx vector index of element ielem
 */
void GenMatVecInterface::vectorSetElementIndex(int ielem, int idx)
{
}

/**
 * Get list of element indices from component
 * @param idx list of indices that component maps onto
 */
void GenMatVecInterface::vectorGetElementIndices(int *idx)
{
}

/**
 * Get a list of vector values contributed by this component and their
 * indices
 * @param values list of vector element values
 * @param idx indices for the vector elements
 */
void GenMatVecInterface::vectorGetElementValues(ComplexType *values, int *idx)
{
}
void GenMatVecInterface::vectorGetElementValues(RealType *values, int *idx)
{
}

/**
 * Transfer vector values to component
 * @param values list of vector element values
 */
void GenMatVecInterface::vectorSetElementValues(ComplexType *values)
{
}
void GenMatVecInterface::vectorSetElementValues(RealType *values)
{
}

/**
 * Return number of rows and columns in matrix from component
 * Number of columns must be the same for all components
 * @return size of block contributed by component
 */
void GenMatVecInterface::slabSize(int *rows, int *cols) const
{
  *rows = 0;
  *cols = 0;
}

/**
 * Set indices corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @param idx row index of row irow
 */
void GenMatVecInterface::slabSetRowIndex(int irow, int idx)
{
}

/**
 * Get list of row indices from component
 * @param idx list of row indices that component maps onto
 */
void GenMatVecInterface::slabGetRowIndices(int *idx)
{
}

/**
 * Get a list of row values contributed by this component and their
 * indices
 * @param values list of values for rows
 * @param idx indices for the matrix rows
 */
void GenMatVecInterface::slabGetValues(std::vector<ComplexType*> &values, int *idx)
{
}
void GenMatVecInterface::slabGetValues(std::vector<RealType*> &values, int *idx)
{
}

/**
 * Transfer slab values to component
 * @param values list of slab values
 */
void GenMatVecInterface::slabSetValues(ComplexType **values)
{
}
void GenMatVecInterface::slabSetValues(RealType **values)
{
}

}  // component
}  // gridpack
