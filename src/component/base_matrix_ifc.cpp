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
  printf(" set matvec indices: %d %d\n",p_idx,p_jdx);
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
  printf(" get matvec indices: %d %d\n",*idx,*jdx);
}

/**
 * Set the internal mode variable
 * @param mode current value of mode
 */
void BaseMatrixInterface::baseSetMode(int mode)
{
  p_mode = mode;
}

}   // component
}   // gridpack
