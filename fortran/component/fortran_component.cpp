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
#include "fortran_component.hpp"
#include <string.h>
#include <stdio.h>

extern "C" void* bus_allocate(void *ptr);
extern "C" void bus_deallocate(void *ptr);
extern "C" bool p_bus_matrix_diag_size(void *ptr, int *isize, int *jsize);
extern "C" void* p_bus_load(void *ptr, gridpack::component::DataCollection *data);
extern "C" void* branch_allocate(void *ptr);
extern "C" void branch_deallocate(void *ptr);
extern "C" bool p_bus_matrix_diag_values(void *ptr,
    gridpack::ComplexType *values);
extern "C" bool p_bus_matrix_forward_size(void *ptr,
    int *isize, int *jsize);
extern "C" bool p_bus_matrix_forward_values(void *ptr,
    gridpack::ComplexType *values);
extern "C" bool p_bus_matrix_reverse_size(void *ptr,
    int *isize, int *jsize);
extern "C" bool p_bus_matrix_reverse_values(void *ptr,
    gridpack::ComplexType *values);
extern "C" bool p_bus_vector_size(void *ptr,  int *size);
extern "C" bool p_bus_vector_values(void *ptr,
    gridpack::ComplexType *values);
extern "C" void  p_bus_set_values(void *ptr,
    gridpack::ComplexType *values);
extern "C" int p_bus_get_xc_buf_size(void *ptr);
extern "C" int p_bus_get_xc_buf(void *ptr, void **buf);
extern "C" void p_bus_set_mode(void *ptr, int mode);
extern "C" bool  p_bus_serial_write(void *ptr, char *string,
    int bufsize, const char* signal, int signal_len, int *write_len);
extern "C" bool p_branch_matrix_diag_size(void *ptr,
    int *isize, int *jsize);
extern "C" bool p_branch_matrix_diag_values(void *ptr,
    gridpack::ComplexType *values);
extern "C" bool p_branch_matrix_forward_size(void *ptr,
    int *isize, int *jsize);
extern "C" bool p_branch_matrix_forward_values(void *ptr,
    gridpack::ComplexType *values);
extern "C" bool p_branch_matrix_reverse_size(void *ptr,
    int *isize, int *jsize);
extern "C" bool p_branch_matrix_reverse_values(void *ptr,
    gridpack::ComplexType *values);
extern "C" bool p_branch_vector_size(void *ptr,  int *size);
extern "C" bool p_branch_vector_values(void *ptr,
    gridpack::ComplexType *values);
extern "C" void  p_branch_set_values(void *ptr,
    gridpack::ComplexType *values);
extern "C" void p_branch_load(void *ptr, gridpack::component::DataCollection
    *data);
extern "C" int p_branch_get_xc_buf_size(void *ptr);
extern "C" int p_branch_get_xc_buf(void *ptr, void **buf);
extern "C" void p_branch_set_mode(void *ptr, int mode);
extern "C" bool  p_branch_serial_write(void *ptr, char *string,
    int bufsize, const char* signal, int signal_len, int *write_len);

// Base implementation of the MatVecInterface. These functions should be
// overwritten in actual components

namespace gridpack {
namespace fortran_component {

// Base implementation for a Fortran bus object.

/**
 * Simple constructor
 */
FortranBusComponent::FortranBusComponent(void)
{
  p_fortran_bus_ptr = bus_allocate(this);
}

/**
 * Simple destructor
 */
FortranBusComponent::~FortranBusComponent(void)
{
  bus_deallocate(p_fortran_bus_ptr);
}

/**
 * Return size of matrix block on the diagonal contributed by component
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool FortranBusComponent::matrixDiagSize(int *isize, int *jsize) const
{
  return p_bus_matrix_diag_size(p_fortran_bus_ptr, isize, jsize);
}

/**
 * Return the values of for a diagonal matrix block. The values are
 * returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool FortranBusComponent::matrixDiagValues(gridpack::ComplexType *values)
{
  return p_bus_matrix_diag_values(p_fortran_bus_ptr, values);
}

/**
 * Return size of off-diagonal matrix block contributed by component.
 * The values are for the forward direction.
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool FortranBusComponent::matrixForwardSize(int *isize, int *jsize) const
{
  return p_bus_matrix_forward_size(p_fortran_bus_ptr, isize, jsize);
}

/**
 * Return size of off-diagonal matrix block contributed by component.
 * The values are for the reverse direction.
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool FortranBusComponent::matrixReverseSize(int *isize, int *jsize) const
{
  return p_bus_matrix_reverse_size(p_fortran_bus_ptr, isize, jsize);
}

/**
 * Return the values of for an off-diagonl matrix block. The values
 * are for the forward direction and are returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool FortranBusComponent::matrixForwardValues(gridpack::ComplexType *values)
{
  return p_bus_matrix_forward_values(p_fortran_bus_ptr, values);
}

/**
 * Return the values of for an off-diagonl matrix block. The values
 * are for the reverse direction and are returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool FortranBusComponent::matrixReverseValues(gridpack::ComplexType *values)
{
  return p_bus_matrix_reverse_values(p_fortran_bus_ptr, values);
}

/** 
 * Return size of vector block contributed by component 
 * @param isize number of vector elements 
 * @return false if network component does not contribute 
 *        vector element 
 */
bool FortranBusComponent::vectorSize(int *isize) const
{
  return p_bus_vector_size(p_fortran_bus_ptr, isize);
}

/**
 * Set values in the bus component based on values in a
 * vector or matrix
 * @param values values in vector or matrix
 */
void FortranBusComponent::setValues(gridpack::ComplexType *values)
{
  p_bus_set_values(p_fortran_bus_ptr, values);
}

/**
 * Return the values of the vector block
 * @param values pointer to vector values
 * @return false if network component does not contribute
 *        vector element
 */
bool FortranBusComponent::vectorValues(gridpack::ComplexType *values)
{
  return p_bus_vector_values(p_fortran_bus_ptr, values);
}

/**
 * Load data from DataCollection object into corresponding
 * component. This needs to be implemented by every component
 */
void FortranBusComponent::load(boost::shared_ptr
    <gridpack::component::DataCollection> data)
{
  gridpack::component::DataCollection *data_ptr = data.get();
  p_bus_load(p_fortran_bus_ptr, data_ptr);
}

/**
 * Return the size of the buffer needed for data exchanges. Note that this
 * must be the same size for all bus and all branch objects (branch buffers
 * do not need to be the same size as bus buffers), even if all objects
 * do not require the same parameters. Thus, the buffer must be big enough
 * to exchange all variables that an object might need, even if individual
 * objects don't need all the variables
 * @return size of buffer
 */
int FortranBusComponent::getXCBufSize(void)
{
  return p_bus_get_xc_buf_size(p_fortran_bus_ptr);
}

/**
 * Return the location of the data exchange buffer
 * @param buf void pointer to exchange buffer
 */
void FortranBusComponent::getXCBuf(void **buf)
{
  p_bus_get_xc_buf(p_fortran_bus_ptr, buf);
}

/**
 * Set an internal variable that can be used to control the behavior of the
 * component. This function doesn't need to be implemented, but if needed,
 * it can be used to change the behavior of the network in different phases
 * of the calculation. For example, if a different matrix needs to be
 * generated at different times, the mode of the calculation can changed to
 * get different values from the MatVecInterface functions
 * @param mode integer indicating which mode should be used
 */
void FortranBusComponent::setMode(int mode)
{
  p_bus_set_mode(p_fortran_bus_ptr, mode);
}

/**
 * Copy a string for output into buffer. The behavior of this method can be
 * altered by inputting different values for the signal string
 * @param string buffer containing string to be written to output
 * @param bufsize size of string buffer in bytes
 * @param signal string to control behavior of routine (e.g. what
 * properties to write
 * @return true if component is writing a contribution, false otherwise
 */
bool FortranBusComponent::serialWrite(char *string, const int bufsize,
    const char *signal)
{
  int signal_len = strlen(signal);
  int write_len;
  bool ret = p_bus_serial_write(p_fortran_bus_ptr, string, bufsize, signal,
      signal_len, &write_len);
  string[write_len] = '\n';
  string[write_len+1] = '\0';
  return ret;
}

/**
 * Set local index
 * @param idx local index of bus
 */
void FortranBusComponent::setLocalIndex(int idx)
{
  p_local_index = idx;
}

/**
 * Get local index
 * @return local index of bus
 */
int FortranBusComponent::getLocalIndex(void) const
{
  return p_local_index;
}


/**
 * Return the number of neighbors that are attached to bus
 * @return number of branches/buses that are attached to calling bus
 */
int FortranBusComponent::getNumNeighbors()
{
  std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > branches;
  getNeighborBuses(branches);
  return branches.size();
}

/**
 * Return pointer to bus to calling program
 * @param idx neighbor index (runs between 0 and number of neighbors - 1)
 * @return pointer to bus object
 */
void* FortranBusComponent::getNeighborBus(int idx) const
{
  std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > buses;
  getNeighborBuses(buses);
  return dynamic_cast<FortranBusComponent*>(buses[idx].get())->getFortranPointer();
  //return buses[idx].get();
}

/**
 * Return pointer to branch to calling program
 * @param idx neighbor index (runs between 0 and number of neighbors - 1)
 * @return pointer to branch object
 */
void* FortranBusComponent::getNeighborBranch(int idx) const
{
  std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > branches;
  getNeighborBranches(branches);
  return dynamic_cast<FortranBranchComponent*>(branches[idx].get())->getFortranPointer();
  //return branches[idx].get();
}

/**
 * Return pointer to imbedded Fortran object
 * @return pointer to Fortran wrapper
 */
void* FortranBusComponent::getFortranPointer() const
{
  return p_fortran_bus_ptr;
}

// Base implementation for a Fortran branch object.

/**
 * Simple constructor
 */
FortranBranchComponent::FortranBranchComponent(void)
{
  p_fortran_branch_ptr = branch_allocate(this);
}

/**
 * Simple destructor
 */
FortranBranchComponent::~FortranBranchComponent(void)
{
  branch_deallocate(p_fortran_branch_ptr);
}

/**
 * Return size of matrix block on the diagonal contributed by component
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool FortranBranchComponent::matrixDiagSize(int *isize, int *jsize) const
{
  return p_branch_matrix_diag_size(p_fortran_branch_ptr, isize, jsize);
}

/**
 * Return the values of for a diagonal matrix block. The values are
 * returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool FortranBranchComponent::matrixDiagValues(gridpack::ComplexType *values)
{
  return p_branch_matrix_diag_values(p_fortran_branch_ptr, values);
}

/**
 * Return size of off-diagonal matrix block contributed by component.
 * The values are for the forward direction.
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool FortranBranchComponent::matrixForwardSize(int *isize, int *jsize) const
{
  return p_branch_matrix_forward_size(p_fortran_branch_ptr, isize, jsize);
}

/**
 * Return size of off-diagonal matrix block contributed by component.
 * The values are for the reverse direction.
 * @param isize,jsize number of rows and columns of matrix
 *        block
 * @return false if network component does not contribute
 *        matrix element
 */
bool FortranBranchComponent::matrixReverseSize(int *isize, int *jsize) const
{
  return p_branch_matrix_reverse_size(p_fortran_branch_ptr, isize, jsize);
}

/**
 * Return the values of for an off-diagonl matrix block. The values
 * are for the forward direction and are returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool FortranBranchComponent::matrixForwardValues(gridpack::ComplexType *values)
{
  return p_branch_matrix_forward_values(p_fortran_branch_ptr, values);
}

/**
 * Return the values of for an off-diagonl matrix block. The values
 * are for the reverse direction and are returned in row-major order.
 * @param values pointer to matrix block values
 * @return false if network component does not contribute
 *        matrix element
 */
bool FortranBranchComponent::matrixReverseValues(gridpack::ComplexType *values)
{
  return p_branch_matrix_reverse_values(p_fortran_branch_ptr, values);
}

/** 
 * Return size of vector block contributed by component 
 * @param isize number of vector elements 
 * @return false if network component does not contribute 
 *        vector element 
 */
bool FortranBranchComponent::vectorSize(int *isize) const
{
  return p_branch_vector_size(p_fortran_branch_ptr, isize);
}

/**
 * Set values in the branch component based on values in a
 * vector or matrix
 * @param values values in vector or matrix
 */
void FortranBranchComponent::setValues(gridpack::ComplexType *values)
{
  p_branch_set_values(p_fortran_branch_ptr, values);
}

/**
 * Return the values of the vector block
 * @param values pointer to vector values
 * @return false if network component does not contribute
 *        vector element
 */
bool FortranBranchComponent::vectorValues(gridpack::ComplexType *values)
{
  return p_branch_vector_values(p_fortran_branch_ptr, values);
}

/**
 * Load data from DataCollection object into corresponding
 * component. This needs to be implemented by every component
 */
void
FortranBranchComponent::load(boost::shared_ptr
    <gridpack::component::DataCollection> data)
{
  gridpack::component::DataCollection *data_ptr = data.get();
  p_branch_load(p_fortran_branch_ptr, data_ptr);
}

/**
 * Return the size of the buffer needed for data exchanges. Note that this
 * must be the same size for all bus and all branch objects (branch buffers
 * do not need to be the same size as bus buffers), even if all objects
 * do not require the same parameters. Thus, the buffer must be big enough
 * to exchange all variables that an object might need, even if individual
 * objects don't need all the variables
 * @return size of buffer
 */
int FortranBranchComponent::getXCBufSize(void)
{
  return p_branch_get_xc_buf_size(p_fortran_branch_ptr);
}

/**
 * Return the location of the data exchange buffer
 * @param buf void pointer to exchange buffer
 */
void FortranBranchComponent::getXCBuf(void **buf)
{
  p_branch_get_xc_buf(p_fortran_branch_ptr, buf);
}

/**
 * Set an internal variable that can be used to control the behavior of the
 * component. This function doesn't need to be implemented, but if needed,
 * it can be used to change the behavior of the network in different phases
 * of the calculation. For example, if a different matrix needs to be
 * generated at different times, the mode of the calculation can changed to
 * get different values from the MatVecInterface functions
 * @param mode integer indicating which mode should be used
 */
void FortranBranchComponent::setMode(int mode)
{
  p_branch_set_mode(p_fortran_branch_ptr, mode);
}

/**
 * Copy a string for output into buffer. The behavior of this method can be
 * altered by inputting different values for the signal string
 * @param string buffer containing string to be written to output
 * @param bufsize size of string buffer in bytes
 * @param signal string to control behavior of routine (e.g. what
 * properties to write
 * @return true if component is writing a contribution, false otherwise
 */
bool FortranBranchComponent::serialWrite(char *string, const int bufsize,
    const char *signal)
{
  int signal_len = strlen(signal);
  int write_len;
  bool ret = p_branch_serial_write(p_fortran_branch_ptr, string, bufsize, signal,
      signal_len, &write_len);
  string[write_len] = '\n';
  string[write_len+1] = '\0';
  return ret;
}

/**
 * Set local index
 * @param network handle for network
 * @param idx local index of branch
 */
void FortranBranchComponent::setLocalIndex(int idx)
{
  p_local_index = idx;
}

/**
 * Get local index
 * @return local index of branch
 */
int FortranBranchComponent::getLocalIndex(void) const
{
  return p_local_index;
}

/**
 * Return pointer to imbedded Fortran object
 * @return pointer to Fortran wrapper
 */
void* FortranBranchComponent::getFortranPointer() const
{
  return p_fortran_branch_ptr;
}


}  // fortran_component
}  // gridpack
