/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   bus_component_c.cpp
 * @author Bruce Palmer
 * @date   2014-08-15 11:05:08 d3g293
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
#include "fortran_component.hpp"

typedef gridpack::fortran_component::FortranBusComponent FortranBus;

struct busWrapper {
  boost::shared_ptr<FortranBus> bus;
};


/**
 * Return size of matrix block on the diagonal contributed by component
 * @param bus GridPACK bus object
 * @param isize,jsize number of rows and columns of matrix block
 * @return false if network component does not contribute matrix element
 */
//extern "C" bool bus_matrix_diag_size(busWrapper *wbus, int *isize, int *jsize)
//{
//  return wbus->bus->matrixDiagSize(isize, jsize);
//}
