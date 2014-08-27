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
  FortranBus *bus;
};

/**
 * @param wbus bus object wrapper
 * @return number of neighboring branches/buses
 */
extern "C" int bus_get_num_neighbors(busWrapper *wbus)
{
  return wbus->bus->getNumNeighbors();
}

/**
 * Get pointer to bus that is attached to calling bus via a branch
 * @param wbus bus object wrapper
 * @param idx index of neighboring bus (value is between 0 and number of
 * neighbors -1)
 * @return pointer to bus wrapper
 */
extern "C" void* bus_get_neighbor_bus(busWrapper *wbus, int idx)
{
  return static_cast<void*>(wbus->bus->getNeighborBus(idx));
}

/**
 * Get pointer to branch that is attached to calling bus
 * @param wbus bus object wrapper
 * @param idx index of neighboring branch (value is between 0 and number of
 * neighbors -1)
 * @return pointer to branch wrapper
 */
extern "C" void* bus_get_neighbor_branch(busWrapper *wbus, int idx)
{
  return static_cast<void*>(wbus->bus->getNeighborBranch(idx));
}
