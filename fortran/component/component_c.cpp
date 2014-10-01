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
#include <stdio.h>

typedef gridpack::fortran_component::FortranBusComponent FortranBus;

struct busWrapper {
  FortranBus *bus;
};

typedef gridpack::fortran_component::FortranBranchComponent FortranBranch;

struct branchWrapper {
  FortranBranch *branch;
};

/**
 * Get the matrix index for component, based on location of
 * component in network
 * @param bus GridPACK bus object
 * @param idx matrix index of bus
 */
extern "C" void p_bus_get_mat_vec_index(FortranBus *bus, int *idx)
{
  return bus->getMatVecIndex(idx);
}

/**
 * Return total number of neighboring branches/buses attached to bus
 * @param bus bus object pointer
 * @return number of neighboring branches/buses
 */
extern "C" int p_bus_get_num_neighbors(FortranBus *bus)
{
  return bus->getNumNeighbors();
}

/**
 * Get pointer to bus that is attached to calling bus via a branch
 * @param bus bus object pointer
 * @param idx index of neighboring bus (value is between 0 and number of
 * neighbors -1)
 * @return pointer to bus
 */
extern "C" void* p_bus_get_neighbor_bus(FortranBus *bus, int idx)
{
  idx--;
  return static_cast<void*>(bus->getNeighborBus(idx));
}

/**
 * Get pointer to branch that is attached to calling bus
 * @param bus bus object pointer
 * @param idx index of neighboring branch (value is between 0 and number of
 * neighbors -1)
 * @return pointer to branch
 */
extern "C" void* p_bus_get_neighbor_branch(FortranBus *bus, int idx)
{
  idx--;
  return static_cast<void*>(bus->getNeighborBranch(idx));
}

/**
 * Clear all pointers to neighboring branches
 * @param bus bus object pointer
 */
extern "C" void p_bus_clear_branches(FortranBus *bus)
{
  bus->clearBranches();
}

/**
 * Clear all pointers to neighboring buses
 * @param bus bus object pointer
 */
extern "C" void p_bus_clear_buses(FortranBus *bus)
{
  bus->clearBuses();
}

/**
 * Set reference bus status
 * @param bus bus object pointer
 * @param status true if bus is reference bus
 */
extern "C" void p_bus_set_reference_bus(FortranBus *bus, bool status)
{
  bus->setReferenceBus(status);
}

/**
 * Get reference bus status
 * @param bus bus object pointer
 * @return true if bus is reference bus
 */
extern "C" bool p_bus_get_reference_bus(FortranBus *bus)
{
  return bus->getReferenceBus();
}

/**
 * Return the original index (from input file) of bus
 * @param bus bus object pointer
 * @return original index
 */
extern "C" int p_bus_get_original_index(FortranBus *bus)
{
  return bus->getOriginalIndex();
}

/**
 * Return the global index of bus
 * @param bus bus object pointer
 * @return global index
 */
extern "C" int p_bus_get_global_index(FortranBus *bus)
{
  return bus->getGlobalIndex();
}

/**
 * Clear all pointers to buses at each end of branch
 * @param branch branch object pointer
 */
extern "C" void p_branch_clear_buses(FortranBranch *branch)
{
  branch->clearBuses();
}

/**
 * Get the matrix index for component, based on location of
 * component in network
 * @param branch GridPACK branch object
 * @param idx matrix row index of bus
 * @param idx matrix column index of bus
 */
extern "C" void p_branch_get_mat_vec_indices(FortranBranch *branch,
    int *idx, int *jdx)
{
  return branch->getMatVecIndices(idx,jdx);
}

/**
 * Get pointer to bus that is attached to "from" end of branch
 * @param branch branch object pointer
 * @return pointer to branch
 */
extern "C" void* p_branch_get_bus1(FortranBranch *branch)
{
  FortranBus *bus = dynamic_cast<FortranBus*>(branch->getBus1().get());
  return bus->getFortranPointer();
}

/**
 * Get pointer to bus that is attached to "to" end of branch
 * @param branch branch object pointer
 * @return pointer to bus
 */
extern "C" void* p_branch_get_bus2(FortranBranch *branch)
{
  FortranBus *bus = dynamic_cast<FortranBus*>(branch->getBus2().get());
  return bus->getFortranPointer();
}

/**
 * Get original index of "from" bus
 * @param branch branch object pointer
 * @return original index from network
 */
extern "C" int p_branch_get_bus1_original_index(FortranBranch *branch)
{
  return branch->getBus1OriginalIndex();
}

/**
 * Get original index of "to" bus
 * @param branch branch object pointer
 * @return original index from network
 */
extern "C" int p_branch_get_bus2_original_index(FortranBranch *branch)
{
  return branch->getBus2OriginalIndex();
}

/**
 * Get original index of "from" bus
 * @param branch branch object pointer
 * @return original index from network
 */
extern "C" int p_branch_get_bus1_global_index(FortranBranch *branch)
{
  return branch->getBus1GlobalIndex();
}

/**
 * Get original index of "to" bus
 * @param branch branch object pointer
 * @return original index from network
 */
extern "C" int p_branch_get_bus2_global_index(FortranBranch *branch)
{
  return branch->getBus2GlobalIndex();
}
