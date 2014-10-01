/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   mapper_c.cpp
 * @author Bruce Palmer
 * @date   2014-09-15 d3g293
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
#include <boost/smart_ptr/shared_ptr.hpp>
#include "gridpack/mapper/full_map.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/gen_matrix_map.hpp"
#include "gridpack/mapper/gen_vector_map.hpp"
#include "../component/fortran_component.hpp"

typedef gridpack::network::BaseNetwork<
  gridpack::fortran_component::FortranBusComponent,
  gridpack::fortran_component::FortranBranchComponent>
  FortranNetwork;

typedef gridpack::mapper::FullMatrixMap<FortranNetwork> FortranFullMatrixMap;
typedef gridpack::mapper::BusVectorMap<FortranNetwork> FortranBusVectorMap;
typedef gridpack::mapper::GenMatrixMap<FortranNetwork> FortranGenMatrixMap;
typedef gridpack::mapper::GenVectorMap<FortranNetwork> FortranGenVectorMap;
typedef gridpack::math::Matrix  FortranMatrix;
typedef gridpack::math::Vector  FortranVector;

struct networkWrapper {
    boost::shared_ptr<FortranNetwork> network;
};

struct matrixWrapper {
  boost::shared_ptr<FortranMatrix> matrix;
};

struct vectorWrapper {
  boost::shared_ptr<FortranVector> vector;
};

/**
 * Create a new FullMatrixMap
 * @param map pointer to Fortran FullMatrixMap object
 * @param network pointer to Fortran network object
 */
extern "C" void p_full_matrix_map_create(FortranFullMatrixMap **map,
    networkWrapper *wnetwork)
{
  *map = new FortranFullMatrixMap(wnetwork->network);
}

/**
 * Destroy old mapper
 * @param map pointer to Fortran FullMatrixMap object
 */
extern "C" void p_full_matrix_map_destroy(FortranFullMatrixMap **map)
{
  delete (*map);
  *map = NULL;
}

/**
 * Create a matrix from the network
 * @param map pointer to Fortran FullMatrixMap object
 * @return pointer to matrix wrapper
 */
extern "C" void* p_full_matrix_map_map_to_matrix(FortranFullMatrixMap *map)
{
  FortranMatrix *matrix = map->intMapToMatrix();
  return matrix;
}

/**
 * Update a matrix from the network
 * @param map pointer to Fortran FullMatrixMap object
 * @param matrix pointer to matrix
 */
extern "C" void p_full_matrix_map_remap_to_matrix(FortranFullMatrixMap *map,
    FortranMatrix *matrix)
{
  map->mapToMatrix(*matrix);
}

/**
 * Overwrite selected elements of existing matrix.
 * @param map pointer to Fortran FullMatrixMap object
 * @param pointer to matrix wrapper
*/
extern "C" void p_full_matrix_map_overwrite_matrix(FortranFullMatrixMap *map,
    matrixWrapper *wmatrix)
{
  map->overwriteMatrix(wmatrix->matrix);
}

/**
 * Increment selected elements of existing matrix.
 * @param map pointer to Fortran FullMatrixMap object
 * @param pointer to matrix wrapper
*/
extern "C" void p_full_matrix_map_increment_matrix(FortranFullMatrixMap *map,
    matrixWrapper *wmatrix)
{
  map->incrementMatrix(wmatrix->matrix);
}

/**
 * Check to see if matrix looks well formed. This method runs through all
 * branches and verifies that the dimensions of the branch contributions match
 * the dimensions of the bus contributions at each end. If there is a
 * discrepancy, an error message is generated.
 * @param map pointer to Fortran FullMatrixMap object
 * @return true if no discrepancy found
 */
extern "C" bool p_full_matrix_map_check(FortranFullMatrixMap *map)
{
  return map->check();
}

/**
 * Create a new BusVectorMap
 * @param map pointer to Fortran BusVectorMap object
 * @param network pointer to Fortran network object
 */
extern "C" void p_bus_vector_map_create(FortranBusVectorMap **map,
    networkWrapper *wnetwork)
{
  *map = new FortranBusVectorMap(wnetwork->network);
}

/**
 * Destroy old mapper
 * @param map pointer to Fortran BusVectorMap object
 */
extern "C" void p_bus_vector_map_destroy(FortranBusVectorMap **map)
{
  delete (*map);
  *map = NULL;
}

/**
 * Create a vector from the network
 * @param map pointer to Fortran BusVectorMap object
 * @return pointer to vector wrapper
 */
extern "C" void* p_bus_vector_map_map_to_vector(FortranBusVectorMap *map)
{
  FortranVector *vector = map->intMapToVector();
//  printf("p[%d] Completed intMapToVector vector: %p map: %p\n",GA_Nodeid(),vector,map);
  return vector;
}

/**
 * Update a vector from the network
 * @param map pointer to Fortran BusVectorMap object
 * @param vector pointer to vector
 */
extern "C" void p_bus_vector_map_remap_to_vector(FortranBusVectorMap *map,
    FortranVector *vector)
{
  map->mapToVector(*vector);
}

/**
 * Push data from a vector to the network buses
 * @param map pointer to Fortran BusVectorMap object
 * @param pointer to vector wrapper
 */
extern "C" void p_bus_vector_map_map_to_bus(FortranBusVectorMap *map,
    FortranVector *vector)
{
//  printf("p[%d] Calling mapToBus vector vector: %p map: %p\n",GA_Nodeid(),vector,map);
  map->mapToBus(*vector);
}

/**
 * Create a new GenMatrixMap
 * @param map pointer to Fortran GenMatrixMap object
 * @param network pointer to Fortran network object
 */
extern "C" void p_gen_matrix_map_create(FortranGenMatrixMap **map,
    networkWrapper *wnetwork)
{
  *map = new FortranGenMatrixMap(wnetwork->network);
}

/**
 * Destroy old mapper
 * @param map pointer to Fortran GenMatrixMap object
 */
extern "C" void p_gen_matrix_map_destroy(FortranGenMatrixMap **map)
{
  delete (*map);
  *map = NULL;
}

/**
 * Create a matrix from the network
 * @param map pointer to Fortran GenMatrixMap object
 * @return pointer to matrix wrapper
 */
extern "C" void* p_gen_matrix_map_map_to_matrix(FortranGenMatrixMap *map)
{
  matrixWrapper *matrix = new matrixWrapper;
  matrix->matrix = map->mapToMatrix();
  return matrix;
}

/**
 * Update a matrix from the network
 * @param map pointer to Fortran GenMatrixMap object
 * @param pointer to matrix wrapper
 */
extern "C" void p_gen_matrix_map_remap_to_matrix(FortranGenMatrixMap *map,
    matrixWrapper *wmatrix)
{
  map->mapToMatrix(wmatrix->matrix);
}

/**
 * Create a new GenVectorMap
 * @param map pointer to Fortran GenVectorMap object
 * @param network pointer to Fortran network object
 */
extern "C" void p_gen_vector_map_create(FortranGenVectorMap **map,
    networkWrapper *wnetwork)
{
  *map = new FortranGenVectorMap(wnetwork->network);
}

/**
 * Destroy old mapper
 * @param map pointer to Fortran GenVectorMap object
 */
extern "C" void p_gen_vector_map_destroy(FortranGenVectorMap **map)
{
  delete (*map);
  *map = NULL;
}

/**
 * Create a vector from the network
 * @param map pointer to Fortran GenVectorMap object
 * @return pointer to vector wrapper
 */
extern "C" void* p_gen_vector_map_map_to_vector(FortranGenVectorMap *map)
{
  vectorWrapper *vector = new vectorWrapper;
  vector->vector = map->mapToVector();
  return vector;
}

/**
 * Update a vector from the network
 * @param map pointer to Fortran GenVectorMap object
 * @param pointer to vector wrapper
 */
extern "C" void p_gen_vector_map_remap_to_vector(FortranGenVectorMap *map,
    vectorWrapper *wvector)
{
  map->mapToVector(wvector->vector);
}

/**
 * Push data from a vector to the network buses
 * @param map pointer to Fortran GenVectorMap object
 * @param pointer to vector wrapper
 */
extern "C" void p_gen_vector_map_map_to_bus(FortranGenVectorMap *map,
    FortranVector *vector)
{
  printf("Function mapToBus function not implemented yet in GenVectorMap\n");
//  map->mapToBus(*vector);
}
