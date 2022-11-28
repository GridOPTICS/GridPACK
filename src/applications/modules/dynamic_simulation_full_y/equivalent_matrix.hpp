/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   equivalent_matrix.hpp
 * @author Bruce Palmer
 * @Last modified:   March 24, 2022
 * 
 * @brief
 * This routine is designed to construct a reduced matrix from the values of a
 * collection of solution vectors. Each vector is added to the matrix one column
 * (or row) at a time. The entire matrix can then be accessed from any process.
 * 
 * 
 */

#ifndef equivalent_matrix_h_
#define equivalent_matrix_h_

#include <vector>
#include <map>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/math/math.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "dsf_factory.hpp"

namespace gridpack {
namespace dynamic_simulation {
class EquivalentMatrix
{
  public:
    /**
     * Basic constructor
     * @param[in] config pointer to configuration object
     * @param[in] network pointer to dynamic simulation network. Network should
     *            already be distributed and initialized.
     * @param[in] flag: 1 for generating negative sequence equivalent matrix
                  flag: 2 for generating zero sequence equivalent matrix
     */
    EquivalentMatrix(gridpack::utility::Configuration *config,
        boost::shared_ptr<DSFullNetwork> &network, int flag = 0);

    /**
     * Basic destructor
     */
    virtual ~EquivalentMatrix();

    /**
     * Get a complete copy of matrix and copy it into a vector
     * @param[out] matrix vector data structure containing full matrix
     * @param[out] n dimension of matrix
     */
    void getMatrix(std::vector<ComplexType> &matrix, int &n);

    /**
     * Print equivalent matrix to standard out
     */
    void printMatrix();

  private:

    /**
     * Add column to matrix
     * @param[in] vector add distributed vector to matrix
     * @param[in] bus_col ID of bus corresponding to column location
     */
    void addColumn(gridpack::math::Vector *vector, int bus_col);

    /**
     * Add row to matrix
     * @param[in] vector add distributed vector to matrix
     * @param[in] bus_row ID of bus corresponding to row location
     */
    void addRow(gridpack::math::Vector *vector, int bus_row);

    int p_ga;  // global array handle used to store matrix
    int p_dim; // Size of matrix
    std::map<int,int> p_map; // map bus ID to matrix index
    boost::shared_ptr<DSFullNetwork> p_network; // pointer to network
    bool p_initialized; // object has been propertly initialized with valid data
    boost::shared_ptr < gridpack::mapper::BusVectorMap<DSFullNetwork> > p_ivecMap;
};
}  // dynamic_simulation
}  // gridpack
#endif
