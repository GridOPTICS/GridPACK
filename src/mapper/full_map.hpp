/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/* *****************************************************************************
 * full_map.hpp
 * gridpack
 * kglass, bjpalmer
 * Jul 22, 2013
 *
 * Documentation status
 * Revision status
 * File contents
 **************************************************************************** */

#ifndef FULLMATRIXMAP_HPP_
#define FULLMATRIXMAP_HPP_

//#define NZ_PER_ROW

#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
#include "gridpack/parallel/parallel.hpp"
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/component/base_component.hpp>
#include <gridpack/component/base_matrix_ifc.hpp>
#include <gridpack/component/matvec_ifc.hpp>
#include <gridpack/network/base_network.hpp>
#include <gridpack/math/matrix.hpp>
#include <gridpack/utilities/exception.hpp>
#include <gridpack/mapper/base_matrix_map.hpp>

#define DBG_CHECK

namespace gridpack {
namespace mapper {

template <class _network>
class FullMatrixMap : BaseMatrixMap<_network> {
  public: 
/**
 * Initialize mapper for the given network and the current mode. Create global
 * arrays that contain offsets that will be used to create matrix from the
 * network component objects
 * @param network network that will generate matrix
 */
FullMatrixMap(boost::shared_ptr<_network> network) 
  : BaseMatrixMap<_network>::BaseMatrixMap(network)
{
  this->setMode(gridpack::component::BaseMatrixInterface::STANDARD);
}

~FullMatrixMap()
{
}

/**
 * Generate matrix from current component state on network
 * @param isDense set to true if creating a dense matrix
 * @return return a pointer to new matrix
 */
boost::shared_ptr<gridpack::math::Matrix> mapToMatrix(bool isDense = false)
{
  return this->baseMapToMatrix(isDense);
}

/**
 * Generate real matrix from current component state on network
 * @param isDense set to true if creating a dense matrix
 * @return return a pointer to new matrix
 */
boost::shared_ptr<gridpack::math::RealMatrix> mapToRealMatrix(bool isDense = false)
{
  return this->baseMapToRealMatrix(isDense);
}

/**
 * Reset existing matrix from current component state on network
 * @param matrix existing matrix (should be generated from same mapper)
 */
void mapToMatrix(gridpack::math::Matrix &matrix)
{
  this->baseMapToMatrix(matrix);
}

/**
 * Reset existing matrix from current component state on network
 * @param matrix existing matrix (should be generated from same mapper)
 */
void mapToRealMatrix(gridpack::math::RealMatrix &matrix)
{
  this->baseMapToRealMatrix(matrix);
}

/**
 * Reset existing matrix from current component state on network
 * @param matrix existing matrix (should be generated from same mapper)
 */
void mapToMatrix(boost::shared_ptr<gridpack::math::Matrix> &matrix)
{
  this->baseMapToMatrix(*matrix);
}

/**
 * Reset existing matrix from current component state on network
 * @param matrix existing matrix (should be generated from same mapper)
 */
void mapToRealMatrix(boost::shared_ptr<gridpack::math::RealMatrix> &matrix)
{
  this->baseMapToRealMatrix(*matrix);
}

/**
 * Overwrite elements of existing matrix. This can be used to overwrite selected
 * elements of a matrix
 * @param matrix existing matrix (should be generated from same mapper)
 */
void overwriteMatrix(gridpack::math::Matrix &matrix)
{
  this->baseOverwriteMatrix(matrix);
}

/**
 * Overwrite elements of existing matrix. This can be used to overwrite selected
 * elements of a matrix
 * @param matrix existing matrix (should be generated from same mapper)
 */
void overwriteMatrix(boost::shared_ptr<gridpack::math::Matrix> &matrix)
{
  this->baseOverwriteMatrix(*matrix);
}

/**
 * Increment elements of existing matrix. This can be used to increment selected
 * elements of a matrix
 * @param matrix existing matrix (should be generated from same mapper)
 */
void incrementMatrix(gridpack::math::Matrix &matrix)
{
  this->baseIncrementMatrix(matrix);
}

/**
 * Increment elements of existing matrix. This can be used to increment selected
 * elements of a matrix
 * @param matrix existing matrix (should be generated from same mapper)
 */
void incrementMatrix(boost::shared_ptr<gridpack::math::Matrix> &matrix)
{
  this->baseIncrementMatrix(*matrix);
}

/**
 * Check to see if matrix looks well formed. This method runs through all
 * branches and verifies that the dimensions of the branch contributions match
 * the dimensions of the bus contributions at each end. If there is a
 * discrepancy, then an error message is generated.
 */

bool check(void)
{
  this->baseCheck();
}

private:

};

} /* namespace mapper */
} /* namespace gridpack */

#endif //FULLMATRIXMAP_HPP_
