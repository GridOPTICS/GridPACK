/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/* *****************************************************************************
 * gen_matrix_map.hpp
 * gridpack
 * Bruce Palmer
 * July 1, 2014
 *
 * Documentation status
 * Revision status
 * File contents
 **************************************************************************** */

#ifndef GENMATRIXMAP_HPP_
#define GENMATRIXMAP_HPP_

#define NZ_PER_ROW

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
#include <gridpack/mapper/base_gen_matrix_map.hpp>

//#define DBG_CHECK

namespace gridpack {
namespace mapper {

template <class _network>
class GenMatrixMap {
  public: 
/**
 * Initialize mapper for the given network and the current mode. Create global
 * arrays that contain offsets that will be used to create matrix from the
 * network component objects
 * @param network network that will generate matrix
 */
GenMatrixMap(boost::shared_ptr<_network> network)
{
  p_base_mapper.reset(new BaseGenMatrixMap<_network>(network));
  p_base_mapper->setMode(gridpack::component::BaseGenMatVecInterface::STANDARD);
}

~GenMatrixMap()
{
}

/**
 * Generate matrix from current component state on network
 * @return return a pointer to new matrix
 */
boost::shared_ptr<gridpack::math::Matrix> mapToMatrix(void)
{
  return p_base_mapper->baseMapToMatrix();
}

/**
 * Reset existing matrix from current component state on network
 * @param matrix existing matrix (should be generated from same mapper)
 */
void mapToMatrix(gridpack::math::Matrix &matrix)
{
  p_base_mapper->baseMapToMatrix(matrix);
}

/**
 * Reset existing matrix from current component state on network
 * @param matrix existing matrix (should be generated from same mapper)
 */
void mapToMatrix(boost::shared_ptr<gridpack::math::Matrix> &matrix)
{
  p_base_mapper->baseMapToMatrix(matrix);
}

/**
 * Overwrite elements of existing matrix. This can be used to overwrite selected
 * elements of a matrix
 * @param matrix existing matrix (should be generated from same mapper)
 */
void overwriteMatrix(gridpack::math::Matrix &matrix)
{
  p_base_mapper->baseOverwriteMatrix(matrix);
}

/**
 * Overwrite elements of existing matrix. This can be used to overwrite selected
 * elements of a matrix
 * @param matrix existing matrix (should be generated from same mapper)
 */
void overwriteMatrix(boost::shared_ptr<gridpack::math::Matrix> &matrix)
{
  p_base_mapper->baseOverwriteMatrix(matrix);
}

/**
 * Increment elements of existing matrix. This can be used to increment selected
 * elements of a matrix
 * @param matrix existing matrix (should be generated from same mapper)
 */
void incrementMatrix(gridpack::math::Matrix &matrix)
{
  p_base_mapper->baseIncrementMatrix(matrix);
}

/**
 * Increment elements of existing matrix. This can be used to increment selected
 * elements of a matrix
 * @param matrix existing matrix (should be generated from same mapper)
 */
void incrementMatrix(boost::shared_ptr<gridpack::math::Matrix> &matrix)
{
  p_base_mapper->baseIncrementMatrix(matrix);
}

private:

boost::shared_ptr<gridpack::mapper::BaseGenMatrixMap<_network> > p_base_mapper;


};

} /* namespace mapper */
} /* namespace gridpack */

#endif //GENMATRIXMAP_HPP_
