// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   linear_matrix_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2013-11-07 12:42:58 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "linear_matrix_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearMatrixSolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// LinearMatrixSolverImplementation:: constructors / destructor
// -------------------------------------------------------------
LinearMatrixSolverImplementation::LinearMatrixSolverImplementation(const Matrix& A)
  : parallel::Distributed(A.communicator()),
    utility::Configurable(),
    utility::Uncopyable(),
    p_A(A.clone())
{
  configurationKey("LinearMatrixSolver");
}

LinearMatrixSolverImplementation::~LinearMatrixSolverImplementation(void)
{
}

} // namespace math
} // namespace gridpack
