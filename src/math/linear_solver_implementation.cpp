// -------------------------------------------------------------
/**
 * @file   linear_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2013-06-14 12:05:23 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created June 14, 2013 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#include "linear_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearSolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// LinearSolverImplementation:: constructors / destructor
// -------------------------------------------------------------
LinearSolverImplementation::LinearSolverImplementation(const Matrix& A)
  : parallel::Distributed(A.communicator()),
    utility::Uncopyable(),
    p_A(A.clone())
{
  
}

LinearSolverImplementation::~LinearSolverImplementation(void)
{
  // empty
}

// -------------------------------------------------------------
// LinearSolverImplementation
// -------------------------------------------------------------



} // namespace math
} // namespace gridpack
