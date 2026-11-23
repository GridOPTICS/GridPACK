// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hat_matrix_diagonal.hpp
 * @brief  Diagonal of the weighted least squares hat matrix
 */
// -------------------------------------------------------------

#ifndef _hat_matrix_diagonal_hpp_
#define _hat_matrix_diagonal_hpp_

#include <string>
#include <vector>
#include "gridpack/math/matrix.hpp"

namespace gridpack {
namespace math {

/// Ways to compute the diagonal
enum HatMatrixDiagonalMethod {
  SelectedInversion,  ///< sparse factorization and selected inverse (needs KLU)
  BlockedSolves       ///< one solve per row, in blocks of right-hand sides
};

/// What was done and how long it took (times are those of process 0)
struct HatMatrixDiagonalInfo {
  std::string method;   ///< method that produced the result
  std::string note;     ///< why the requested method was not used, if so
  int rows;             ///< rows of H
  int cols;             ///< columns of H
  double gatherTime;    ///< collecting H on process 0
  double factorTime;    ///< forming and factorizing H^T W H
  double diagonalTime;  ///< computing the diagonal from the factors
};

/// Diagonal of H (H^T W H)^-1 H^T W for the rows owned by this process
/**
 * H is collected on process 0 and factorized there, so it must fit in the
 * memory of one process. Only the real parts of H are used. Collective.
 *
 * @param H matrix with one row per equation and full column rank
 * @param localWeights diagonal of W for the rows owned by this process
 * @param localDiagonal result for the rows owned by this process
 * @param info method used and timings
 * @param method requested method
 */
extern void hatMatrixDiagonal(const Matrix& H,
                              const std::vector<double>& localWeights,
                              std::vector<double>& localDiagonal,
                              HatMatrixDiagonalInfo& info,
                              HatMatrixDiagonalMethod method = SelectedInversion);

} // namespace math
} // namespace gridpack

#endif
