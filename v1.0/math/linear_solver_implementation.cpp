// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   linear_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2013-10-11 11:17:14 d3g096
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

#include <boost/format.hpp>

#include "linear_solver_implementation.hpp"
#include "gridpack/utilities/exception.hpp"

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
    utility::Configurable("LinearSolver"),
    utility::Uncopyable(),
    p_A(A.clone())
{
  
}

LinearSolverImplementation::~LinearSolverImplementation(void)
{
  // empty
}

// -------------------------------------------------------------
// LinearSolverImplementation::solve
// -------------------------------------------------------------
Matrix *
LinearSolverImplementation::solve(const Matrix& B) const
{
  if (p_A->rows() != B.rows()) {
    std::string msg =
      boost::str(boost::format("LinearSolver::solve(Matrix): matrix size mismatch (%dx%d) and (%dx%d)") %
                 p_A->rows() % p_A->cols() % B.rows() % B.cols());
    throw Exception(msg);
  }

  Vector b(B.communicator(), B.localRows());
  Vector X(B.communicator(), B.localRows());
  Matrix *result(new Matrix(B.communicator(), B.localRows(), B.cols(), Matrix::Dense));

  int ilo, ihi;
  X.localIndexRange(ilo, ihi);
  // std::vector<ComplexType> locX(X.localSize());

  for (int j = 0; j < B.cols(); ++j) {
    column(B, j, b);
    X.zero();
    X.ready();
    this->solve(b, X);
    // std::cout << X.processor_rank() << ": X: " << ilo << "-" << ihi << std::endl;
    // X.print();
    // X.getElementRange(ilo, ihi, &locX[0]);
    for (int i = ilo; i < ihi; ++i) {
      ComplexType v;
      X.getElement(i, v);
      result->setElement(i, j, v);
    }
  }
  
  result->ready();
  return result;
}



} // namespace math
} // namespace gridpack
