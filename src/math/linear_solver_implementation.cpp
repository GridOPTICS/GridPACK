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
 * @date   2014-10-22 09:12:58 d3g096
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
LinearSolverImplementation::LinearSolverImplementation(const parallel::Communicator& comm)
  : parallel::Distributed(comm),
    utility::Configurable("LinearSolver"),
    utility::Uncopyable(),
    p_solutionTolerance(1.0e-06),
    p_relativeTolerance(p_solutionTolerance),
    p_maxIterations(100)
{
  
}

LinearSolverImplementation::~LinearSolverImplementation(void)
{
  // empty
}

// -------------------------------------------------------------
// LinearSolverImplementation::p_solve
// -------------------------------------------------------------
Matrix *
LinearSolverImplementation::p_solve(const Matrix& B) const
{
  Vector b(B.communicator(), B.localRows());
  Vector X(B.communicator(), B.localRows());
  Matrix *result(new Matrix(B.communicator(), B.localRows(), B.localCols(), Matrix::Dense));

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

// -------------------------------------------------------------
// LinearSolverImplementation::p_configure
// -------------------------------------------------------------
void
LinearSolverImplementation::p_configure(utility::Configuration::CursorPtr props)
{
  if (props) {
    p_solutionTolerance = props->get("SolutionTolerance", p_solutionTolerance);
    p_relativeTolerance = props->get("RelativeTolerance", p_solutionTolerance);
    p_maxIterations = props->get("MaxIterations", p_maxIterations);
  }
}


// -------------------------------------------------------------
// LinearSolverImplementation::p_tolerance
// -------------------------------------------------------------
double
LinearSolverImplementation::p_tolerance(void) const
{
  return p_solutionTolerance;
}

void
LinearSolverImplementation::p_tolerance(const double& tol)
{
  p_solutionTolerance = tol;
}

// -------------------------------------------------------------
// LinearSolverImplementation::p_maximumIterations
// -------------------------------------------------------------
int
LinearSolverImplementation::p_maximumIterations(void) const
{
  return p_maxIterations;
}

void
LinearSolverImplementation::p_maximumIterations(const int& n) 
{
  p_maxIterations = n;
}



} // namespace math
} // namespace gridpack
