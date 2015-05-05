// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   linear_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2015-03-05 12:49:11 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  1, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _linear_solver_implementation_hpp_
#define _linear_solver_implementation_hpp_

#include <boost/scoped_ptr.hpp>
#include <gridpack/math/linear_solver_interface.hpp>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/configuration/configurable.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearSolverImplementation
// -------------------------------------------------------------
template <typename T, typename I>
class LinearSolverImplementation 
  : public parallel::Distributed,
    public utility::Configurable,
    private utility::Uncopyable,
    public BaseLinearSolverInterface<T, I>
{
public:

  typedef T TheType;
  typedef I IdxType;
  typedef typename BaseLinearSolverInterface<T, I>::MatrixType MatrixType;
  typedef typename BaseLinearSolverInterface<T, I>::VectorType VectorType;
  
  /// Default constructor.
  LinearSolverImplementation(const parallel::Communicator& comm)
    : parallel::Distributed(comm),
      utility::Configurable("LinearSolver"),
      utility::Uncopyable(),
      p_solutionTolerance(1.0e-06),
      p_relativeTolerance(p_solutionTolerance),
      p_maxIterations(100)
  {
  }

  /// Destructor
  ~LinearSolverImplementation(void)
  {
    // empty
  }

protected:

  /// The solution residual norm tolerance
  /**
   * This is the absolute solution tolerance. The linear system is
   * considered solved when the solution residual norm is below this value.
   * This is just handed to the underlying math library, so 
   * 
   */
  double p_solutionTolerance;

  /// The relative solution tolerance.
  /**
   * The linear system is considered solved when the solution residual
   * norm is @em reduced by the specified amount.
   * 
   */
  double p_relativeTolerance;

  /// The maximum number of iterations to perform
  int p_maxIterations;

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props)
  {
    if (props) {
      p_solutionTolerance = props->get("SolutionTolerance", p_solutionTolerance);
      p_relativeTolerance = props->get("RelativeTolerance", p_solutionTolerance);
      p_maxIterations = props->get("MaxIterations", p_maxIterations);
    }
  }

  /// Get the solution tolerance (specialized)
  double p_tolerance(void) const
  {
    return p_solutionTolerance;
  }

  /// Set the solver tolerance (specialized)
  void p_tolerance(const double& tol)
  {
    p_solutionTolerance = tol;
  }

  /// Get the maximum iterations (specialized)
  int p_maximumIterations(void) const
  {
    return p_maxIterations;
  }

  /// Set the maximum solution iterations
  void p_maximumIterations(const int& n)
  {
    p_maxIterations = n;
  }

  /// Solve multiple systems w/ each column of the Matrix a single RHS
  MatrixType *p_solve(const MatrixType& B) const
  {
    VectorType b(B.communicator(), B.localRows());
    VectorType X(B.communicator(), B.localRows());
    MatrixType *result(new MatrixType(B.communicator(), B.localRows(), B.localCols(), Dense));

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
      for (IdxType i = ilo; i < ihi; ++i) {
        TheType v;
        X.getElement(i, v);
        result->setElement(i, j, v);
      }
    }
  
    result->ready();
    return result;
  }

};

} // namespace math
} // namespace gridpack



#endif
