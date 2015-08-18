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
 * @date   2015-08-18 13:40:52 d3g096
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
  LinearSolverImplementation(MatrixType& A)
    : parallel::Distributed(A.communicator()),
      utility::Configurable("LinearSolver"),
      utility::Uncopyable(),
      p_matrix(A),
      p_solutionTolerance(1.0e-06),
      p_relativeTolerance(p_solutionTolerance),
      p_maxIterations(100),
      p_doSerial(false),
      p_constSerialMatrix(),
      p_guessZero(false),
      p_serialSolution()
  {
  }

  /// Destructor
  ~LinearSolverImplementation(void)
  {
    // empty
  }

protected:

  /// A reference to coefficient matrix
  MatrixType& p_matrix;

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

  /// Use a serial solver, even in parallel environment
  bool p_doSerial;

  /// Assume the coefficient matrix is constant (only if p_doSerial)
  /**
   * When a serial solver is used in a parallel environment
   * (::p_doSerial is true), the coefficient matrix is collected from
   * all processes for each solve.  This is unnecessary if the
   * coefficient matrix is constant. If this flag is true, the
   * coefficient matrix will be collected only once (and stored as
   * ::p_serialMatrix). 
   * 
   */
  bool p_constSerialMatrix;

  /// A place to keep the (constant) "serial" matrix if called for
  mutable boost::scoped_ptr<MatrixType> p_serialMatrix;

  /// A place to keep the collected "serial" RHS vector, if called for
  mutable boost::scoped_ptr<VectorType> p_serialRHS;

  /// Assume the initial guess is zero
  bool p_guessZero;

  /// A place to keep the "serial" solution vector, if called for
  /**
   * When a serial solver is used in a parallel environment
   * (::p_doSerial is true), the RHS and solution (guess) vectors are
   * collected from all processors. If the solution guess is zero
   * (::p_guessZero is true), then collecting the solution vector is
   * unnecessary.  This provides a place to put a serial solution
   * vector.
   * 
   */
  mutable boost::scoped_ptr<VectorType> p_serialSolution;

  /// A buffer to use for value transfer
  mutable std::vector<TheType> p_valueBuffer;

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props)
  {
    if (props) {
      p_solutionTolerance = props->get("SolutionTolerance", p_solutionTolerance);
      p_relativeTolerance = props->get("RelativeTolerance", p_relativeTolerance);
      p_maxIterations = props->get("MaxIterations", p_maxIterations);

      p_doSerial = props->get("ForceSerial", p_doSerial);
      p_constSerialMatrix = props->get("SerialMatrixConstant", p_constSerialMatrix);

      // SerialOnly has no effect unless parallel
      p_doSerial = (p_doSerial && (this->processor_size() > 1));

      p_guessZero = props->get("InitialGuessZero", p_guessZero);
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

  /// Solve the specified system w/ RHS and estimate (implementation)  
  virtual void p_solveImpl(MatrixType& A, const VectorType& b, VectorType& x) const = 0;

  /// Solve the system again w/ RHS and estimate (implementation)
  virtual void p_resolveImpl(const VectorType& b, VectorType& x) const = 0;

  /// Gather the RHS and initial estimate vectors
  void p_serialSolvePrep(const VectorType& b, VectorType& x) const
  {
    // make a buffer for value transfer between vectors
    if (p_valueBuffer.empty()) {
      p_valueBuffer.resize(b.size());
    }
    
    // always collect the RHS vector 
    if (!p_serialRHS) {
      p_serialRHS.reset(b.localClone());
    } else {
      b.getAllElements(&p_valueBuffer[0]);
      p_serialRHS->setElementRange(0, b.size() - 1, &p_valueBuffer[0]);
    }

    // Collect the initial guess, if it's not zero
    if (!p_serialSolution) {
      if (!p_guessZero) {
        p_serialSolution.reset(x.localClone());
      } else {
        p_serialSolution.reset(p_serialRHS->clone());
      }
    } else {
      if (!p_guessZero) {
        x.getAllElements(&p_valueBuffer[0]);
        p_serialSolution->setElementRange(0, x.size() - 1, &p_valueBuffer[0]);
      }
    }
  }

  /// After serial solve, distribute solution
  void p_serialSolutionDist(VectorType& x) const
  {
    BOOST_ASSERT(!p_valueBuffer.empty());
    IdxType lo, hi;
    x.localIndexRange(lo, hi);
    p_serialSolution->getElementRange(lo, hi, &p_valueBuffer[0]);
    x.setElementRange(lo, hi, &p_valueBuffer[0]);
    x.ready();
  }

  /// Solve w/ the specified RHS and estimate (result in x)
  void p_solve(const VectorType& b, VectorType& x) const
  {
    if (p_doSerial) {

      // collect the coefficient matrix to the local processor, if not
      // already here or it's not constant

      if (!p_serialMatrix ) {
          p_serialMatrix.reset(p_matrix.localClone());
      } else if (!p_constSerialMatrix) {
          p_serialMatrix.reset(p_matrix.localClone());
      }
      
      this->p_serialSolvePrep(b, x);

      // solve the system serially

      this->p_solveImpl(*p_serialMatrix, *p_serialRHS, *p_serialSolution);

      // each process puts its share of the solution into the parallel
      // solution vector

      this->p_serialSolutionDist(x);


    } else {

      this->p_solveImpl(p_matrix, b, x);

    }
  }

  /// Solve again w/ the specified RHS, put result in specified vector (specialized)
  void p_resolve(const VectorType& b, VectorType& x) const
  {
    if (p_doSerial) {

      this->p_serialSolvePrep(b, x);

      // solve the system serially

      this->p_resolveImpl(*p_serialRHS, *p_serialSolution);

      // each process puts its share of the solution into the parallel
      // solution vector

      this->p_serialSolutionDist(x);

    } else {

      this->p_resolveImpl(b, x);

    }
  }

  /// Solve multiple systems w/ each column of the Matrix a single RHS
  MatrixType *p_solve(const MatrixType& B) const
  {
    VectorType b(B.communicator(), B.localRows());
    VectorType X(B.communicator(), B.localRows());
    MatrixType *result(new MatrixType(B.communicator(), B.localRows(), B.localCols(), Dense));

    int ilo, ihi;
    X.localIndexRange(ilo, ihi);
    int nloc(X.localSize());
    std::vector<IdxType> iidx;
    iidx.reserve(nloc);
    for (IdxType i = ilo; i < ihi; ++i) { iidx.push_back(i); }
    std::vector<IdxType> jidx(nloc);
    std::vector<TheType> locX(nloc);

    for (int j = 0; j < B.cols(); ++j) {
      column(B, j, b);
      X.zero();
      X.ready();
      if (j == 0) {
        this->solve(b, X);
      } else {
        this->resolve(b, X);
      }
      std::fill(jidx.begin(), jidx.end(), j);
      X.getElements(nloc, &iidx[0], &locX[0]);
      result->setElements(nloc, &iidx[0], &jidx[0], &locX[0]);
    }
  
    result->ready();
    return result;
  }

};

} // namespace math
} // namespace gridpack



#endif
