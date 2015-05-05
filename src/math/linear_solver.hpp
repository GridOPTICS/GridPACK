// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   linear_solver.hpp
 * @author William A. Perkins
 * @date   2015-03-05 13:09:23 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _linear_solver_hpp_
#define _linear_solver_hpp_

#include <boost/scoped_ptr.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/math/linear_solver_implementation.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearSolver
// -------------------------------------------------------------
/// A solver for system(s) of linear equations in parallel
/**
 * A LinearSolver is used to solve a system of linear equations:
 * \f[
 *   \left[ \mathbf{A} \right] \mathbf{x} ~ = ~ \mathbf{b}
 * \f]
 * where \f$\left[ \mathbf{A} \right]\f$ is the coefficient matrix,
 * and \f$\mathbf{x}\f$ and \f$\mathbf{b}\f$ are vectors.
 *
 * A LinearSolver is instantiated using an existing coefficient
 * Matrix.  It must be intantiated simutaneously on all processes
 * involved in the \ref parallel::Communicator "communicator" used by
 * the coefficent Matrix. It similarly must be destroyed
 * simultaneously on all processes.
 *
 * A particular system is solved using the solve() method, which
 * requires the \f$\mathbf{b}\f$ Vector and a Vector that serves as
 * the intial estimate of \f$\mathbf{x}\f$. Multiple systems can be
 * solved using the same coefficient matrix by repeatedly calling
 * solve() with different \f$\mathbf{b}\f$ and \f$\mathbf{x}\f$
 * vectors.  If the coefficient matrix changes in any way (even just
 * the coefficients) between system solutions, the setMatrix() method
 * must be called before solve();
 * 
 * This class encapuslates the linear system solver of the underlying
 * math library. The Pimpl idiom is used for \ref
 * LinearSolverImplementation "implementation", so user code is
 * completely independent of the underlying library. This class simply
 * provides an interface to a specific \ref LinearSolverImplementation
 * "implementation".
 * 
 */
template <typename T, typename I = int>
class LinearSolverT 
  : public parallel::WrappedDistributed,
    public utility::WrappedConfigurable,
    private utility::Uncopyable, 
    public BaseLinearSolverInterface<T, I>
{
public:
  
  typedef typename BaseLinearSolverInterface<T, I>::MatrixType MatrixType;
  typedef typename BaseLinearSolverInterface<T, I>::VectorType VectorType;

  /// Default constructor.
  /** 
   * @e Collective
   *
   * This will create a LinearSolver instance on the same \ref
   * parallel::Communicator "communicator" used by the coefficient
   * Matrix @c A.  All processes in the \ref parallel::Communicator
   * "communicator" must call this simultaneously.
   *
   * The Matrix @c A must exist for the life of this instance. 
   * 
   * @param A existing, filled coefficient matrix
   * 
   * @return new LinearSolver instance
   */
  LinearSolverT(MatrixType& A);
  
  /// Destructor
  /** 
   * @e Collective
   * 
   * This destroys a LinearSolver. It must be called simutaneously by
   * all processes involved in calling the constructor.
   * 
   */
  ~LinearSolverT(void)
  {}

protected:

  /// Where the work really happens
  /**
   * The Pimpl idiom is used for \ref LinearSolverImplementation
   * "implementation", so user code is completely independent of the
   * underlying library. This class simply provides an interface on a
   * specific \ref LinearSolverImplementation "implementation".
   */
  boost::scoped_ptr< LinearSolverImplementation<T, I> > p_solver;

  /// Get the solution tolerance (specialized)
  /** 
   * 
   * 
   * 
   * @return current solution tolerance
   */
  double p_tolerance(void) const
  {
    return p_solver->tolerance();
  }

  /// Set the solver tolerance (specialized)
  /** 
   * 
   * 
   * @param tol new solution tolerance
   */
  void p_tolerance(const double& tol)
  {
    p_solver->tolerance(tol);
  }

  /// Get the maximum iterations (specialized)
  /** 
   * 
   * 
   * 
   * @return current maximum number of solution iterations
   */
  int p_maximumIterations(void) const
  {
    return p_solver->maximumIterations();
  }


  /// Set the maximum solution iterations (specialized)
  /** 
   * 
   * 
   * @param n new maximum number of iterations
   */
  void p_maximumIterations(const int& n)
  {
    p_solver->maximumIterations(n);
  }

  /// Solve w/ the specified RHS, put result in specified vector
  /** 
   * @e Collective.
   *
   * Solve a linear system of equations. When called, Vector @c x
   * should contain the initial solution estimate.  The final solution
   * is returned in @c x.  
   *
   * The \ref parallel::Communicator "communicator" @c x and @c b must
   * be the same and match that of the coefficient Matrix used for
   * construction. The length of both @c x and @c b must the number of
   * columns in the coeffienct Matrix used for constructor or passed
   * the last call to setMatrix().  If these conditions are not met,
   * an \ref Exception "exception" is thrown.
   *  
   * @param b Vector containing right hand side of linear system
   * @param x solution Vector
   */
  void p_solve(const VectorType& b, VectorType& x) const
  {
    p_solver->solve(b, x);
  }

  /// Solve multiple systems w/ each column of the Matrix a single RHS
  /** 
   * 
   * 
   * @param B RHS matrix -- each column is used as a RHS Vector
   * 
   * @return @e dense solution Matrix -- each column is the solution for the corresponding column in @c B
   */
  MatrixType *p_solve(const MatrixType& B) const
  {
    return p_solver->solve(B);
  }


};

typedef LinearSolverT<ComplexType> ComplexLinearSolver;
typedef LinearSolverT<RealType> RealLinearSolver;

typedef ComplexLinearSolver LinearSolver;

} // namespace math
} // namespace gridpack

#endif
