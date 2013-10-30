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
 * @date   2013-10-11 10:35:20 d3g096
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
class LinearSolver 
  : public parallel::WrappedDistributed,
    private utility::Uncopyable
{
public:
  
  /// Default constructor.
  /** 
   * @e Collective
   *
   * This will create a LinearSolver instance on the same \ref
   * parallel::Communicator "communicator" used by the coefficient
   * Matrix @c A.  All processes in the \ref parallel::Communicator
   * "communicator" must call this simultaneously.
   * 
   * A copy is made of @c A so it's contents are in no way affected by
   * this call.
   * 
   * @param A existing, filled coefficient matrix
   * 
   * @return new LinearSolver instance
   */
  LinearSolver(const Matrix& A);
  
  /// Destructor
  /** 
   * @e Collective
   * 
   * This destroys a LinearSolver. It must be called simutaneously by
   * all processes involved in calling the constructor.
   * 
   */
  ~LinearSolver(void);

  /// Set the configuration key for this instance
  void configurationKey(const std::string& s)
  {
    p_solver->configurationKey(s);
  }

  /// Configure and do whatever is necessary to make this instance ready
  void configure(utility::Configuration::Cursor *props)
  {
    p_solver->configure(props);
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
  void solve(const Vector& b, Vector& x) const
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
  Matrix *solve(const Matrix& B) const
  {
    return p_solver->solve(B);
  }

  /// Use different coefficient matrix (i.e. new values)
  /** 
   * @e Collective
   * 
   * This should be called if the coefficient matrix changes in any
   * way between calls to solve(). Typically, this would just involve
   * changes to the Matrix used in \ref LinearSolver() "construction",
   * but a completely different Matrix instance can be used here. The
   * \ref parallel::Communicator "communicator" used by @c A should be
   * the same at \ref LinearSolver() "construction".
   *
   * A copy is made of @c A so it's contents are in no way affected by
   * this call.
   * 
   * @param A existing, filled coefficient matrix
   */
  void setMatrix(const Matrix& A)
  {
    p_solver->setMatrix(A);
  }
  

  //! @cond DEVDOC

  /// Allow visits by implemetation visitor
  /** 
   * This simply passes an ImplementationVisitor on to the \ref
   * LinearSolverImplementation "implementation".
   * 
   * @param visitor 
   */
  void accept(ImplementationVisitor& visitor)
  {
    p_solver->accept(visitor);
  }

  /// Allow visits by implemetation visitor
  /** 
   * This simply passes a ConstImplementationVisitor on to the \ref
   * LinearSolverImplementation "implementation".
   * 
   * @param visitor 
   */
  void accept(ConstImplementationVisitor& visitor) const
  {
    p_solver->accept(visitor);
  }

  //! @endcond

protected:

  /// Where the work really happens
  /**
   * The Pimpl idiom is used for \ref LinearSolverImplementation
   * "implementation", so user code is completely independent of the
   * underlying library. This class simply provides an interface on a
   * specific \ref LinearSolverImplementation "implementation".
   */
  boost::scoped_ptr<LinearSolverImplementation> p_solver;

};

} // namespace math
} // namespace gridpack

#endif
