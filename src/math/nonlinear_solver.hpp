// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver.hpp
 * @author William A. Perkins
 * @date   2015-03-25 15:54:13 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _nonlinear_solver_hpp_
#define _nonlinear_solver_hpp_

#include <boost/scoped_ptr.hpp>
#include "gridpack/math/nonlinear_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NonlinearSolver
// -------------------------------------------------------------
/// A solver for a system of nonlinear system of equations in parallel
/**
 * This class is used to solve a system of nonlinear equations in the form

 * \f[
 * \left[ \mathbf{J}\left( \mathbf{x} \right) \right] \Delta \mathbf{x} ~ = ~ 
 *    -\mathbf{F}\left( \mathbf{x} \right) 
 * \f]

 * where \f$\mathbf{J}\left( \mathbf{x} \right)\f$ is the Jacobian
 * matrix, \f$\mathbf{x}\f$ is the solution Vector, and
 * \f$\mathbf{F}\left( \mathbf{x} \right)\f$ is some Vector function
 * of \f$\mathbf{x}\f$. 
 *
 * Users of this class must specify functions or functors that build
 * the \ref JacobianBuilder "Jacobian" Matrix and the \ref
 * FunctionBuilder "right hand side" Vector. Typically, it's best to
 * use functor classes or structs, since extra required information
 * can be available to the Matrix/Vector construction.
 *
 * Implementation ...
 */
template <typename T, typename I = int>
class NonlinearSolverT 
  : public NonlinearSolverInterface<T, I>,
    public parallel::WrappedDistributed,
    public utility::WrappedConfigurable,
    private utility::Uncopyable 
{
public:

  typedef typename NonlinearSolverImplementation<T, I>::VectorType VectorType;
  typedef typename NonlinearSolverImplementation<T, I>::MatrixType MatrixType;
  typedef typename NonlinearSolverImplementation<T, I>::JacobianBuilder JacobianBuilder;
  typedef typename NonlinearSolverImplementation<T, I>::FunctionBuilder FunctionBuilder;

  /// Default constructor.
  /** 
   * @e Collective.
   *
   * A NonlinearSolverT must be constructed simultaneously on all
   * processes involved in @c comm.  
   * 
   * @param comm communicator on which the instance is to exist
   * @param local_size number Jacobian rows and Vector entries to be owned by this process
   * @param form_jacobian function to fill the Jacobian Matrix, \f$\left[ \mathbf{J}\left( \mathbf{x} \right) \right]\f$
   * @param form_function function to fill the RHS function Vector, \f$\mathbf{F}\left( \mathbf{x} \right)\f$
   * 
   * @return new NonlinearSolverT instance
   */
  NonlinearSolverT(const parallel::Communicator& comm,
                  const int& local_size,
                  JacobianBuilder form_jacobian,
                  FunctionBuilder form_function);


  NonlinearSolverT(MatrixType& J,
                   JacobianBuilder form_jacobian,
                   FunctionBuilder form_function);

  /// Destructor
  /**
   * This must be called simultaneously by all processes involved in
   * the \ref parallel::Communicator "communicator" used for \ref
   * NonlinearSolverT() "construction".
   */
  ~NonlinearSolverT(void)
  { }

protected:

  /// A constuctor for children that avoids instantiating an implementation
  NonlinearSolverT(void)
    : NonlinearSolverInterface<T, I>(),
      parallel::WrappedDistributed(),
      utility::WrappedConfigurable(),
      utility::Uncopyable()
  {}


  /// Where things really happen
  /**
   * The Pimpl idiom is used for \ref NonlinearSolverImplementation
   * "implementation", so user code is completely independent of the
   * underlying library. This class simply provides an interface to a
   * specific \ref NonlinearSolverImplementation "implementation".
   * 
   */
  boost::scoped_ptr< NonlinearSolverImplementation<T, I> > p_impl;

  /// Get the solution tolerance (specialized)
  double p_tolerance(void) const
  {
    return p_impl->tolerance();
  }

  /// Set the solver tolerance (specialized)
  void p_tolerance(const double& tol)
  {
    p_impl->tolerance(tol);
  }

  /// Get the maximum iterations (specialized)
  int p_maximumIterations(void) const
  {
    return p_impl->maximumIterations();
  }

  /// Set the maximum solution iterations  (specialized)
  void p_maximumIterations(const int& n) 
  {
    p_impl->maximumIterations(n);
  }

  /// Solve w/ the specified initial estimated, put result in same vector
  void p_solve(VectorType& x)
  {
    p_impl->solve(x);
  }

  /// Set the implementation
  /** 
   * Does what is necessary to set the \ref
   * NonlinearSolverImplementation "implementation".  Subclasses are
   * required to call this at construction.
   * 
   * @param impl specific nonlinear solver implementation to use
   */
  void p_setImpl(NonlinearSolverImplementation<T, I> *impl)
  {
    p_impl.reset(impl);
    p_setDistributed(p_impl.get());
    p_setConfigurable(p_impl.get());
  }

};

typedef NonlinearSolverT<ComplexType> ComplexNonlinearSolver;
typedef NonlinearSolverT<RealType> RealNonlinearSolver;
typedef ComplexNonlinearSolver NonlinearSolver;


} // namespace math
} // namespace gridpack


#endif
