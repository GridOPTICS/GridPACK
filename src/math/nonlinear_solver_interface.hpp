// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver_interface.hpp
 * @author William A. Perkins
 * @date   2013-11-07 12:33:31 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _nonlinear_solver_interface_hpp_
#define _nonlinear_solver_interface_hpp_

#include <boost/scoped_ptr.hpp>
#include <gridpack/math/nonlinear_solver_implementation.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NonlinearSolverInterface
// -------------------------------------------------------------
/// Interface to a solve system of nonlinear equations in parallel
/**
 * This serves as a base for classes that solve a system of nonlinear
 * equations. While not strictly abstract, it has no function if
 * instantiated on its own.
 * 
 * It encapuslates the nonlinear system solver of some
 * underlying implementation. The Pimpl idiom is used for \ref
 * NonlinearSolverImplementation "implementation", so user code is
 * completely independent of the underlying library. This class simply
 * provides an interface to a specific \ref
 * NonlinearSolverImplementation "implementation".  Subclasses are
 * required to call p_set_impl() at construction to set the \ref
 * NonlinearSolverImplementation "implementation".
 */
class NonlinearSolverInterface 
  : public parallel::WrappedDistributed,
    public utility::WrappedConfigurable,
    private utility::Uncopyable 
{
public:

  /// Default constructor.
  NonlinearSolverInterface();

  /// Destructor
  ~NonlinearSolverInterface(void);

  /// Solve w/ the specified initial estimated, put result in same vector
  /** 
   * This solves the system of nonlinear equations using the contents
   * of @c x as an initial solution estimate.  The final result is
   * placed back in @c x upon completion.
   * 
   * @param x solution Vector
   */
  void solve(Vector& x)
  {
    p_impl->solve(x);
  }

protected:

  /// Where things really happen
  /**
   * The Pimpl idiom is used for \ref NonlinearSolverImplementation
   * "implementation", so user code is completely independent of the
   * underlying library. This class simply provides an interface to a
   * specific \ref NonlinearSolverImplementation "implementation".
   * 
   */
  boost::scoped_ptr<NonlinearSolverImplementation> p_impl;
  
  /// Set the implementation
  /** 
   * Does what is necessary to set the \ref
   * NonlinearSolverImplementation "implementation".  Subclasses are
   * required to call this at construction.
   * 
   * @param impl specific nonlinear solver implementation to use
   */
  void p_set_impl(NonlinearSolverImplementation *impl);
};


} // namespace math
} // namespace gridpack


#endif
