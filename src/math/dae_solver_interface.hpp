// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   dae_solver_interface.hpp
 * @author William A. Perkins
 * @date   2015-05-05 11:44:17 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _dae_solver_interface_hpp_
#define _dae_solver_interface_hpp_

#include <gridpack/math/vector.hpp>
#include <gridpack/math/matrix.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class DAESolverInterface
// -------------------------------------------------------------
template <typename T, typename I = int>
class DAESolverInterface 
{
public:

  typedef VectorT<T, I> VectorType;
  typedef MatrixT<T, I> MatrixType;

  /// Default constructor.
  DAESolverInterface(void)
  {
  }

  /// Destructor
  virtual ~DAESolverInterface(void)
  {
  }

  /// Initialize the solver
  /** 
   * 
   * 
   * @param t0 start time, corresponding to @c x0
   * @param deltat0 initial time step
   * @param x0 initial solution corresponding to @c t0
   */
  void initialize(const double& t0,
                  const double& deltat0,
                  VectorType& x0)
  {
    this->p_initialize(t0, deltat0, x0);
  }

  /// Solve the system to when @c end_time or @c maxsteps is exceeded
  /** 
   * This solves the system from the initial time specified by
   * ::initialize() or continues the solution from the end time of the
   * previous call to ::solve().  After the solution completes, the
   * solution is placed back in the Vector passed to ::initialize()
   * 
   * @param maxtime time at which solution is to end (in); actual time
   * at which solution ended (out)
   *
   * @param maxsteps maximum number of time steps to perform (in);
   * actual number of time steps performed
   *
   * @param solution system solution at @c maxtime (out)
   */
  void solve(double& maxtime, int& maxsteps)
  {
    this->p_solve(maxtime, maxsteps);
  }

protected:

  /// Initialize the system (specialized)
  virtual void p_initialize(const double& t0,
                            const double& deltat0,
                            VectorType& x0) = 0;
                       

  /// Solve the system (specialized)
  virtual void p_solve(double& maxtime, int& maxsteps) = 0;
};



} // namespace math
} // namespace gridpack

#endif
