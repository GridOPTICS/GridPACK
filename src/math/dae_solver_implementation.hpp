// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dae_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2013-11-14 11:39:40 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _dae_solver_implementation_hpp_
#define _dae_solver_implementation_hpp_

#include <boost/shared_ptr.hpp>
#include <gridpack/math/dae_solver_functions.hpp>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/configuration/configurable.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class DAESolverImplementation
// -------------------------------------------------------------
class DAESolverImplementation 
  : public parallel::Distributed,
    public utility::Configurable,
    private utility::Uncopyable
{
public:

  /// Default constructor.
  DAESolverImplementation(const parallel::Communicator& comm, 
                          const int local_size,
                          DAEJacobianBuilder& jbuilder,
                          DAEFunctionBuilder& fbuilder);

  /// Destructor
  ~DAESolverImplementation(void);

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
                  Vector& x0)
  {
    this->p_initialize(t0, deltat0, x0);
  }

  /// Solve the system
  void solve(double& maxtime, int& maxsteps)
  {
    this->p_solve(maxtime, maxsteps);
  }

protected:

  /// A matrix to hold the Jacobian
  Matrix p_J;

  /// A function to build the RHS vector
  DAEFunctionBuilder p_Fbuilder;

  /// A function to build the Jacobian
  DAEJacobianBuilder p_Jbuilder;

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props);

  /// Initialize the system (specialized)
  virtual void p_initialize(const double& t0,
                            const double& deltat0,
                            Vector& x0) = 0;
                       

  /// Solve the system (specialized)
  virtual void p_solve(double& maxtime, int& maxsteps) = 0;

};



} // namespace math
} // namespace gridpack

#endif
