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
 * @date   2013-11-13 09:55:43 d3g096
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

  /// Solve the system
  void solve(const double& time,
             const double& deltat0,
             double& maxtime,
             int& maxsteps,
             Vector& solution)
  {
    this->p_solve(time, deltat0, maxtime, maxsteps, solution);
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

  /// Solve the system (specialized)
  virtual void p_solve(const double& time,
                       const double& deltat0,
                       double& maxtime,
                       int& maxsteps,
                       Vector& solution) = 0;

};



} // namespace math
} // namespace gridpack

#endif
