// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dae_solver.hpp
 * @author William A. Perkins
 * @date   2013-11-13 09:54:32 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _dae_solver_hpp_
#define _dae_solver_hpp_

#include <gridpack/math/dae_solver_implementation.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class DAESolver
// -------------------------------------------------------------
class DAESolver 
  : public parallel::WrappedDistributed,
    public utility::WrappedConfigurable,
    private utility::Uncopyable
{
public:

  /// Default constructor.
  DAESolver(const parallel::Communicator& comm, 
            const int local_size,
            DAEJacobianBuilder& jbuilder,
            DAEFunctionBuilder& fbuilder);

  /// Destructor
  ~DAESolver(void);

  /// Solve the system from @c time to @c end_time
  void solve(const double& time,
             const double& deltat0,
             double& maxtime,
             int& maxsteps,
             Vector& solution)
  {
    p_impl->solve(time, deltat0, maxtime, maxsteps, solution);
  }

protected:
  
  /// Where things really happen
  boost::scoped_ptr<DAESolverImplementation> p_impl;

  /// Set the implementation
  void p_setImpl(DAESolverImplementation *impl);
};



} // namespace math
} // namespace gridpack

#endif
