// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   dae_solver.cpp
 * @author William A. Perkins
 * @date   2015-05-05 10:13:34 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "dae_solver.hpp"
#include "petsc/petsc_dae_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class DAESolver
// -------------------------------------------------------------

// -------------------------------------------------------------
// DAESolver:: constructors / destructor
// -------------------------------------------------------------
DAESolver::DAESolver(const parallel::Communicator& comm, 
                     const int local_size,
                     DAESolverImplementation::JacobianBuilder& jbuilder,
                     DAESolverImplementation::FunctionBuilder& fbuilder)
  : parallel::WrappedDistributed(),
    utility::WrappedConfigurable(),
    utility::Uncopyable(),
    p_impl()
{
  p_setImpl(new PETScDAESolverImplementation(comm, local_size,
                                             jbuilder, fbuilder));
}

} // namespace math
} // namespace gridpack
