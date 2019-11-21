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
 * @date   2019-11-21 07:41:15 d3g096
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
template <typename T, typename I>
DAESolverT<T, I>::DAESolverT(const parallel::Communicator& comm, 
                             const int local_size,
                             DAESolverT<T, I>::JacobianBuilder& jbuilder,
                             DAESolverT<T, I>::FunctionBuilder& fbuilder,
                             DAESolverT<T, I>::EventManagerPtr eman)
  : parallel::WrappedDistributed(),
    utility::WrappedConfigurable(),
    utility::Uncopyable(),
    p_impl()
{
  p_setImpl(new PETScDAESolverImplementation<T, I>(comm, local_size, jbuilder, fbuilder, eman));
}



template 
DAESolverT<ComplexType>::DAESolverT(const parallel::Communicator& comm, 
                                    const int local_size,
                                    DAESolverT<ComplexType>::JacobianBuilder& jbuilder,
                                    DAESolverT<ComplexType>::FunctionBuilder& fbuilder,
                                    DAESolverT<ComplexType>::EventManagerPtr eman);

template 
DAESolverT<RealType>::DAESolverT(const parallel::Communicator& comm, 
                                 const int local_size,
                                 DAESolverT<RealType>::JacobianBuilder& jbuilder,
                                 DAESolverT<RealType>::FunctionBuilder& fbuilder,
                                 DAESolverT<RealType>::EventManagerPtr eman);

template <typename T, typename I>
DAESolverT<T, I>::DAESolverT(const parallel::Communicator& comm, 
                             const int local_size,
                             DAESolverT<T, I>::JacobianBuilder& jbuilder,
                             DAESolverT<T, I>::FunctionBuilder& fbuilder)
  : parallel::WrappedDistributed(),
    utility::WrappedConfigurable(),
    utility::Uncopyable(),
    p_impl()
{
  p_setImpl(new PETScDAESolverImplementation<T, I>(comm, local_size, jbuilder, fbuilder,
                                                   DAESolverT<T, I>::EventManagerPtr()));
}

template 
DAESolverT<ComplexType>::DAESolverT(const parallel::Communicator& comm, 
                                    const int local_size,
                                    DAESolverT<ComplexType>::JacobianBuilder& jbuilder,
                                    DAESolverT<ComplexType>::FunctionBuilder& fbuilder);

template 
DAESolverT<RealType>::DAESolverT(const parallel::Communicator& comm, 
                                 const int local_size,
                                 DAESolverT<RealType>::JacobianBuilder& jbuilder,
                                 DAESolverT<RealType>::FunctionBuilder& fbuilder);

} // namespace math
} // namespace gridpack
