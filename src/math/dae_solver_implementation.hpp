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
 * @date   2015-05-05 10:08:32 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _dae_solver_implementation_hpp_
#define _dae_solver_implementation_hpp_

#include <boost/shared_ptr.hpp>
#include <gridpack/math/dae_solver_interface.hpp>
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
  : public DAESolverInterface<ComplexType>,
    public parallel::Distributed,
    public utility::Configurable,
    private utility::Uncopyable
{
public:

  typedef typename DAESolverInterface<ComplexType>::VectorType VectorType;
  typedef typename DAESolverInterface<ComplexType>::MatrixType MatrixType;
  typedef typename DAEBuilder<ComplexType>::Jacobian JacobianBuilder;
  typedef typename DAEBuilder<ComplexType>::Function FunctionBuilder;
  

  /// Default constructor.
  DAESolverImplementation(const parallel::Communicator& comm, 
                          const int local_size,
                          JacobianBuilder& jbuilder,
                          FunctionBuilder& fbuilder);

  /// Destructor
  ~DAESolverImplementation(void);

protected:

  /// A matrix to hold the Jacobian
  MatrixType p_J;

  /// A function to build the RHS vector
  FunctionBuilder p_Fbuilder;

  /// A function to build the Jacobian
  JacobianBuilder p_Jbuilder;

  /// Specialized way to configure from property tree
  void p_configure(utility::Configuration::CursorPtr props);

};



} // namespace math
} // namespace gridpack

#endif
