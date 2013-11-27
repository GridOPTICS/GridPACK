// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   newton_raphson_solver.cpp
 * @author William A. Perkins
 * @date   2013-11-08 11:51:29 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/assert.hpp>
#include "newton_raphson_solver.hpp"
#include "newton_raphson_solver_implementation.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class NewtonRaphsonSolver
// -------------------------------------------------------------

// -------------------------------------------------------------
// NewtonRaphsonSolver:: constructors / destructor
// -------------------------------------------------------------
NewtonRaphsonSolver::NewtonRaphsonSolver(const parallel::Communicator& comm,
                                         const int& local_size,
                                         JacobianBuilder form_jacobian,
                                         FunctionBuilder form_function)
  : NonlinearSolverInterface()
{
  p_setImpl(new NewtonRaphsonSolverImplementation(comm, local_size,
                                                   form_jacobian,
                                                   form_function));
  p_setDistributed(p_impl.get());
}

NewtonRaphsonSolver::~NewtonRaphsonSolver(void)
{
}

} // namespace math
} // namespace gridpack
