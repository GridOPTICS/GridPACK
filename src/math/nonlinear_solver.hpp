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
 * @date   2013-10-09 12:28:22 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _nonlinear_solver_hpp_
#define _nonlinear_solver_hpp_

#include <boost/scoped_ptr.hpp>
#include "gridpack/math/nonlinear_solver_interface.hpp"

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
class NonlinearSolver 
  : public NonlinearSolverInterface
{
public:

  /// Default constructor.
  /** 
   * @e Collective.
   *
   * A NonlinearSolver must be constructed simultaneously on all
   * processes involved in @c comm.  
   * 
   * @param comm communicator on which the instance is to exist
   * @param local_size number Jacobian rows and Vector entries to be owned by this process
   * @param form_jacobian function to fill the Jacobian Matrix, \f$\left[ \mathbf{J}\left( \mathbf{x} \right) \right]\f$
   * @param form_function function to fill the RHS function Vector, \f$\mathbf{F}\left( \mathbf{x} \right)\f$
   * 
   * @return new NonlinearSolver instance
   */
  NonlinearSolver(const parallel::Communicator& comm,
                  const int& local_size,
                  JacobianBuilder form_jacobian,
                  FunctionBuilder form_function);

  /// Destructor
  /**
   * This must be called simultaneously by all processes involved in
   * the \ref parallel::Communicator "communicator" used for \ref
   * NonlinearSolver() "construction".
   */
  ~NonlinearSolver(void);

};


} // namespace math
} // namespace gridpack


#endif
