// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   nonlinear_solver_c.cpp
 * @author William A. Perkins
 * @date   2015-05-06 09:41:07 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created September  2, 2014 by William A. Perkins
// -------------------------------------------------------------

#include <boost/utility.hpp>
#include <boost/scoped_ptr.hpp>
#include "configuration/cursor_wrapper.hpp"
#include <gridpack/math/nonlinear_solver.hpp>
#include <gridpack/math/newton_raphson_solver.hpp>

// callable functions in the Fortran module

/// Call the Fortran thing to build the Jacobian
extern "C" void
builder_build_jacobian(void *bptr, 
                       gridpack::math::Matrix *J,
                       const gridpack::math::Vector *x);

/// Call the Fortran thing to build the RHS
extern "C" void
builder_build_function(void *bptr, 
                       gridpack::math::Vector *F,
                       const gridpack::math::Vector *x);


// -------------------------------------------------------------
//  class FortranBuilder
// -------------------------------------------------------------
class FortranBuilder 
  : private boost::noncopyable {
public:

  /// Default constructor.
  FortranBuilder(void *fbuilder)
    : p_fortranBuilder(fbuilder)
  {}

  /// Destructor
  ~FortranBuilder(void)
  {}

  /// Operator to build Jacobian
  void operator()(const gridpack::math::Vector& x, gridpack::math::Matrix& J)
  {
    builder_build_jacobian(p_fortranBuilder, &J, &x);
    x.print();
    J.print();
  }

  /// Operator to build RHS
  void operator()(const gridpack::math::Vector& x, gridpack::math::Vector& F) 
  {
    builder_build_function(p_fortranBuilder, &F, &x);
    F.print();
  }

protected:

  /// A pointer to the Fortran thing used to build the Jacobian and RHS
  void *p_fortranBuilder;

};


// -------------------------------------------------------------
// struct NonlinearSolverWrapper
// -------------------------------------------------------------
struct NonlinearSolverWrapper {
  boost::scoped_ptr< gridpack::math::NonlinearSolver > impl;
  boost::scoped_ptr<FortranBuilder> bldr;
};




extern "C" NonlinearSolverWrapper *
nonlinear_solver_create(gridpack::parallel::Communicator *comm,
                        CursorWrapper *conf,
                        int local_size,
                        void *builder)
{
  NonlinearSolverWrapper *result(new NonlinearSolverWrapper);
  result->bldr.reset(new FortranBuilder(builder));
  gridpack::math::NonlinearSolver::JacobianBuilder j = boost::ref(*(result->bldr));
  gridpack::math::NonlinearSolver::FunctionBuilder f = boost::ref(*(result->bldr));
  result->impl.reset(new gridpack::math::NonlinearSolver(*comm, local_size, j, f));
  gridpack::utility::Configuration::CursorPtr cursor;
  if (conf != NULL) {
    cursor = conf->cursor;
  } else {
    cursor = gridpack::utility::Configuration::configuration()->getCursor("");
  }
  result->impl->configure(cursor);
  return result;
}

extern "C" NonlinearSolverWrapper *
newton_solver_create(gridpack::parallel::Communicator *comm,
                     CursorWrapper *conf,
                     int local_size,
                     void *builder)
{
  NonlinearSolverWrapper *result(new NonlinearSolverWrapper);
  result->bldr.reset(new FortranBuilder(builder));
  gridpack::math::NewtonRaphsonSolver::JacobianBuilder j = boost::ref(*(result->bldr));
  gridpack::math::NewtonRaphsonSolver::FunctionBuilder f = boost::ref(*(result->bldr));
  result->impl.reset(new gridpack::math::NewtonRaphsonSolver(*comm, local_size, j, f));
  gridpack::utility::Configuration::CursorPtr cursor;
  if (conf != NULL) {
    cursor = conf->cursor;
  } else {
    cursor = gridpack::utility::Configuration::configuration()->getCursor("");
  }
  result->impl->configure(cursor);
  return result;
}
                          
extern "C" void
nonlinear_solver_destroy(NonlinearSolverWrapper **ptr)
{
  delete *ptr;
  *ptr = NULL;
}

extern "C" void
nonlinear_solver_solve(NonlinearSolverWrapper *ptr,
                       gridpack::math::Vector *x)
{
  ptr->impl->solve(*x);
}
