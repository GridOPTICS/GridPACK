// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   linear_solverc.cpp
 * @author William A. Perkins
 * @date   2014-05-29 14:30:03 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// Created May 29, 2014 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#include <gridpack/math/linear_solver.hpp>


extern "C" void
linear_solver_initialize(gridpack::math::LinearSolver **s,
                         gridpack::math::Matrix *A)
{
  *s = new gridpack::math::LinearSolver(*A);
  gridpack::utility::Configuration::CursorPtr 
    empty(new gridpack::utility::Configuration);
  (*s)->configure(empty);
}

extern "C" void
linear_solver_finalize(gridpack::math::LinearSolver **s)
{
  delete *s;
  *s = NULL;
}

extern "C" void
linear_solver_solve(gridpack::math::LinearSolver *s,
                    gridpack::math::Vector *b,
                    gridpack::math::Vector *x)
{
  s->solve(*b, *x);
}


extern "C" void
linear_solver_set_tolerance(gridpack::math::LinearSolver *s,
                            double t)
{
  s->tolerance(t);
}

extern "C" void
linear_solver_set_iterations(gridpack::math::LinearSolver *s,
                             int n)
{
  s->maximumIterations(n);
}
