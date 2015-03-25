// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_nonlinear_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2015-03-25 12:15:10 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <boost/format.hpp>

#include "petsc/petsc_nonlinear_solver_implementation.hpp"

namespace gridpack {
namespace math {

PetscErrorCode  
MonitorNorms(SNES snes, PetscInt its, PetscReal fgnorm, void *dummy)
{
  PetscViewer    viewer = PETSC_VIEWER_STDOUT_WORLD;

  Vec dx;
  SNESGetSolutionUpdate(snes, &dx);
  PetscReal dxnorm;
  VecNorm(dx, NORM_2, &dxnorm);
  PetscInt tablevel;
  PetscObjectGetTabLevel(((PetscObject)snes), &tablevel);
  PetscViewerASCIIAddTab(viewer, tablevel);
  PetscViewerASCIIPrintf(viewer,"%3D SNES Function norm %14.12e, Solution residual norm %14.12e\n",
                         its, (double)fgnorm, (double)dxnorm);
  PetscViewerASCIISubtractTab(viewer, tablevel);
  return(0);
}

} // namespace math
} // namespace gridpack
