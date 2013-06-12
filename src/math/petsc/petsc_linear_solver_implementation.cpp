/**
 * @file   petsc_linear_solver_implementation.cpp
 * @author William A. Perkins
 * @date   2013-06-11 14:08:23 d3g096
 * 
 * @brief  
 * 
 * 
 */

#include "petsc_linear_solver_implementation.hpp"

namespace gridpack {
namespace math {


// -------------------------------------------------------------
//  class PETScLinearSolverImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// PETScLinearSolverImplementation:: constructors / destructor
// -------------------------------------------------------------
PETScLinearSolverImplementation::PETScLinearSolverImplementation(const Matrix& A)
  : 
{
  
}

PETScLinearSolverImplementation::~PETScLinearSolverImplementation(void)
{
  PetscErrorCode ierr;
  try  {
    PetscBool ok;
    ierr = PetscInitialized(&ok);
    if (ok) {
      ierr = KSPDestroy(&p_KSP);
    }
  } catch (...) {
    // just eat it
  }
}

// -------------------------------------------------------------
// PETScLinearSolverImplementation::p_accept
// -------------------------------------------------------------
void 
PETScLinearSolverImplementation::p_accept(ImplementationVisitor& visitor)
{
  visitor.visit(*this);
}

void 
PETScLinearSolverImplementation::p_accept(ConstImplementationVisitor& visitor) const
{
  visitor.visit(*this);
}




} // namespace math
} // namespace gridpack
