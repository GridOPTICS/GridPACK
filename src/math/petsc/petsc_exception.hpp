// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_exception.hpp
 * @author William A. Perkins
 * @date   2013-10-09 13:23:43 d3g096
 * 
 * @brief 
 * 
 * 
 */
// -------------------------------------------------------------

#include <string>
#include <petscsys.h>
#include <petscsys.hh>
#include "gridpack/utilities/exception.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScException
// -------------------------------------------------------------
class PETScException : public gridpack::Exception {
public:

  /// Construct with a PETSc error code.
  explicit PETScException(const PetscErrorCode ierr);

  /// Construct from a PETSc exception
  explicit PETScException(const PETSc::Exception& e);

  /// Construct from a PETSc exception and error code
  PETScException(const PetscErrorCode& ierr, const PETSc::Exception& e);

  /// Copy constructor
  PETScException(const PETScException& old);

  /// Destructor
  ~PETScException(void) throw();

protected:
  
  /// PETSc error code causing the problem
  PetscErrorCode petsc_err_;

  
};





} // namespace math
} // namespace gridpack
