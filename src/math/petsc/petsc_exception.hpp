/**
 * @file   petsc_exception.hpp
 * @author William A. Perkins
 * @date   2013-05-08 13:44:13 d3g096
 * 
 * @brief 
 * 
 * 
 */

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
