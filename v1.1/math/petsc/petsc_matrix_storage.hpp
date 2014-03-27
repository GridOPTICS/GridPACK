// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_storage.hpp
 * @author William A. Perkins
 * @date   2013-10-31 10:02:22 d3g096
 * 
 * @brief Helper functions to convert GridPACK/PETSc matrix storage
 * types
 * 
 * 
 */

// -------------------------------------------------------------

#ifndef _petsc_matrix_storage_hpp_
#define _petsc_matrix_storage_hpp_

#include <boost/assert.hpp>
#include <petscmat.h>
#include "matrix.hpp"


namespace gridpack {
namespace math {

// -------------------------------------------------------------
// petscStorageType
// -------------------------------------------------------------
/// 
inline MatType
petscStorageType(const Communicator& comm, const Matrix::StorageType& gtype)
{
  int nproc(comm.size());
  MatType result;
  switch (new_type) {
  case (Matrix::Dense):
    if (nproc > 1) {
      new_mat_type = MATMPIDENSE;
    } else {
      new_mat_type = MATSEQDENSE;
    } 
    break;
  case (Matrix::Sparse):
    if (nproc > 1) {
      new_mat_type = MATMPIAIJ;
    } else {
      new_mat_type = MATSEQAIJ;
    } 
    break;
   default:
     BOOST_ASSERT(false); 
  }
}


} // namespace math
} // namespace gridpack
#endif
