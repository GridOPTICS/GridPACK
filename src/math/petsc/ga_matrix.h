/* -------------------------------------------------------------
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 * ------------------------------------------------------------- */
/**
 * @file   ga_matrix.h
 * @author William A. Perkins
 * @date   2015-04-29 14:28:44 d3g096
 * 
 * @brief  
 * 
 * 
 */

#ifndef _ga_matrix_h_
#define _ga_matrix_h_

#include "petscmat.h" 


/** 
 * Create a dense PETSc matrix that uses Global Arrays for data
 * storage.
 * 
 * @param comm 
 * @param m local number of rows (or PETSC_DECIDE)
 * @param n local number of columns (or PETSC_DECIDE)
 * @param M total number of rows (or PETSC_DETERMINE)
 * @param N total number of columns (or PETSC_DETERMINE)
 * @param A new PETSc matrix instance 
 * 
 * @return 
 */
PETSC_EXTERN
PetscErrorCode
MatCreateDenseGA(MPI_Comm comm, 
                 PetscInt m,PetscInt n,PetscInt M,PetscInt N, 
                 Mat *A);

#endif
