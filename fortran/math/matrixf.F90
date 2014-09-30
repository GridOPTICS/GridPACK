! ----------------------------------------------------------------
! file: matrixf.f90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created May 27, 2014 by William A. Perkins
! Last Change: 2014-08-27 12:19:51 d3g096
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! MODULE gridpack_matrix
! ----------------------------------------------------------------
MODULE gridpack_matrix

  USE iso_c_binding
  USE gridpack_communicator
  USE gridpack_vector

  IMPLICIT NONE
  
  PRIVATE

  TYPE, PUBLIC :: matrix
     TYPE (c_ptr) :: mat
   CONTAINS
     PROCEDURE :: initialize
     PROCEDURE :: finalize      ! should be FINAL
     PROCEDURE :: rows
     PROCEDURE :: local_rows
     PROCEDURE :: local_row_range
     PROCEDURE :: cols
     PROCEDURE :: set_element
     PROCEDURE :: add_element
     PROCEDURE :: get_element
     PROCEDURE :: identity
     PROCEDURE :: real
     PROCEDURE :: imaginary
     PROCEDURE :: conjugate
     PROCEDURE :: ready
     PROCEDURE :: clone
     PROCEDURE :: print
     PROCEDURE :: scale
     PROCEDURE :: add
     PROCEDURE :: zero
  END type matrix

  INTERFACE
     SUBROUTINE matrix_initialize(mat, comm, local_rows, local_cols, nz_per_row) BIND(c)
       USE iso_c_binding, ONLY: c_int, c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(INOUT) :: mat
       TYPE (c_ptr), VALUE, INTENT(IN) :: comm
       INTEGER(c_int), VALUE, INTENT(IN) :: local_rows, local_cols, nz_per_row
     END SUBROUTINE matrix_initialize

     SUBROUTINE matrix_finalize(mat) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(INOUT) :: mat
     END SUBROUTINE matrix_finalize

     INTEGER(c_int) FUNCTION matrix_rows(mat) BIND(c)
       USE iso_c_binding, ONLY: c_int, c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
     END FUNCTION matrix_rows
     
     INTEGER(c_int) FUNCTION matrix_local_rows(mat) BIND(c)
       USE iso_c_binding, ONLY: c_int, c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
     END FUNCTION matrix_local_rows

     SUBROUTINE matrix_local_row_range(mat, lo, hi) BIND(c)
       USE iso_c_binding, ONLY: c_int, c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
       INTEGER (c_int), INTENT(OUT) :: lo, hi
     END SUBROUTINE matrix_local_row_range
     
     INTEGER(c_int) FUNCTION matrix_cols(mat) BIND(c)
       USE iso_c_binding, ONLY: c_int, c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
     END FUNCTION matrix_cols

     SUBROUTINE matrix_set_element(mat, i, j, x) BIND(c)
       USE iso_c_binding, ONLY: c_int, c_ptr, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
       INTEGER (c_int), VALUE, INTENT(IN) :: i, j
       COMPLEX(c_double_complex), INTENT(IN) :: x
     END SUBROUTINE matrix_set_element

     SUBROUTINE matrix_add_element(mat, i, j, x) BIND(c)
       USE iso_c_binding, ONLY: c_int, c_ptr, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
       INTEGER (c_int), VALUE, INTENT(IN) :: i, j
       COMPLEX(c_double_complex), INTENT(IN) :: x
     END SUBROUTINE matrix_add_element
     
     SUBROUTINE matrix_get_element(mat, i, j, x) BIND(c)
       USE iso_c_binding, ONLY: c_int, c_ptr, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
       INTEGER (c_int), VALUE, INTENT(IN) :: i, j
       COMPLEX(c_double_complex), INTENT(OUT) :: x
     END SUBROUTINE matrix_get_element

     SUBROUTINE matrix_identity(mat) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
     END SUBROUTINE matrix_identity

     SUBROUTINE matrix_real(mat) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
     END SUBROUTINE matrix_real

     SUBROUTINE matrix_imaginary(mat) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
     END SUBROUTINE matrix_imaginary

     SUBROUTINE matrix_conjugate(mat) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
     END SUBROUTINE matrix_conjugate
     
     SUBROUTINE matrix_ready(mat) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
     END SUBROUTINE matrix_ready

     TYPE (c_ptr) FUNCTION matrix_clone(mat) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
     END FUNCTION matrix_clone

     SUBROUTINE matrix_print(mat) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
     END SUBROUTINE matrix_print

     SUBROUTINE matrix_scale(mat, x) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
       COMPLEX(c_double_complex), INTENT(IN) :: x
     END SUBROUTINE matrix_scale

     SUBROUTINE matrix_multiply_diagonal(mat, x) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
       COMPLEX(c_double_complex), INTENT(IN) :: x
     END SUBROUTINE matrix_multiply_diagonal

     SUBROUTINE matrix_add(mat, a) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat, a
     END SUBROUTINE matrix_add

     SUBROUTINE matrix_zero(mat) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: mat
     END SUBROUTINE matrix_zero

  END INTERFACE

  TYPE, PUBLIC :: matrix_wrap
     CLASS(matrix), POINTER :: matrix
  END type matrix_wrap

CONTAINS

  ! ----------------------------------------------------------------
  ! SUBROUTINE initialize
  ! ----------------------------------------------------------------
  SUBROUTINE initialize(this, comm, local_rows, local_cols, nz_per_row)
    IMPLICIT NONE
    CLASS (matrix), INTENT(INOUT) :: this
    CLASS (communicator), INTENT(IN) :: comm
    INTEGER, INTENT(IN) :: local_rows, local_cols, nz_per_row
    INTEGER (c_int) :: crows, ccols, cnz
    crows = local_rows
    ccols = local_cols
    cnz = nz_per_row
    CALL matrix_initialize(this%mat, comm%comm, crows, ccols, cnz)
  END SUBROUTINE initialize

  ! ----------------------------------------------------------------
  ! SUBROUTINE finalize
  ! ----------------------------------------------------------------
  SUBROUTINE finalize(this)
    IMPLICIT NONE
    CLASS (matrix), INTENT(INOUT) :: this
    CALL matrix_finalize(this%mat)
  END SUBROUTINE finalize



  ! ----------------------------------------------------------------
  ! INTEGER FUNCTION rows
  ! ----------------------------------------------------------------
  INTEGER FUNCTION rows(this)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    rows = matrix_rows(this%mat)
  END FUNCTION rows

  ! ----------------------------------------------------------------
  ! INTEGER FUNCTION local_rows
  ! ----------------------------------------------------------------
  INTEGER FUNCTION local_rows(this)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    local_rows = matrix_local_rows(this%mat)
  END FUNCTION local_rows

  ! ----------------------------------------------------------------
  ! SUBROUTINE local_row_range
  ! ----------------------------------------------------------------
  SUBROUTINE local_row_range(this, lo, hi)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    INTEGER, INTENT(OUT) :: lo, hi
    INTEGER (c_int) :: clo, chi
    CALL matrix_local_row_range(this%mat, clo, chi)
    lo = clo
    hi = chi
  END SUBROUTINE local_row_range


  ! ----------------------------------------------------------------
  ! INTEGER FUNCTION cols
  ! ----------------------------------------------------------------
  INTEGER FUNCTION cols(this)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    cols = matrix_cols(this%mat)
  END FUNCTION cols

  ! ----------------------------------------------------------------
  ! SUBROUTINE set_element
  ! ----------------------------------------------------------------
  SUBROUTINE set_element(this, i, j, x)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: i, j
    COMPLEX(c_double_complex), INTENT(IN) :: x
    INTEGER(c_int) :: ci, cj
    ci = i
    cj = j
    CALL matrix_set_element(this%mat, ci, cj, x)
  END SUBROUTINE set_element

  ! ----------------------------------------------------------------
  ! SUBROUTINE add_element
  ! ----------------------------------------------------------------
  SUBROUTINE add_element(this, i, j, x)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: i, j
    COMPLEX(c_double_complex), INTENT(IN) :: x
    INTEGER(c_int) :: ci, cj
    ci = i
    cj = j
    CALL matrix_add_element(this%mat, ci, cj, x)
  END SUBROUTINE add_element

  ! ----------------------------------------------------------------
  ! SUBROUTINE get_element
  ! ----------------------------------------------------------------
  SUBROUTINE get_element(this, i, j, x)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: i, j
    COMPLEX(c_double_complex), INTENT(OUT) :: x
    INTEGER(c_int) :: ci, cj
    ci = i
    cj = j
    CALL matrix_get_element(this%mat, ci, cj, x)
  END SUBROUTINE get_element

  ! ----------------------------------------------------------------
  ! SUBROUTINE identity
  ! ----------------------------------------------------------------
  SUBROUTINE identity(this)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    CALL matrix_identity(this%mat)
  END SUBROUTINE identity

  ! ----------------------------------------------------------------
  ! SUBROUTINE real
  ! ----------------------------------------------------------------
  SUBROUTINE real(this)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    CALL matrix_real(this%mat)
  END SUBROUTINE real

  ! ----------------------------------------------------------------
  ! SUBROUTINE imaginary
  ! ----------------------------------------------------------------
  SUBROUTINE imaginary(this)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    CALL matrix_imaginary(this%mat)
  END SUBROUTINE imaginary

  ! ----------------------------------------------------------------
  ! SUBROUTINE conjugate
  ! ----------------------------------------------------------------
  SUBROUTINE conjugate(this)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    CALL matrix_conjugate(this%mat)
  END SUBROUTINE conjugate

  ! ----------------------------------------------------------------
  ! SUBROUTINE ready
  ! ----------------------------------------------------------------
  SUBROUTINE ready(this)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    CALL matrix_ready(this%mat)
  END SUBROUTINE ready

  ! ----------------------------------------------------------------
  ! TYPE (C_PTR) FUNCTION clone
  ! ----------------------------------------------------------------
  TYPE (matrix) FUNCTION clone(this)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    clone%mat = matrix_clone(this%mat)
  END FUNCTION clone

  ! ----------------------------------------------------------------
  ! SUBROUTINE print
  ! ----------------------------------------------------------------
  SUBROUTINE print(this)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    CALL matrix_print(this%mat)
  END SUBROUTINE print

  ! ----------------------------------------------------------------
  ! SUBROUTINE scale
  ! ----------------------------------------------------------------
  SUBROUTINE scale(this, x)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    COMPLEX (c_double_complex), INTENT(IN) :: x
    CALL matrix_scale(this%mat, x)
  END SUBROUTINE scale


  ! ----------------------------------------------------------------
  ! SUBROUTINE add
  ! ----------------------------------------------------------------
  SUBROUTINE add(this, a)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this, a
    CALL matrix_add(this%mat, a%mat)
  END SUBROUTINE add


  ! ----------------------------------------------------------------
  ! SUBROUTINE zero
  ! ----------------------------------------------------------------
  SUBROUTINE zero(this)
    IMPLICIT NONE
    CLASS (matrix), INTENT(IN) :: this
    CALL matrix_zero(this%mat)
  END SUBROUTINE zero


  

END MODULE gridpack_matrix
