! ----------------------------------------------------------------
! file: vectorf.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created May 15, 2014 by William A. Perkins
! Last Change: 2014-10-01 07:31:27 d3g096
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! MODULE gridpack_vector
! ----------------------------------------------------------------
MODULE gridpack_vector

  USE iso_c_binding
  USE gridpack_parallel

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: vector
     TYPE (c_ptr) :: vec
   CONTAINS
     PROCEDURE :: initialize 
     PROCEDURE :: finalize ! should be FINAL
     PROCEDURE :: size
     PROCEDURE :: local_size
     PROCEDURE :: local_index_range
     PROCEDURE :: set_element
     PROCEDURE :: add_element
     PROCEDURE :: get_element
     PROCEDURE :: get_all_elements
     PROCEDURE :: zero
     PROCEDURE :: fill 
     PROCEDURE :: norm1
     PROCEDURE :: norm2
     PROCEDURE :: norm_infinity
     PROCEDURE :: ready
     PROCEDURE :: clone
     PROCEDURE :: print
  END type vector

  INTERFACE
     SUBROUTINE vector_initialize(vec, comm, local_size) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_int
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(INOUT) :: vec
       TYPE (c_ptr), VALUE, INTENT(IN) :: comm
       INTEGER (c_int), VALUE, INTENT(IN) :: local_size
     END SUBROUTINE vector_initialize

     SUBROUTINE vector_finalize(vec) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(INOUT) :: vec
     END SUBROUTINE vector_finalize

     INTEGER(c_int) FUNCTION vector_size(vec) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_int
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
     END FUNCTION vector_size

     INTEGER(c_int) FUNCTION vector_local_size(vec) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_int
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
     END FUNCTION vector_local_size

     SUBROUTINE vector_local_index_range(vec, lo, hi) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_int
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
       INTEGER (c_int), INTENT(OUT) :: lo, hi
     END SUBROUTINE vector_local_index_range

     SUBROUTINE vector_set_element(vec, i, x) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_int, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
       INTEGER (c_int), VALUE, INTENT(IN) :: i
       COMPLEX(c_double_complex), INTENT(IN) :: x
     END SUBROUTINE vector_set_element

     SUBROUTINE vector_add_element(vec, i, x) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_int, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
       INTEGER (c_int), VALUE, INTENT(IN) :: i
       COMPLEX(c_double_complex), INTENT(IN) :: x
     END SUBROUTINE vector_add_element

     SUBROUTINE vector_get_element(vec, i, x) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_int, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
       INTEGER (c_int), VALUE, INTENT(IN) :: i
       COMPLEX(c_double_complex), INTENT(OUT) :: x
     END SUBROUTINE vector_get_element

     SUBROUTINE vector_get_all_elements(vec, x) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
       TYPE (c_ptr), VALUE, INTENT(IN) :: x
     END SUBROUTINE vector_get_all_elements

     SUBROUTINE vector_zero(vec) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
     END SUBROUTINE vector_zero

     SUBROUTINE vector_fill(vec, x) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
       COMPLEX(c_double_complex), VALUE, INTENT(IN) :: x
     END SUBROUTINE vector_fill

     COMPLEX(c_double_complex) FUNCTION vector_norm1(vec) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
     END FUNCTION vector_norm1

     COMPLEX(c_double_complex) FUNCTION vector_norm2(vec) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
     END FUNCTION vector_norm2

     COMPLEX(c_double_complex) FUNCTION vector_norm_infinity(vec) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_double_complex
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
     END FUNCTION vector_norm_infinity

     SUBROUTINE vector_ready(vec) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
     END SUBROUTINE vector_ready

     SUBROUTINE vector_print(vec) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
     END SUBROUTINE vector_print

     TYPE (c_ptr) FUNCTION vector_clone(vec) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: vec
     END FUNCTION vector_clone

  END INTERFACE

  TYPE, PUBLIC :: vector_wrap
     CLASS(vector), POINTER :: vector
  END type vector_wrap

CONTAINS

  ! ----------------------------------------------------------------
  ! SUBROUTINE initialize
  ! ----------------------------------------------------------------
  SUBROUTINE initialize(this, comm, local_size)
    IMPLICIT NONE
    CLASS (vector), INTENT(INOUT) :: this
    CLASS (communicator), INTENT(IN) :: comm
    INTEGER, INTENT(IN) :: local_size
    INTEGER(c_int) :: lsize
    lsize = local_size
    CALL vector_initialize(this%vec, comm%comm, lsize)
  END SUBROUTINE initialize

  ! ----------------------------------------------------------------
  ! SUBROUTINE finalize
  ! ----------------------------------------------------------------
  SUBROUTINE finalize(this)
    IMPLICIT NONE
    CLASS (vector), INTENT(INOUT) :: this
    CALL vector_finalize(this%vec)
  END SUBROUTINE finalize

  ! ----------------------------------------------------------------
  ! INTEGER FUNCTION size
  ! ----------------------------------------------------------------
  INTEGER FUNCTION size(this)
    IMPLICIT NONE
    CLASS (vector), INTENT(IN) :: this
    size = vector_size(this%vec)
  END FUNCTION size

  ! ----------------------------------------------------------------
  ! INTEGER FUNCTION local_size
  ! ----------------------------------------------------------------
  INTEGER FUNCTION local_size(this)
    IMPLICIT NONE
    CLASS (vector), INTENT(IN) :: this
    local_size = vector_local_size(this%vec)
  END FUNCTION local_size

  ! ----------------------------------------------------------------
  ! SUBROUTINE local_index_range
  ! ----------------------------------------------------------------
  SUBROUTINE local_index_range(this, lo, hi)
    IMPLICIT NONE
    CLASS (vector), INTENT(IN) :: this
    INTEGER, INTENT(OUT) :: lo, hi
    INTEGER(c_int) :: clo, chi

    clo = 0
    chi = 0
    CALL vector_local_index_range(this%vec, clo, chi)
    lo = clo
    hi = chi
    ! WRITE(*,*) clo, chi, lo, hi

  END SUBROUTINE local_index_range

  ! ----------------------------------------------------------------
  ! SUBROUTINE set_element
  ! ----------------------------------------------------------------
  SUBROUTINE set_element(this, i, x)
    IMPLICIT NONE
    CLASS (vector), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: i
    COMPLEX(c_double_complex), INTENT(IN) :: x
    INTEGER(c_int) ci
    ci = i
    CALL vector_set_element(this%vec, ci, x)
  END SUBROUTINE set_element

  ! ----------------------------------------------------------------
  ! SUBROUTINE add_element
  ! ----------------------------------------------------------------
  SUBROUTINE add_element(this, i, x)
    IMPLICIT NONE
    CLASS (vector), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: i
    COMPLEX(c_double_complex), INTENT(IN) :: x
    INTEGER(c_int) ci
    ci = i
    CALL vector_add_element(this%vec, ci, x)
  END SUBROUTINE add_element

  ! ----------------------------------------------------------------
  ! SUBROUTINE get_element
  ! ----------------------------------------------------------------
  SUBROUTINE get_element(this, i, x)
    IMPLICIT NONE
    CLASS (vector), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: i
    COMPLEX(c_double_complex), INTENT(OUT) :: x
    INTEGER(c_int) ci
    ci = i
    CALL vector_get_element(this%vec, ci, x)
  END SUBROUTINE get_element

  ! ----------------------------------------------------------------
  ! SUBROUTINE get_all_elements
  ! ----------------------------------------------------------------
  SUBROUTINE get_all_elements(this, xout)
    IMPLICIT NONE
    CLASS (vector), INTENT(IN) :: this
    COMPLEX(c_double_complex), TARGET, INTENT(OUT) :: xout(*)
    TYPE (c_ptr) :: xptr
    xptr = c_loc(xout)
    CALL vector_get_all_elements(this%vec, xptr)
  END SUBROUTINE get_all_elements


  ! ----------------------------------------------------------------
  ! SUBROUTINE zero
  ! ----------------------------------------------------------------
  SUBROUTINE zero(this)
    IMPLICIT NONE
    CLASS (vector), INTENT(IN) :: this
    CALL vector_zero(this%vec)
  END SUBROUTINE zero

  ! ----------------------------------------------------------------
  ! SUBROUTINE fill
  ! ----------------------------------------------------------------
  SUBROUTINE fill(this, x)
    IMPLICIT NONE
    CLASS (vector), INTENT(IN) :: this
    COMPLEX(c_double_complex), INTENT(IN) :: x
    CALL vector_fill(this%vec, x)
  END SUBROUTINE fill

  ! ----------------------------------------------------------------
  ! COMPLEX(c_double_complex) FUNCTION norm1
  ! ----------------------------------------------------------------
  COMPLEX(c_double_complex) FUNCTION norm1(this)
    IMPLICIT NONE
    CLASS (vector), INTENT(IN) :: this
    norm1 = vector_norm1(this%vec)
  END FUNCTION norm1

  ! ----------------------------------------------------------------
  ! COMPLEX(c_double_complex) FUNCTION norm1
  ! ----------------------------------------------------------------
  COMPLEX(c_double_complex) FUNCTION norm2(this)
    IMPLICIT NONE
    CLASS (vector), INTENT(IN) :: this
    norm2 = vector_norm2(this%vec)
  END FUNCTION norm2

  ! ----------------------------------------------------------------
  ! COMPLEX(c_double_complex) FUNCTION norm_infinity
  ! ----------------------------------------------------------------
  COMPLEX(c_double_complex) FUNCTION norm_infinity(this)
    IMPLICIT NONE
    CLASS (vector), INTENT(IN) :: this
    norm_infinity = vector_norm_infinity(this%vec)
  END FUNCTION norm_infinity


  ! ----------------------------------------------------------------
  ! SUBROUTINE ready
  ! ----------------------------------------------------------------
  SUBROUTINE ready(this)
    IMPLICIT NONE
    CLASS (vector), INTENT(INOUT) :: this
    CALL vector_ready(this%vec)
  END SUBROUTINE ready

  ! ----------------------------------------------------------------
  ! TYPE (vector) FUNCTION clone
  ! ----------------------------------------------------------------
  TYPE (vector) FUNCTION clone(this)
    IMPLICIT NONE
    CLASS (vector), INTENT(INOUT) :: this
    clone%vec = vector_clone(this%vec)
  END FUNCTION clone

  ! ----------------------------------------------------------------
  ! SUBROUTINE print
  ! ----------------------------------------------------------------
  SUBROUTINE print(this)
    IMPLICIT NONE
    CLASS (vector), INTENT(INOUT) :: this
    CALL vector_print(this%vec)
  END SUBROUTINE print

END MODULE gridpack_vector
  
