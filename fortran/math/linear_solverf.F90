! ----------------------------------------------------------------
! file: linear_solverf.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created May 29, 2014 by William A. Perkins
! Last Change: 2014-08-22 07:46:14 d3g096
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! MODULE gridpack_linear_solver
! ----------------------------------------------------------------
MODULE gridpack_linear_solver

  USE iso_c_binding
  USE gridpack_configuration
  USE gridpack_vector
  USE gridpack_matrix

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: linear_solver
     TYPE (c_ptr) :: solvr
   CONTAINS
     PROCEDURE :: initialize
     PROCEDURE :: finalize      ! should be FINAL
     PROCEDURE :: solve
  END type linear_solver

  INTERFACE
     SUBROUTINE linear_solver_initialize(solvr, a, conf) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       TYPE (c_ptr), INTENT(INOUT) :: solvr
       TYPE (c_ptr), VALUE, INTENT(IN) :: a
       TYPE (c_ptr), VALUE, INTENT(IN) :: conf
     END SUBROUTINE linear_solver_initialize

     SUBROUTINE linear_solver_finalize(solvr) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       TYPE (c_ptr), INTENT(INOUT) :: solvr
     END SUBROUTINE linear_solver_finalize

     SUBROUTINE linear_solver_set_tolerance(solvr, tol) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_double
       TYPE (c_ptr), VALUE, INTENT(IN) :: solvr
       REAL(c_double), VALUE, INTENT(IN) :: tol
     END SUBROUTINE linear_solver_set_tolerance

     SUBROUTINE linear_solver_set_iterations(solvr, iter) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_int
       TYPE (c_ptr), VALUE, INTENT(IN) :: solvr
       INTEGER(c_int), VALUE, INTENT(IN) :: iter
     END SUBROUTINE linear_solver_set_iterations

     SUBROUTINE linear_solver_solve(solvr, b, x) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       TYPE (c_ptr), VALUE, INTENT(IN) :: solvr
       TYPE (c_ptr), VALUE, INTENT(IN) :: b, x
     END SUBROUTINE linear_solver_solve
  END INTERFACE
CONTAINS

  ! ----------------------------------------------------------------
  ! SUBROUTINE initialize
  ! ----------------------------------------------------------------
  SUBROUTINE initialize(this, mat, conf)
    IMPLICIT NONE
    CLASS (linear_solver), INTENT(INOUT) :: this
    CLASS (matrix), INTENT(IN) :: mat
    CLASS (cursor), INTENT(IN), OPTIONAL :: conf
    TYPE(c_ptr) :: confptr
    confptr = C_NULL_PTR
    IF (PRESENT(conf)) THEN
       confptr = conf%impl
    END IF
    CALL linear_solver_initialize(this%solvr, mat%mat, confptr)
  END SUBROUTINE initialize

  ! ----------------------------------------------------------------
  ! SUBROUTINE finalize
  ! ----------------------------------------------------------------
  SUBROUTINE finalize(this)
    IMPLICIT NONE
    CLASS (linear_solver), INTENT(INOUT) :: this
    CALL linear_solver_finalize(this%solvr)
  END SUBROUTINE finalize

  ! ----------------------------------------------------------------
  ! SUBROUTINE tolerance
  ! ----------------------------------------------------------------
  SUBROUTINE tolerance(this, tol)
    IMPLICIT NONE
    CLASS (linear_solver), INTENT(INOUT) :: this
    DOUBLE PRECISION, INTENT(IN) :: tol
    REAL(c_double) :: ctol
    ctol = tol
    CALL linear_solver_set_tolerance(this%solvr, ctol)
  END SUBROUTINE tolerance


  ! ----------------------------------------------------------------
  ! SUBROUTINE iterations
  ! ----------------------------------------------------------------
  SUBROUTINE iterations(this, maxiter)
    IMPLICIT NONE
    CLASS (linear_solver), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: maxiter
    INTEGER(c_int) :: cmaxiter
    cmaxiter = maxiter
    CALL linear_solver_set_iterations(this%solvr, cmaxiter)
  END SUBROUTINE iterations


  ! ----------------------------------------------------------------
  ! SUBROUTINE solve
  ! ----------------------------------------------------------------
  SUBROUTINE solve(this, b, x)
    IMPLICIT NONE
    CLASS (linear_solver), INTENT(IN) :: this
    CLASS (vector), INTENT(IN) :: b, x
    CALL linear_solver_solve(this%solvr, b%vec, x%vec)
  END SUBROUTINE solve


END MODULE gridpack_linear_solver
  
