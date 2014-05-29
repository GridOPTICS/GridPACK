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
! Last Change: 2014-05-29 14:03:59 d3g096
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! MODULE gridpack_linear_solver
! ----------------------------------------------------------------
MODULE gridpack_linear_solver

  USE iso_c_binding
  USE gridpack_vector
  USE gridpack_matrix

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: linear_solver
     TYPE (c_ptr) :: solvr
   CONTAINS
     PROCEDURE :: initialize
     PROCEDURE :: finalize
     PROCEDURE :: solve
  END type linear_solver

  INTERFACE
     SUBROUTINE linear_solver_initialize(solvr, a) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       TYPE (c_ptr), INTENT(INOUT) :: solvr
       TYPE (c_ptr), VALUE, INTENT(IN) :: a
     END SUBROUTINE linear_solver_initialize

     SUBROUTINE linear_solver_finalize(solvr) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       TYPE (c_ptr), INTENT(INOUT) :: solvr
     END SUBROUTINE linear_solver_finalize

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
  SUBROUTINE initialize(this, mat)
    IMPLICIT NONE
    CLASS (linear_solver), INTENT(INOUT) :: this
    CLASS (matrix), INTENT(IN) :: mat
    CALL linear_solver_initialize(this%solvr, mat%mat)
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
  ! SUBROUTINE linear_solver_solve
  ! ----------------------------------------------------------------
  SUBROUTINE solve(this, b, x)
    IMPLICIT NONE
    CLASS (linear_solver), INTENT(IN) :: this
    CLASS (vector), INTENT(IN) :: b, x
    CALL linear_solver_solve(this%solvr, b%vec, x%vec)
  END SUBROUTINE solve


END MODULE gridpack_linear_solver
  
