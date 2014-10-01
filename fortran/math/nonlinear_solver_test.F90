! ----------------------------------------------------------------
! file: nonlinear_solver_test.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created September 15, 2014 by William A. Perkins
! Last Change: 2014-10-01 08:35:03 d3g096
! ----------------------------------------------------------------

MODULE tiny
  USE gridpack_nonlinear_solver
  IMPLICIT NONE
  TYPE, EXTENDS(builder) :: tiny_builder
   CONTAINS
     PROCEDURE :: jacobian => tiny_jacobian
     PROCEDURE :: function => tiny_function
  END type tiny_builder
CONTAINS
  ! ----------------------------------------------------------------
  ! SUBROUTINE tiny_jacobian
  ! ----------------------------------------------------------------
  SUBROUTINE tiny_jacobian(this, x, J)
    USE iso_c_binding, ONLY : c_double_complex
    USE gridpack_matrix
    USE gridpack_vector
    IMPLICIT NONE
    CLASS (tiny_builder), INTENT(IN) :: this
    TYPE (vector), INTENT(IN) :: x
    TYPE (matrix), INTENT(INOUT) :: J
    COMPLEX(c_double_complex) :: x0, x1, y
    INTEGER :: lo, hi
    
    ! assume the same for both x and J rows
    CALL x%local_index_range(lo, hi)
    
    CALL x%get_element(lo+0, x0)
    CALL x%get_element(lo+1, x1)
    CALL J%set_element(lo+0, lo+0, 2.0*x0-2.0)
    y = -1.0
    CALL J%set_element(lo+0, lo+1, y)
    CALL J%set_element(lo+1, lo+0, 2.0*x0)
    CALL J%set_element(lo+1, lo+1, 8.0*x1)
    Call J%ready()
  END SUBROUTINE tiny_jacobian

  ! ----------------------------------------------------------------
  ! SUBROUTINE tiny_function
  ! ----------------------------------------------------------------
  SUBROUTINE tiny_function(this, x, F)
    USE iso_c_binding, ONLY : c_double_complex
    USE gridpack_vector
    IMPLICIT NONE
    CLASS (tiny_builder), INTENT(IN) :: this
    TYPE (vector), INTENT(IN) :: x
    TYPE (vector), INTENT(INOUT) :: F
    COMPLEX(c_double_complex) :: x0, x1
    INTEGER :: lo, hi

    CALL x%local_index_range(lo, hi)

    CALL x%get_element(lo+0, x0)
    CALL x%get_element(lo+1, x1)
    
    CALL F%set_element(lo+0, x0*x0 - 2.0*x0 - x1 + 0.5)
    CALL F%set_element(lo+1, x0*x0 + 4.0*x1*x1 - 4.0)
    CALL F%ready()
  END SUBROUTINE tiny_function
END MODULE tiny

PROGRAM nonlinear_solver_test
  USE gridpack_parallel
  USE gridpack_math
  USE tiny
  IMPLICIT NONE

  TYPE (communicator) :: comm
  INTEGER(c_int) :: ma_stack = 200000, ma_heap = 200000
  TYPE (cursor) :: conf
  CLASS (builder), POINTER :: bldr
  CLASS (nonlinear_solver), POINTER :: solver
  TYPE (vector) :: x
  COMPLEX(c_double_complex) :: myx
  INTEGER :: lo, hi

  CALL gridpack_initialize_parallel(ma_stack, ma_heap)
  CALL gridpack_initialize_math()
  CALL comm%initialize()

  IF (conf%open(comm, "gridpack.xml")) THEN
     CALL conf%set_path("GridPACK.MathTests")
  ELSE
     WRITE(*,*) "ERROR: cannot open configuration"
     STOP
  ENDIF

  ALLOCATE(tiny_builder::bldr)
  CALL x%initialize(comm, 2)

#if NEWTON
  ALLOCATE(newton_raphson_solver::solver)
#else
  ALLOCATE(nonlinear_solver::solver)
#endif

  CALL solver%initialize(comm, conf, 2, bldr)
  
  CALL x%local_index_range(lo, hi)
  myx = 2.00
  CALL x%set_element(lo+0, myx)
  myx = 0.25
  CALL x%set_element(lo+1, myx)
  CALL x%ready()
  CALL solver%solve(x)
  CALL x%print()

  CALL solver%finalize()
  DEALLOCATE(solver)
  CALL x%finalize()
  DEALLOCATE(bldr)
  CALL conf%finalize()
  CALL comm%finalize()
  CALL gridpack_finalize_math()
  CALL gridpack_finalize_parallel()

END PROGRAM nonlinear_solver_test

