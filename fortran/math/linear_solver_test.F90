! ----------------------------------------------------------------
! file: linear_solver_test.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created May 29, 2014 by William A. Perkins
! Last Change: 2014-08-22 08:31:33 d3g096
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! SUBROUTINE assemble
! ----------------------------------------------------------------
SUBROUTINE assemble(imax, jmax, a, b)

  USE iso_c_binding
  USE gridpack_parallel
  USE gridpack_math

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: imax, jmax
  CLASS (matrix), INTENT(INOUT) :: a
  CLASS (vector), INTENT(INOUT) :: b

  DOUBLE PRECISION, PARAMETER :: k = 1000 ! conductivity, W/m/K 
  DOUBLE PRECISION, PARAMETER :: t = 0.01 ! plate thickness, m
  DOUBLE PRECISION, PARAMETER :: W = 0.3  ! plate width, m
  DOUBLE PRECISION, PARAMETER :: H = 0.4  ! plate height, m 
  DOUBLE PRECISION :: dx, dy

  INTEGER :: ilo, ihi, i, j
  COMPLEX(c_double_complex) :: ap, aw, ae, as, an, bp
  INTEGER :: iP, iN, iS, iE, iW

  dx = W/REAL(imax)
  dy = H/REAL(jmax)
  

  CALL b%local_index_range(ilo, ihi)

  DO i = 0, imax-1
     DO j = 0, jmax-1
        iP = i*jmax + j
        IF (ilo .LE. iP .AND. iP .LT. ihi) THEN
           iE = (i+1)*jmax + j
           iW = (i-1)*jmax + j
           iN = i*jmax + (j+1)
           iS = i*jmax + (j-1)

           bp = 0.0
           ap = 0.0
           
           IF (j .EQ. 0) THEN      ! insulated south boundary
              as = 0.0
              bp = bp + 0.0
              ap = ap - 0.0
           ELSE 
              as = (k/dx)*(dx*t)
           END IF
           
           
           IF (j .EQ. jmax - 1) THEN ! constant temperature (100C) north boundary
              an = 0.0
              bp = bp + 2*k/dy*(dy*t)*100.0
              ap = ap - (-2*k/dy*(dy*t))
           ELSE 
              an = (k/dx)*(dx*t)
           END IF
           
           IF (i .EQ. 0) THEN    ! constant flux (500kw/m2) west boundary
              aw = 0.0
              bp = bp + 500000.0*dy*t
              ap = ap - 0.0
           ELSE
              aw = (k/dx)*(dx*t)
           END IF
           
           IF (i .EQ. imax - 1) THEN ! insulated east boundary 
              ae = 0.0
              bp = bp + 0.0 
              ap = ap - 0.0
           ELSE 
              ae = (k/dx)*(dx*t)
           END IF
           
           ap = ap + as + an + aw + ae
           
           CALL a%set_element(iP, iP, ap)
           
           IF (an .NE. 0.0) CALL a%set_element(iP, iN, -an)
           IF (as .NE. 0.0) CALL a%set_element(iP, iS, -as)
           IF (ae .NE. 0.0) CALL a%set_element(iP, iE, -ae)
           IF (aw .NE. 0.0) CALL a%set_element(iP, iW, -aw)
           CALL b%set_element(iP, bp)
        END IF
     END DO
  END DO

END SUBROUTINE assemble


! ----------------------------------------------------------------
! PROGRAM linear_solver_test
! ----------------------------------------------------------------
  
PROGRAM linear_solver_test

  USE iso_c_binding
  USE gridpack_parallel
  USE gridpack_configuration
  USE gridpack_math

  IMPLICIT NONE

  TYPE (communicator) :: comm
  INTEGER(c_int) :: ma_stack = 200000, ma_heap = 200000

  TYPE (cursor) :: conf

  TYPE (matrix) :: a
  TYPE (vector) :: b, x
  TYPE (linear_solver) :: solver

  INTEGER :: me, nproc
  INTEGER :: imax, jmax, global_size, local_size

  ! This is needed w/ GNU Fortran 4.8 for some reason?
  INTERFACE
     SUBROUTINE assemble(imax, jmax, a, b)
       USE iso_c_binding
       USE gridpack_parallel
       USE gridpack_math
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: imax, jmax
       CLASS (matrix), INTENT(INOUT) :: a
       CLASS (vector), INTENT(INOUT) :: b
     END SUBROUTINE assemble
  END INTERFACE

  CALL gridpack_initialize_parallel(ma_stack, ma_heap)
  CALL gridpack_initialize_math()
  CALL comm%initialize()

  IF (conf%open(comm, "gridpack.xml")) THEN
     CALL conf%set_path("GridPACK.MathTests")
  ELSE
     WRITE(*,*) "ERROR: cannot open configuration"
     STOP
  ENDIF

  me = comm%rank()
  nproc = comm%size()

  imax = 3*nproc
  jmax = 4*nproc
  global_size = imax*jmax
  local_size = global_size/nproc

  CALL a%initialize(comm, local_size, local_size, 13)
  CALL b%initialize(comm, local_size)
  CALL x%initialize(comm, local_size)

  CALL assemble(imax, jmax, a, b)
  CALL a%ready()
  CALL b%ready()
  CALL x%zero()

  CALL solver%initialize(a, conf)
  CALL solver%solve(b, x)

  CALL x%print()
  
  CALL solver%finalize()
  CALL a%finalize()
  CALL b%finalize()
  CALL x%finalize()
  CALL comm%finalize()
  CALL gridpack_finalize_math()
  CALL gridpack_finalize_parallel()

  
END PROGRAM linear_solver_test
