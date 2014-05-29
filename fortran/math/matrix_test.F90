! ----------------------------------------------------------------
! file: matrix_test.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created May 29, 2014 by William A. Perkins
! Last Change: 2014-05-29 10:38:18 d3g096
! ----------------------------------------------------------------
PROGRAM matrix_test

  USE gridpack_parallel
  USE gridpack_math

  IMPLICIT NONE

  TYPE (communicator) :: comm
  INTEGER(c_int) :: ma_stack = 200000, ma_heap = 200000

  TYPE (matrix) :: mat, mat2
  INTEGER :: me, lo, hi, i, j, rows, cols
  COMPLEX(c_double_complex) :: myx
  COMPLEX(c_double_complex) :: norm1, norm2, norminf

  CALL gridpack_initialize_parallel(ma_stack, ma_heap)
  CALL gridpack_initialize_math()
  CALL comm%initialize()

  me = comm%rank()

  CALL mat%initialize(comm, 5, 5, 5)
  rows = mat%rows()
  cols = mat%cols()

  CALL mat%local_row_range(lo, hi)
  WRITE(*, '(I2, ": I own ", I2, " of ", I2, " rows: ", I2, "-", I2)') &
       &me, mat%local_rows(), rows, lo, hi

  DO i = lo, hi-1
     DO j = MAX(i-2, 0), MIN(i+2, hi-1)
        myx = me
        CALL mat%set_element(i, j, myx)
     END DO
  END DO
  CALL mat%ready()
  CALL mat%print()

  mat2 = mat%clone()
  CALL mat2%add(mat)
  CALL mat2%print()

  myx = cmplx(3.0, 3.0)
  CALL mat%scale(myx)
  CALL mat%print()

  CALL mat%finalize()
  CALL comm%finalize()
  CALL gridpack_finalize_math()
  CALL gridpack_finalize_parallel()


END PROGRAM matrix_test
