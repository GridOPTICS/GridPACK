! ----------------------------------------------------------------
! file: hello.f90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created May 14, 2014 by William A. Perkins
! Last Change: 2014-05-15 14:07:10 d3g096
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! PROGRAM hello
! ----------------------------------------------------------------
PROGRAM hello

  USE gridpack_parallel
  
  IMPLICIT NONE

  TYPE (communicator) :: comm
  INTEGER :: ierror, rank, size, p
  INTEGER(c_int) :: ma_stack = 200000, ma_heap = 200000

  CALL gridpack_initialize_parallel(ma_stack, ma_heap)
  CALL comm%initialize()

  rank = comm%rank()
  size = comm%size()
  DO p = 0, size
     IF (rank .EQ. p) THEN 
        print*, 'node', rank, 'of', size, ': Hello world'
     END IF
     CALL comm%barrier()
  END DO
  CALL comm%finalize()
  CALL gridpack_finalize_parallel()
END PROGRAM hello
