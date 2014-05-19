! ----------------------------------------------------------------
! file: hello_math.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created May 14, 2014 by William A. Perkins
! Last Change: 2014-05-15 14:05:56 d3g096
! ----------------------------------------------------------------
PROGRAM hello_math

  USE gridpack_parallel
  USE gridpack_math

  IMPLICIT NONE

  TYPE (communicator) :: comm
  INTEGER(c_int) :: ma_stack = 200000, ma_heap = 200000

  CALL gridpack_initialize_parallel(ma_stack, ma_heap)
  CALL gridpack_initialize_math()
  CALL comm%initialize()

  IF (gridpack_math_initialized()) THEN
     WRITE (*, *) 'GridPACK math library initialized'
  ELSE 
     WRITE (*, *) 'GridPACK math library not initialized'
  END IF

  CALL comm%finalize()
  CALL gridpack_finalize_math()
  CALL gridpack_finalize_parallel()

END PROGRAM hello_math
