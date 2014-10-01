! ----------------------------------------------------------------
! file: parallel.f90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created May 14, 2014 by William A. Perkins
! Last Change: 2014-05-15 14:07:25 d3g096
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! MODULE gridpack_parallel
! ----------------------------------------------------------------
MODULE gridpack_parallel
  
  USE iso_c_binding, ONLY: c_int
  USE gridpack_communicator

  IMPLICIT NONE

  PUBLIC
  INTERFACE
     SUBROUTINE gridpack_initialize_parallel(stack, heap) BIND(c)
       USE iso_c_binding, ONLY: c_int
       IMPLICIT NONE
       INTEGER(c_int), INTENT(IN) :: heap, stack
     END SUBROUTINE gridpack_initialize_parallel
     SUBROUTINE gridpack_finalize_parallel() BIND(c)
     END SUBROUTINE gridpack_finalize_parallel
  END INTERFACE

END MODULE gridpack_parallel
