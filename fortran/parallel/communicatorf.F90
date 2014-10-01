! ----------------------------------------------------------------
! file: communicatorf.f90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created May 14, 2014 by William A. Perkins
! Last Change: 2014-05-15 15:23:32 d3g096
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! MODULE gridpack_communicator
! ----------------------------------------------------------------
MODULE gridpack_communicator

  USE iso_c_binding

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: communicator
     TYPE (c_ptr) :: comm
   CONTAINS
     PROCEDURE :: initialize 
     PROCEDURE :: finalize 
     PROCEDURE :: size
     PROCEDURE :: rank
     PROCEDURE :: barrier
     PROCEDURE :: sync
  END TYPE communicator

  INTERFACE 
     SUBROUTINE communicator_initialize(comm) BIND(c)
       USE iso_c_binding, ONLY : c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(INOUT) :: comm
     END SUBROUTINE communicator_initialize
     SUBROUTINE communicator_finalize(comm) BIND(c)
       USE iso_c_binding, ONLY : c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(INOUT) :: comm
     END SUBROUTINE communicator_finalize
     INTEGER (c_int) FUNCTION communicator_size(comm) BIND(c)
       USE iso_c_binding, ONLY : c_ptr, c_int
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: comm
     END FUNCTION communicator_size
     INTEGER (c_int) FUNCTION communicator_rank(comm) BIND(c)
       USE iso_c_binding, ONLY : c_ptr, c_int
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: comm
     END FUNCTION communicator_rank
     INTEGER (c_int) FUNCTION communicator_world_rank(comm) BIND(c)
       USE iso_c_binding, ONLY : c_ptr, c_int
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: comm
     END FUNCTION communicator_world_rank
     SUBROUTINE communicator_barrier(comm) BIND(c)
       USE iso_c_binding, ONLY : c_ptr, c_int
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: comm
     END SUBROUTINE communicator_barrier
     SUBROUTINE communicator_sync(comm) BIND(c)
       USE iso_c_binding, ONLY : c_ptr, c_int
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: comm
     END SUBROUTINE communicator_sync
  END INTERFACE

CONTAINS

  ! ----------------------------------------------------------------
  ! SUBROUTINE initialize
  ! ----------------------------------------------------------------
  SUBROUTINE initialize(this)
    IMPLICIT NONE
    CLASS (communicator), INTENT(INOUT) :: this
    CALL communicator_initialize(this%comm)
  END SUBROUTINE initialize

  ! ----------------------------------------------------------------
  ! SUBROUTINE finalize
  ! ----------------------------------------------------------------
  SUBROUTINE finalize(this)
    IMPLICIT NONE
    CLASS (communicator), INTENT(INOUT) :: this
    CALL communicator_finalize(this%comm)
  END SUBROUTINE finalize

  ! ----------------------------------------------------------------
  ! INTEGER FUNCTION size
  ! ----------------------------------------------------------------
  INTEGER FUNCTION size(this)
    IMPLICIT NONE
    CLASS (communicator), INTENT(IN) :: this
    size = communicator_size(this%comm)
  END FUNCTION size

  ! ----------------------------------------------------------------
  ! INTEGER FUNCTION rank
  ! ----------------------------------------------------------------
  INTEGER FUNCTION rank(this)
    IMPLICIT NONE
    CLASS (communicator), INTENT(IN) :: this
    rank = communicator_rank(this%comm)
  END FUNCTION rank

  ! ----------------------------------------------------------------
  ! SUBROUTINE barrier
  ! ----------------------------------------------------------------
  SUBROUTINE barrier(this)
    IMPLICIT NONE
    CLASS (communicator), INTENT(IN) :: this
    CALL communicator_barrier(this%comm)
  END SUBROUTINE barrier

  ! ----------------------------------------------------------------
  ! SUBROUTINE sync
  ! ----------------------------------------------------------------
  SUBROUTINE sync(this)
    IMPLICIT NONE
    CLASS (communicator), INTENT(IN) :: this
    CALL communicator_sync(this%comm)
  END SUBROUTINE sync


END MODULE gridpack_communicator
