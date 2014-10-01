! ----------------------------------------------------------------
! file: mpitest.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created August  6, 2014 by William A. Perkins
! Last Change: 2014-08-06 12:04:33 d3g096
! ----------------------------------------------------------------
PROGRAM mpitest
  USE MPI

  IMPLICIT NONE

  INTEGER :: ierr
  INTEGER :: nproc
  INTEGER :: me

  CALL MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)

  WRITE ( *, '(A,I8,A)' ) '  Process ', me, ' says "Hello, world!"'

  call MPI_Finalize (ierr)

END PROGRAM mpitest
