! ----------------------------------------------------------------
! file: parser_test.f90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created September 23, 2014 by Bruce Palmer
! Last Change: 2014-09-23 14:07:10 d3g293
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! PROGRAM parser
! ----------------------------------------------------------------
PROGRAM parser_test
  USE, intrinsic :: iso_c_binding
  USE gridpack_network
  USE gridpack_parser
  USE gridpack_parallel
  USE gridpack_data_collection
  USE application_components
  USE application_factory
!
  implicit none
#include "mpif.h"
  INTEGER(c_int) :: ma_stack = 200000, ma_heap = 200000
  TYPE (communicator) :: comm
  TYPE(network) :: grid
  TYPE(data_collection) :: data_ptr
  TYPE(pti23_parser) parser
  TYPE(application_bus), pointer :: bus_ptr
  TYPE(application_branch), pointer :: branch_ptr
  integer nbus, nbranch, i, idx
  integer schk, rchk, ierr
  double precision rval
  character(len=32) key, busname
  logical ok
  integer me, nprocs
!
! Initialize gridpack and create world communicator
!
  CALL gridpack_initialize_parallel(ma_stack, ma_heap)
  CALL comm%initialize()
  me = comm%rank()
  nprocs = comm%size()
!
!  Initialize network
!
  call grid%create(comm)
!
  call parser%create(grid)
  call parser%parse("parser_data.raw")
  call grid%partition()
!
! Check to see if contents of data_collection objects are correct
!
  schk = 0
  nbus = grid%num_buses()
  do i = 1, nbus
    call data_ptr%bind(grid%get_bus_data(i))
    key = 'LOAD_PL'
    ok = data_ptr%get_double_value(key,rval)
    idx = grid%get_global_bus_index(i) - 1
    if (rval.ne.dble(mod(idx,100))) then
      write(6,'(a,f6.1,a,f6.1)') 'Mismatch for LOAD_PL expected: ', &
         dble(mod(idx,100)),' actual: ',rval
      schk = 1
    endif
    key = 'LOAD_QL'
    ok = data_ptr%get_double_value(key,rval)
    if (rval.ne.dble(mod(idx+1,100))) then
      write(6,'(a,f6.1,a,f6.1)') 'Mismatch for LOAD_QL expected: ', &
         dble(mod(idx+1,100)),' actual: ',rval
      schk = 1
    endif
    key = 'BUS_NAME'
    ok = data_ptr%get_string_value(key,busname)
    write(6,'(a,a)') 'Bus name: ',trim(busname)
  end do
!
  nbranch = grid%num_branches()
  do i = 1, nbranch
    call data_ptr%bind(grid%get_branch_data(i))
    key = 'BRANCH_RATING_A'
    ok = data_ptr%get_double_indexed_value(key,rval,1)
    idx = grid%get_global_branch_index(i) - 1
    if (rval.ne.dble(mod(idx,100))) then
      write(6,'(a,f6.1,a,f6.1)') 'Mismatch for BRANCH_RATING expected: ', &
         dble(mod(idx+1,100)),' actual: ',rval
      schk = 1
    endif
  end do
  call mpi_allreduce(schk,rchk,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (rchk.eq.0.and.me.eq.0) then
    write(6,*)
    write(6,'(a)') 'Parsing of test configuration is ok'
  else if (me.eq.0) then
    write(6,*)
    write(6,'(a)') 'Error in parsing of test configuration'
  endif
!
  call grid%destroy()
  CALL comm%finalize()
  CALL gridpack_finalize_parallel()
  stop
END PROGRAM parser_test
