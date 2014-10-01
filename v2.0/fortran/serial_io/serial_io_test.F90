! ----------------------------------------------------------------
! file: serial_io_test.f90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created September 24, 2014 by Bruce Palmer
! Last Change: 2014-09-30 07:31:40 d3g096
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! PROGRAM serial_io_test
! ----------------------------------------------------------------
#define XDIM 10
#define YDIM 10
PROGRAM serial_io_test
  USE, intrinsic :: iso_c_binding
  USE gridpack_network
  USE gridpack_parallel
  USE gridpack_serial_io
  USE application_components
  USE application_factory

  IMPLICIT NONE
  include 'mpif.h'
  TYPE (communicator) :: comm
  INTEGER :: rank, size, p
  INTEGER(c_int) :: ma_stack = 200000, ma_heap = 200000
  TYPE(network) :: grid
  integer ierr, me, nprocs
  integer ipx, ipy, pdx, pdy
  integer ixmin, ixmax, iymin, iymax
  integer iaxmin, iaxmax, iaymin, iaymax
  integer ncnt, n, ix, iy, nx, ny, i, j, k
  integer n1, n2, lx, ly, ldx, bsize
  integer oks, okr
  logical ok, t_ok
  integer nbus, nbranch, nnghbr
  integer, allocatable :: nghbr(:), nghbr_bus(:)
  integer, allocatable :: bus_index(:)
  type(bus_serial_io) bus_io
  type(branch_serial_io) branch_io
  type(app_factory) :: io_factory
  logical bus_added
  integer one, chk, idx, jdx, isize, jsize
  double complex v
  double precision rv
  integer lo, hi, rlo, rhi
  character(32) signal
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
  call factor_grid(nprocs,XDIM,YDIM,pdx,pdy)
  if (me.eq.0) then
    write(6,*)
    write(6,'(a,i6,a,i6)') 'Processor configuration is ',pdx,' X ',pdy
  endif
  ipx = mod(me,pdx)
  ipy = (me-ipx)/pdx

  ixmin = int(dble(ipx*XDIM)/dble(pdx)) + 1
  ixmax = int(dble((ipx+1)*XDIM)/dble(pdx))
  iymin = int(dble(ipy*YDIM)/dble(pdy)) + 1
  iymax = int(dble((ipy+1)*YDIM)/dble(pdy))

  iaxmin = ixmin - 1
  if (ixmin.eq.1) iaxmin = 1
  iaxmax = ixmax + 1
  if (ixmax.eq.XDIM) iaxmax = XDIM

  iaymin = iymin - 1
  if (iymin.eq.1) iaymin = 1
  iaymax = iymax + 1
  if (iymax.eq.YDIM) iaymax = YDIM
!
! Add buses to network
!
  ncnt = 1
  nx = iaxmax - iaxmin + 1
  ny = iaymax - iaymin + 1
!
! Use bus_index aray to keep track of local index of buses
!
  allocate(bus_index(nx*ny))
  do j = 0, ny-1 
    iy = j + iaymin
    do i = 0, nx-1
      ix = i + iaxmin
      n = (iy-1)*XDIM + ix-1
      n = 2*n   ! Provide original index that is not equal to global index
      bus_added = .false.
      if (ix.eq.iaxmax.and.iy.eq.iaymax) then
        if (iaxmax.eq.ixmax.or.iaymax.eq.iymax) then
          call grid%add_bus(n+1)
          bus_added = .true.
        endif
      else if (ix.eq.iaxmin.and.iy.eq.iaymax) then
        if (iaxmin.eq.ixmin.or.iaymax.eq.iymax) then
          call grid%add_bus(n+1)
          bus_added = .true.
        endif
      else if (ix.eq.iaxmax.and.iy.eq.iaymin) then
        if (iaxmax.eq.ixmax.or.iaymin.eq.iymin) then
          call grid%add_bus(n+1)
          bus_added = .true.
        endif
      else if (ix.eq.iaxmin.and.iy.eq.iaymin) then
        if (iaxmin.eq.ixmin.or.iaymin.eq.iymin) then
          call grid%add_bus(n+1)
          bus_added = .true.
        endif
      else
        call grid%add_bus(n+1)
        bus_added = .true.
      endif
!
!  Set active flag for network buses
!
      if (bus_added) then
        if (ix.ge.ixmin.and.ix.le.ixmax.and.iy.ge.iymin.and.iy.le.iymax) then
          ok = grid%set_active_bus(ncnt,.true.)
        else
          ok = grid%set_active_bus(ncnt,.false.)
        endif
        n = n/2 + 1
        ok = grid%set_global_bus_index(ncnt,n)
        if (ix.eq.1.and.iy.eq.1) then
          call grid%set_reference_bus(ncnt)
        endif
        bus_index(j*nx+i+1) = ncnt
        ncnt = ncnt + 1
      endif
    end do
  end do
!
!  Add branches to network. Start with branches connecting buses in the
!  i-direction
!
  ncnt = 1
  nx = iaxmax-iaxmin
  ny = iymax -iymin + 1
  do j = 0, ny-1
    iy = j + iymin
    do i = 0, nx-1
      ix = i + iaxmin
      n1 = (iy-1)*XDIM+ix-1
      n1 = 2*n1
      n2 = (iy-1)*XDIM+ix
      n2 = 2*n2
      call grid%add_branch(n1+1, n2+1)
      n1 = n1/2 + 1
      n2 = n2/2 + 1
      ok = grid%set_global_bus_index1(ncnt, n1)
      ok = grid%set_global_bus_index1(ncnt, n2)
      n = (iy-1)*(XDIM-1) + ix
      ok = grid%set_global_branch_index(ncnt,n)
!
!  Figure out local indices of buses
!
      lx = ix-iaxmin
      ly = iy-iaymin
      n1 = ly*(iaxmax-iaxmin+1) + lx + 1
      n1 = bus_index(n1)
      n2 = ly*(iaxmax-iaxmin+1) + lx + 2
      n2 = bus_index(n2)
      ok = grid%set_local_bus_index1(ncnt, n1)
      ok = grid%set_local_bus_index2(ncnt, n2)
      ok = grid%add_branch_neighbor(n1, ncnt)
      ok = grid%add_branch_neighbor(n2, ncnt)
!
!  Determine which branches are locally held. Use the rule that if bus 1
!  is local, then the branch belongs to this processor
!
      if (ix.ge.ixmin.and.ix.le.ixmax.and.iy.ge.iymin.and.iy.le.iymax) then
        ok = grid%set_active_branch(ncnt,.true.)
      else
        ok = grid%set_active_branch(ncnt,.false.)
      endif
      ncnt = ncnt + 1
    end do
  end do
!
!  Add branches connecting buses in the j-direction
!
  nx = ixmax-ixmin+1
  ny = iaymax-iaymin
  do j = 0, ny-1
    iy = j + iaymin
    do i = 0, nx-1
      ix = i + ixmin
      n1 = iy*XDIM + ix
      n1 = 2*n1
      n2 = (iy+1)*XDIM + ix
      n2 = 2*n2
      call grid%add_branch(n1+1, n2+1)
      n1 = n1/2 + 1
      n2 = n2/2 + 1
      ok = grid%set_global_bus_index1(ncnt, n1)
      ok = grid%set_global_bus_index2(ncnt, n2)
      n = (iy-1)*XDIM + ix + (XDIM-1)*YDIM
      ok = grid%set_global_branch_index(ncnt, n)
!
!  Figure out local indices of buses
!
      lx = ix -iaxmin
      ly = iy -iaymin
      n1 = ly*(iaxmax-iaxmin+1) + lx + 1
      n1 = bus_index(n1)
      n2 = (ly+1)*(iaxmax-iaxmin+1) + lx + 1
      n2 = bus_index(n2)
      ok = grid%set_local_bus_index1(ncnt, n1)
      ok = grid%set_local_bus_index2(ncnt, n2)
      ok = grid%add_branch_neighbor(n1, ncnt)
      ok = grid%add_branch_neighbor(n2, ncnt)
!
!  Determine which branches are locally held. Use the rule that if bus 1
!  is local, then the branch belongs to this processor
!
      if (ix.ge.ixmin.and.ix.le.ixmax.and.iy.ge.iymin.and.iy.le.iymax) then
        ok = grid%set_active_branch(ncnt,.true.)
      else
        ok = grid%set_active_branch(ncnt,.false.)
      endif
      ncnt = ncnt + 1
    end do
  end do
  deallocate(bus_index)
!
!  Set up remaining properties of the network and network components so that
!  matrix-vector interface is ready to go
!
  call io_factory%create(grid)
  call io_factory%set_components()
!
  signal = ''
  call bus_io%bus_create(128,grid)
  call bus_io%bus_header('')
  call bus_io%bus_header('  Bus Properties');
  call bus_io%bus_header('')
  call bus_io%bus_header('       Original   Global');
  call bus_io%bus_write(signal)
!
  call branch_io%branch_create(128,grid)
  call branch_io%branch_header('')
  call branch_io%branch_header('  Branch Properties');
  call branch_io%branch_header('')
  call branch_io%branch_header('         Original1  Original2  Global1 Global2');
  call branch_io%branch_write(signal)
!
  call bus_io%bus_destroy()
  call branch_io%branch_destroy()
  call io_factory%destroy()
  call grid%destroy()
  CALL comm%finalize()
  CALL gridpack_finalize_parallel()
END PROGRAM serial_io_test
!
!  Subroutine to factor grid
!
subroutine factor_grid(nproc, xsize, ysize, pdx, pdy)
  implicit none
  integer nproc, xsize, ysize
  integer pdx, pdy
  integer i,j,it,ip,ifac,pmax,prime(0:999),chk
  integer idx, idy
  integer xtmp, ytmp
  integer fac(0:999)

  ip = nproc
!
!  Factor nproc completely, First, find all prime numbers less than or equal to
!  nproc
!
  pmax = 0
  do i = 2, nproc
    chk = 0
    j = 0
    do j = 0, pmax-1
      if (mod(i,prime(j)).eq.0) then
        chk = 1
        exit
      endif
    end do
    if (chk.eq.0) then
      prime(pmax) = i
      pmax = pmax + 1
    endif
  end do
!
!  Find all prime factors of nproc
!
  ifac = 0
  do i = 0, pmax-1
    do while (mod(ip,prime(i)).eq.0.and.ip.gt.1)
      ifac = ifac + 1
      fac(ifac) = prime(i) 
      ip = ip / prime(i)
    end do
  end do
!
!  Find three factors of nproc such that the resulting three dimensions of the
!  simulation cell on each processor are as close as possible to being the same
!  size
!
  xtmp = xsize
  ytmp = ysize
  idx = 1
  idy = 1
  do i = ifac, 1, -1
    if (xtmp.ge.ytmp) then
      idx = fac(i)*idx
      xtmp = xtmp/fac(i)
    else
      idy = fac(i)*idy
      ytmp = ytmp/fac(i)
    endif
  end do
  pdx = idx
  pdy = idy
!
  return
end subroutine factor_grid
!
!  Subroutine to check if something is true on all processors
!
subroutine check_ok(flag)
  implicit none
include 'mpif.h'
  logical, intent(inout) :: flag
  integer oks, okr, ierr
  if (flag) then
    oks = 1
  else
    oks = 0
  endif
  call MPI_Allreduce(oks, okr, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD, ierr)
  if (okr.eq.1) then
    flag = .true.
  else
    flag = .false.
  endif
end subroutine check_ok
