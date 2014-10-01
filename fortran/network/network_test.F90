! ----------------------------------------------------------------
! file: network_test.f90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created August 4, 2014 by Bruce Palmer
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! PROGRAM network_test
! ----------------------------------------------------------------
#define XDIM 20
#define YDIM 20
PROGRAM network_test
  USE gridpack_network
  USE gridpack_parallel
  USE application_components
  USE, intrinsic :: iso_c_binding

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
  class(application_bus), pointer :: bus_ptr
  class(application_branch), pointer :: branch_ptr
  type(C_PTR) xc_ptr
!
!  Example data types. Replace with application-specific values
!
    integer(C_INT) int_reg
    integer(C_LONG) int_long
    real(C_FLOAT) real_s
    real(C_DOUBLE) real_d
    complex(C_FLOAT_COMPLEX) complex_s
    complex(C_DOUBLE_COMPLEX) complex_d
    logical(C_BOOL) log_reg
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
  do j = 0, ny-1 
    iy = j + iaymin
    do i = 0, nx-1
      ix = i + iaxmin
      n = (iy-1)*XDIM + ix-1
      n = 2*n   ! Provide original index that is not equal to global index
      call grid%add_bus(n+1)
!
!  Set active flag for network buses
!
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
      ncnt = ncnt + 1
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
      n2 = ly*(iaxmax-iaxmin+1) + lx + 2
      ok = grid%set_local_bus_index1(ncnt, n1)
      ok = grid%set_local_bus_index2(ncnt, n2)
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
      n2 = (ly+1)*(iaxmax-iaxmin+1) + lx + 1
      ok = grid%set_local_bus_index1(ncnt, n1)
      ok = grid%set_local_bus_index2(ncnt, n2)
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
!  Check that the number of buses and branches match the expected values
!
  n = (iaxmax-iaxmin+1)*(iaymax-iaymin+1)
  ok = .true.
  if (grid%num_buses().ne.n) then
    write(6,'(a,i4,a,i6,a,i6)') 'p[',me,'] Number of buses: ', &
      grid%num_buses(),' expected: ',n
    ok = .false.
  endif
!
  call check_ok(ok)
!
  if (me.eq.0.and.ok) then
    write(6,*)
    write(6,'(a)') 'Number of buses on each process ok'
  endif
!
  ok = .true.
  n = grid%total_buses()
  ncnt = XDIM*YDIM
  if (n.ne.ncnt) then
    write(6,'(a,i4,a,i6,a,i6)') 'p[',me,'] Total number of buses: ', &
      n,' expected: ',ncnt
    ok = .false.
  endif
!
  call check_ok(ok)
!
  if (me.eq.0.and.ok) then
    write(6,*)
    write(6,'(a)') 'Total number of buses ok'
  endif
!
  ok = .true.
  n = (iaxmax-iaxmin)*(iymax-iymin+1)+(ixmax-ixmin+1)*(iaymax-iaymin)
  if (grid%num_branches().ne.n) then
    write(6,'(a,i4,a,i6,a,i6)') 'p[',me,'] Number of branches: ', &
      grid%num_branches(),' expected: ',n
    ok = .false.
  endif
!
  call check_ok(ok)
!
  if (me.eq.0.and.ok) then
    write(6,*)
    write(6,'(a)') 'Number of branches on each process ok'
  endif
!
  ok = .true.
  n = grid%total_branches()
  ncnt = (XDIM-1)*YDIM+XDIM*(YDIM-1)
  if (ncnt.ne.n) then
    ok = .false.
  endif
!
  call check_ok(ok)
  if (me.eq.0.and.ok) then
    write(6,*)
    write(6,'(a)') 'Total number of branches process ok'
  end if
!
!  Test location of reference bus
!
  n = grid%get_reference_bus()
  if ((.not.(me.eq.0.and.n.eq.1)).and.(.not.(me.ne.0.and.n.eq.-1))) then
    write(6,'(a,i4,a,i6)') 'p[', me, '] Reference bus error: ', n
  else if (me.eq.0.and.n.eq.1) then
    write(6,*)
    write(6,'(a)') 'Reference bus ok'
  endif
!
!  Set up number of branches attached to bus
!
  nbus = grid%num_buses()
  do i = 1, nbus
    ok = grid%clear_branch_neighbors(i);
    if (.not.ok) then
      write(6,'(a,i4,a,i6)') 'p[', me, &
        '] clear_branch_neighbors failed on bus: ', i
    endif
  end do
!
!  loop over all branches
!
  nbranch = grid%num_branches()
  do i = 1, nbranch
    call grid%get_branch_endpoints(i, n1, n2)
    if (grid%get_active_bus(n1).or.grid%get_active_bus(n2)) then
      ok = grid%add_branch_neighbor(n1,i)
      if (.not.ok) then
        write(6,'(a,i4,a,i6)') 'p[', me, &
          '] add_branch_neighbor failed for bus: ', n1
      endif
      ok = grid%add_branch_neighbor(n2,i)
      if (.not.ok) then
        write(6,'(a,i4,a,i6)') 'p[', me, &
          '] add_branch_neighbor failed for bus: ', n2
      endif
    endif
  end do
!
!  Test active buses
!
  ldx = iaxmax-iaxmin+1
  ok = .true.
  do i = 1, nbus
    ix = mod(i-1,ldx);
    iy = (i-1-ix)/ldx;
    ix = (ix + iaxmin);
    iy = (iy + iaymin);
    if (ix.ge.ixmin.and.ix.le.ixmax.and.iy.ge.iymin.and.iy.le.iymax) then
      if (.not.grid%get_active_bus(i)) then
        write(6,'(a,i4,a,i6)') 'p[', me, &
          '] inactive bus error for bus: ', i
        ok = .false.
      endif
    else
      if (grid%get_active_bus(i)) then
        write(6,'(a,i4,a,i6)') 'p[', me, &
          '] active bus error for bus: ', i
        ok = .false.
      endif
    endif 
  end do
!
  call check_ok(ok)
!
  if (me.eq.0.and.ok) then
    write(6,*)
    write(6,'(a)') 'Active bus settings ok'
  endif
!
!  Test active branches
!
  ok = .true.
  do i = 1, nbranch
    call grid%get_branch_endpoints(i,n1,n2)
    ix = mod(n1-1,ldx)
    iy = (n1-1-ix)/ldx
    ix = ix + iaxmin
    iy = iy + iaymin
    if (ix.ge.ixmin.and.ix.le.ixmax.and.iy.ge.iymin.and.iy.le.iymax) then
      if (.not.grid%get_active_branch(i)) then
        write(6,'(a,i4,a,i6)') 'p[', me, &
          '] inactive branch error for branch: ', i
        ok = .false.
      endif
    else
      if (grid%get_active_branch(i)) then
        write(6,'(a,i4,a,i6)') 'p[', me, &
          '] active branch error for branch: ', i
        ok = .false.
      endif
    endif
  end do
!
  call check_ok(ok)
!
  if (me.eq.0.and.ok) then
    write(6,*)
    write(6,'(a)') 'Active branch settings ok'
  endif
!
!  Check neighbors of buses
!
  ok = .true.
  do i = 1, nbus
    if (grid%get_active_bus(i)) then
      ix = mod(i-1,ldx)
      iy = (i-1-ix)/ldx
      ix = ix + iaxmin
      iy = iy + iaymin
      n = 0
      if (ix.gt.iaxmin) n = n+1
      if (ix.lt.iaxmax) n = n+1
      if (iy.gt.iaymin) n = n+1
      if (iy.lt.iaymax) n = n+1
      nnghbr = grid%get_num_connected_branches(i)
      allocate(nghbr(nnghbr))
      allocate(nghbr_bus(nnghbr))
      call grid%get_connected_branches(i,nghbr)
      if (n.ne.nnghbr) then
        write(6,'(a,i4,a,i6)') 'p[', me, &
          '] incorrect neighbor branches on bus ', i
        ok = .false.
      endif
      do j = 1, nnghbr
        call grid%get_branch_endpoints(nghbr(j),n1,n2)
        if (n1.ne.i.and.n2.ne.i) then
          write(6,'(a,i4,a,i6)') 'p[', me, &
            '] incorrectly assigned neighbors on bus ', i
          ok = .false.
        endif
        if (n1.ne.i) then
          nghbr_bus(j) = n1
        else if (n2.ne.i) then
          nghbr_bus(j) = n2
        endif
      end do
      call grid%get_connected_buses(i,nghbr)
!
!  Double loop, not very efficient but lists are small
!
      t_ok = .false.
      do j = 1, nnghbr
        do k = 1, nnghbr
          if (nghbr(j).eq.nghbr_bus(k)) then
            t_ok = .true.
            exit
          endif
        end do
        if (t_ok) exit
      end do
      if (.not.t_ok) ok = .false.
      deallocate(nghbr)
      deallocate(nghbr_bus)
    endif
  end do
!
  call check_ok(ok)
!
  if (me.eq.0.and.ok) then
    write(6,*)
    write(6,'(a)') 'Bus neighbors are ok'
  endif
!
! Test ghost updates. Start by setting exchange data equal to original bus index
! active buses
! 
  bus_ptr => bus_cast(grid%get_bus(1))
  bsize = bus_ptr%bus_get_xc_buf_size()
  call grid%alloc_xc_bus_pointers(bsize)
  do i = 1, nbus
    bus_ptr => bus_cast(grid%get_bus(i))
    xc_ptr = bus_ptr%xc_ptr
    call grid%set_xc_bus_buffer(i,xc_ptr)
    if (grid%get_active_bus(i)) then
      bus_ptr%xc_buf%idx = grid%get_original_bus_index(i) 
    else
      bus_ptr%xc_buf%idx = -1
    endif
  end do
  branch_ptr => branch_cast(grid%get_branch(1))
  bsize = branch_ptr%branch_get_xc_buf_size()
  call grid%alloc_xc_branch_pointers(bsize)
  do i = 1, nbranch
    branch_ptr => branch_cast(grid%get_branch(i))
    xc_ptr = branch_ptr%xc_ptr
    call grid%set_xc_branch_buffer(i,xc_ptr)
    if (grid%get_active_branch(i)) then
      call grid%get_branch_endpoints(i,n1,n2)
      branch_ptr%xc_buf%idx1 = grid%get_original_bus_index(n1)
      branch_ptr%xc_buf%idx2 = grid%get_original_bus_index(n2)
    else
      branch_ptr%xc_buf%idx1 = -1
      branch_ptr%xc_buf%idx2 = -1
    endif
  end do
!
! Call update operation
!
  call grid%init_bus_update
  call grid%update_buses
#if 1
  call grid%init_branch_update
  call grid%update_branches
#endif
!
! Check updates
!
  ok = .true.
  do i = 1, nbus
    bus_ptr => bus_cast(grid%get_bus(i))
    if (bus_ptr%xc_buf%idx.ne.grid%get_original_bus_index(i)) ok = .false.
  end do
#if 1
  do i = 1, nbranch
    branch_ptr => branch_cast(grid%get_branch(i))
    call grid%get_branch_endpoints(i,n1,n2)
    n1 = grid%get_original_bus_index(n1)
    n2 = grid%get_original_bus_index(n2)
    if (branch_ptr%xc_buf%idx1.ne.n1.or.branch_ptr%xc_buf%idx2.ne.n2) ok = .false.
  end do
#endif
  call check_ok(ok) 
  if (me.eq.0.and.ok) then
    write(6,*)
    write(6,'(a)') 'Update operation is ok'
  else if (me.eq.0) then
    write(6,*)
    write(6,'(a)') 'Update operation failed'
  endif
#if 1
!
!  Test clean function
!
  call grid%clean()
  n = (ixmax-ixmin+1)*(iymax-iymin+1)
  ok = .true.
  if (grid%num_buses().ne.n) then
    write(6,'(a,i4,a,i6,a,i6)') 'p[', me, &
      '] Number of buses after clean: ', grid%num_buses(), &
      ' expected: ',n
    ok = .false.
  endif
!
  oks = n
  call MPI_Allreduce(oks, okr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
!
  if (okr.ne.XDIM*YDIM) then
    ok = .false.
  endif
!
  call check_ok(ok)
  if (me.eq.0.and.ok) then
    write(6,*)
    write(6,'(a)') 'Number of buses after clean ok'
  endif
!
  n = (iaxmax-ixmin)*(iymax-iymin+1)+(ixmax-ixmin+1)*(iaymax-iymin)
  ok = .true.
  if (grid%num_branches().ne.n) then
    write(6,'(a,i4,a,i6,a,i6)') 'p[', me, &
      '] Number of branches after clean: ', grid%num_branches(), &
      ' expected: ',n
    ok = .false.
  endif
!
  oks = n
  call MPI_Allreduce(oks, okr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
!
  if (okr.ne.(XDIM-1)*YDIM+XDIM*(YDIM-1)) then
    ok = .false.
  endif
!
  call check_ok(ok)
  if (me.eq.0.and.ok) then
    write(6,*)
    write(6,'(a)') 'Number of branches after clean ok'
  endif
!
!  Check neighbors of buses after performing clean operation
!
  ok = .true.
  nbus = grid%num_buses()
  ldx = ixmax-ixmin+1
  do i = 1, nbus
    if (grid%get_active_bus(i)) then
      ix = mod(i-1,ldx)
      iy = (i-1-ix)/ldx
      ix = ix + ixmin
      iy = iy + iymin
      n = 0
      if (ix.gt.ixmin) n = n+1
      if (ix.lt.iaxmax) n = n+1
      if (iy.gt.iymin) n = n+1
      if (iy.lt.iaymax) n = n+1
      nnghbr = grid%get_num_connected_branches(i)
      allocate(nghbr(nnghbr))
      allocate(nghbr_bus(nnghbr))
      call grid%get_connected_branches(i,nghbr)
      if (n.ne.nnghbr) then
        write(6,'(a,i4,a,i6)') 'p[', me, &
          '] incorrect neighbor branches on bus ', i
        ok = .false.
      endif
      do j = 1, nnghbr
        call grid%get_branch_endpoints(nghbr(j),n1,n2)
        if (n1.ne.i.and.n2.ne.i) then
          write(6,'(a,i4,a,i6)') 'p[', me, &
            '] incorrectly assigned neighbors on bus ', i
          ok = .false.
        endif
        if (n1.ne.i) then
          nghbr_bus(j) = n1
        else if (n2.ne.i) then
          nghbr_bus(j) = n2
        endif
      end do
      call grid%get_connected_buses(i,nghbr)
!
!  Double loop, not very efficient but lists are small
!
      t_ok = .false.
      do j = 1, nnghbr
        do k = 1, nnghbr
          if (nghbr(j).eq.nghbr_bus(k)) then
            t_ok = .true.
            exit
          endif
        end do
        if (t_ok) exit
      end do
      if (.not.t_ok) ok = .false.
      deallocate(nghbr)
      deallocate(nghbr_bus)
    endif
  end do
!
  call check_ok(ok)
!
  if (me.eq.0.and.ok) then
    write(6,*)
    write(6,'(a)') 'Bus and branches are ok after clean operation'
  endif
#endif
  call grid%destroy()
  CALL comm%finalize()
  CALL gridpack_finalize_parallel()
END PROGRAM network_test
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
