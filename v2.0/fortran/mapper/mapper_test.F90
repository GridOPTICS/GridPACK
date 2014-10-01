! ----------------------------------------------------------------
! file: mapper_test.f90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created September 15, 2014 by Bruce Palmer
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! PROGRAM mapper_test
! ----------------------------------------------------------------
#define XDIM 10
#define YDIM 10
PROGRAM mapper_test
  USE, intrinsic :: iso_c_binding
  USE gridpack_network
  USE gridpack_full_matrix_map
  USE gridpack_bus_vector_map
  USE gridpack_gen_matrix_map
  USE gridpack_gen_vector_map
  USE gridpack_parallel
  USE gridpack_math
  USE gridpack_vector
  USE gridpack_matrix
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
  type(application_bus), pointer :: bus_ptr
  type(application_branch), pointer :: branch_ptr
  type(app_factory) :: map_factory
  type(matrix) :: mmatrix
  type(vector) :: vvector
  type(full_matrix_map) :: fm_map
  type(bus_vector_map) :: bv_map
  type(C_PTR) xc_ptr
  logical bus_added
  integer one, chk, idx, jdx, isize, jsize
  double complex v
  double precision rv
  integer lo, hi, rlo, rhi
!
! Initialize gridpack and create world communicator
!
  CALL gridpack_initialize_parallel(ma_stack, ma_heap)
  CALL comm%initialize()
  CALL gridpack_initialize_math()
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
  call map_factory%create(grid)
  call map_factory%set_components()
!
  if (me.eq.0) then
    write(6,*)
    write(6,'(a)') 'Testing full_matrix_map'
  endif
  call fm_map%create(grid)
  mmatrix = fm_map%map_to_matrix()
  one = 1
  chk = 0
  nbus = grid%num_buses()
  nbranch = grid%num_branches()
  do i = 1, nbus
    if (grid%get_active_bus(i)) then
      bus_ptr => bus_cast(grid%get_bus(i))
      if (bus_ptr%bus_matrix_diag_size(isize,jsize).and. &
          isize.gt.0.and.jsize.gt.0) then
        bus_ptr => bus_cast(grid%get_bus(i))
        call bus_ptr%bus_get_mat_vec_index(idx)
        idx = idx - 1
        call mmatrix%get_element(idx, idx, v)
        rv = real(v)
        if (rv.ne.-4.0d00) then
          write(6,'(a,i5,a,i5,a,i5,a,e12.4)') 'p[',me, &
           '] Diagonal matrix error i: ',idx,' j: ',idx,' v: ',rv
          chk = 1
        endif
      endif
    endif
  end do
!
! Get min and max row indices
!
  rhi = 0
  rlo = XDIM*YDIM
  do i = 1, nbus
    if (grid%get_active_bus(i)) then
      bus_ptr => bus_cast(grid%get_bus(i))
      call bus_ptr%bus_get_mat_vec_index(idx)
      if (rhi.lt.idx) rhi = idx
      if (rlo.gt.idx) rhi = idx
    endif
  end do
  do i = 1, nbranch
    branch_ptr => branch_cast(grid%get_branch(i))
    if (branch_ptr%branch_matrix_forward_size(isize,jsize).and. &
        isize.gt.0.and.jsize.gt.0) then
      call branch_ptr%branch_get_mat_vec_indices(idx,jdx)
      idx = idx - 1
      jdx = jdx - 1
      if (idx.ge.rlo-1.and.idx.le.rhi-1) then
        call mmatrix%get_element(idx,jdx,v);
        rv = real(v)
        if (rv.ne.1.0d00) then
          write(6,'(a,i5,a,i5,a,i5,a,e12.4)') 'p[',me, &
            '] Forward matrix error i: ',idx,' j: ',jdx,' v: ',rv
          chk = 1
        endif
      endif
      if (jdx.ge.rlo-1.and.idx.le.rhi-1) then
        call mmatrix%get_element(jdx,idx,v);
        rv = real(v)
        if (rv.ne.1.0d00) then
          write(6,'(a,i5,a,i5,a,i5,a,e12.4)') 'p[',me, &
            '] Reverse matrix error i: ',jdx,' j: ',idx,' v: ',rv
          chk = 1
        endif
      endif
    endif
  end do
  call ga_igop(1,chk,1,"+")
  if (me.eq.0) then
    write(6,*)
    if (chk.eq.0) then
      write(6,'(a)') 'Matrix elements are ok'
    else
      write(6,'(a)') 'Error found in matrix elements'
    endif
  endif
  if (me.eq.0) then
    write(6,*)
    write(6,'(a)') 'Testing bus_vector_map'
  endif
  call bv_map%create(grid)
  vvector = bv_map%map_to_vector()
!
! Check to see if vector has correct values
!
  chk = 0
  do i = 1, nbus
    if (grid%get_active_bus(i)) then
      bus_ptr => bus_cast(grid%get_bus(i))
      if (bus_ptr%bus_vector_size(isize)) then
        call bus_ptr%bus_get_mat_vec_index(idx)
        idx = idx - 1
        call vvector%get_element(idx,v)
        rv = real(v)
        if (rv.ne.dble(idx+1)) then
          write(6,'(a,i5,a,i5,a,e12.4)') 'p[',me, &
            '] vector error i: ',idx,' v: ',rv
          chk = 1
        endif
      endif
    endif
  end do
  call ga_igop(2,chk,1,"+")
!  call bv_map%destroy()
  if (me.eq.0) then
    write(6,*)
    if (chk.eq.0) then
      write(6,'(a)') 'Vector elements are ok'
    else
      write(6,'(a)') 'Error found in vector elements'
    endif
  endif
  if (me.eq.0) then
    write(6,*)
    write(6,'(a)') 'Testing map_to_bus'
  endif
!
! Multiply values in vector by factor of 2
! 
  call vvector%local_index_range(lo,hi)
  do i = lo, hi-1
    call vvector%get_element(i,v)
    v = 2.0d00*v
    call vvector%set_element(i,v)
  end do
!
! Push values back onto buses
!
!  call bv_map%destroy()
!  call bv_map%create(grid)
  call bv_map%map_to_bus(vvector)
  chk = 0
  do i = 1, nbus
    if (grid%get_active_bus(i)) then
      bus_ptr => bus_cast(grid%get_bus(i))
      if (bus_ptr%bus_vector_size(isize)) then
        call bus_ptr%bus_get_mat_vec_index(idx)
        rv = bus_ptr%bus_get_value()
        if (rv.ne.dble(2*idx)) then
          write(6,'(a,i4,a,i4,a,f12.4,a,f12.4)') 'p[',me,'] Bus error i: ', &
           idx,' v: ',rv,' expected: ',dble(2*idx)
        endif
      endif
    endif
  end do
  if (me.eq.0) then
    write(6,*)
    if (chk.eq.0) then
      write(6,'(a)') 'Bus values are ok'
    else
      write(6,'(a)') 'Error found in bus values'
    endif
  endif
  call ga_igop(3,chk,1,"+")
  call map_factory%destroy()
  call grid%destroy()
  CALL comm%finalize()
  CALL gridpack_finalize_math()
  CALL gridpack_finalize_parallel()
END PROGRAM mapper_test
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
