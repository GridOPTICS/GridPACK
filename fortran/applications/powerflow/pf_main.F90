! ----------------------------------------------------------------
! file: pf_main.f90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created August 27, 2014 by Bruce Palmer
! Last Change: 2014-09-27 14:07:10 d3g293
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! PROGRAM pf_main
! ----------------------------------------------------------------
PROGRAM powerflow
  USE, intrinsic :: iso_c_binding
  USE gridpack_communicator
  USE gridpack_configuration
  USE gridpack_network
  USE gridpack_parser
  USE gridpack_full_matrix_map
  USE gridpack_bus_vector_map
  USE gridpack_serial_io
  USE gridpack_parallel
  USE gridpack_math
  USE gridpack_vector
  USE gridpack_matrix
  USE gridpack_linear_solver
  USE application_components
  USE application_factory
  USE gridpack_component

  IMPLICIT NONE
  TYPE (communicator) :: comm
  INTEGER :: rank, size, p
  INTEGER(c_int) :: ma_stack = 200000, ma_heap = 200000
  TYPE(cursor) :: config, crs
  TYPE(network) :: grid
  type(app_factory) :: factory
  type(pti23_parser) :: parser
  type(matrix) :: y, j
  type(vector) :: pq, x
  type(full_matrix_map) :: j_map, y_map
  type(bus_vector_map) :: v_map
  type(bus_serial_io) bus_io
  type(branch_serial_io) branch_io
  type(linear_solver) solver
  logical ok
  character(512) filename, inputfile
  character(512) iobuf
  character(512) dbug
  double precision tolerance
  integer max_iteration, iter
  double complex tol
!
! Initialize gridpack and create world communicator
!
  CALL gridpack_initialize_parallel(ma_stack, ma_heap)
  CALL comm%initialize()
  CALL gridpack_initialize_math()
!
!  Initialize network
!
  call grid%create(comm)  
!
  call config%initialize()
  call get_command_argument(1,inputfile)
  ok = .false.
  if (len_trim(inputfile).gt.0) then
    ok = config%open(comm,trim(inputfile))
  endif
  if (.not.ok) ok = config%open(comm,"input.xml")
  if (.not.ok) stop
  crs = config%get_cursor("Configuration.Powerflow")
  if (.not.crs%get_string("networkConfiguration",filename)) then
    write(6,'(a)') 'No network configuration file specified'
    stop
  endif
!
  ok = crs%get_double("tolerance",tolerance)
  ok = crs%get_int("maxIteration",max_iteration)
!
! Set up parser
!
  call parser%create(grid)
  call parser%parse(trim(filename))
  call parser%destroy()
!
  call bus_io%bus_create(512,grid)
  call bus_io%bus_header("")
  write(iobuf,'(a,i6)') 'Maximum number of iterations: ',max_iteration
  call bus_io%bus_header(iobuf)
  call bus_io%bus_header("")
  write(iobuf,'(a,e12.4)') 'Convergence tolerance: ',tolerance
  call bus_io%bus_header(iobuf)
!
! Partition network
!
  call grid%partition()
!
! Create factory
!
  call factory%create(grid)
  call factory%load()
!
! Set network components
!
  call factory%set_components()
!
! Configure exchange buffers
!
  call factory%set_exchange()
!
! Initialize bus data exchange
!
  call grid%init_bus_update() 
!
! Set up y-matrix components
!
  call factory%set_mode(YBUS)
  call factory%set_y_bus()
!
!  call y_map%create(grid)
!  y = y_map%map_to_matrix()
!  call bus_io%bus_header("")
!  call bus_io%bus_header("Y-matrix")
!  call y%print()
!
  call factory%set_s_bus()
  call bus_io%bus_header("")
  call bus_io%bus_header("Iteration 0")
!
  call factory%set_mode(RHS)
  call v_map%create(grid)
  pq = v_map%map_to_vector()
!  call bus_io%bus_header("")
!  call bus_io%bus_header("PQ vector")
!  call pq%print()
!
  call factory%set_mode(JACOBIAN)
  call j_map%create(grid)
  j = j_map%map_to_matrix()
!  call bus_io%bus_header("")
!  call bus_io%bus_header("Jacobian matrix")
!  call j%print()
!
  x = pq%clone()
!
  call solver%initialize(j,crs)
!
  tol = dcmplx(2.0d00*tolerance,0.0d00)
  iter = 0
!
! First iteration
!
  call x%zero()
  call bus_io%bus_header("")
  call bus_io%bus_header("Calling solver")
  call solver%solve(pq,x)
!  call bus_io%bus_header("")
!  call bus_io%bus_header("Solution vector")
!  call x%print
  tol = pq%norm_infinity()
!
  do while (real(tol).gt.tolerance.and.iter.lt.max_iteration)
    call factory%set_mode(RHS)
    call v_map%map_to_bus(x)
    call grid%update_buses()
    call v_map%remap_to_vector(pq)
!  call bus_io%bus_header("")
!  call bus_io%bus_header("new PQ vector")
!  call pq%print
    call factory%set_mode(JACOBIAN)
    call j_map%remap_to_matrix(j)
    call x%zero()
    call solver%solve(pq,x)
    tol = pq%norm_infinity()
    call bus_io%bus_header("")
    write(iobuf,'(a,i4,a,e12.6)') 'Iteration ',iter+1,' Tol: ',real(tol)
    call bus_io%bus_header(iobuf)
    iter = iter + 1
  end do
!
! Push final results back onto buses
!
  call factory%set_mode(RHS)
  call v_map%map_to_bus(x)
  call grid%update_buses()
!
  call branch_io%branch_create(512,grid)
  call branch_io%branch_header("")
  call branch_io%branch_header("   Branch Power Flow")
  call branch_io%branch_header("")
  call branch_io%branch_header("        Bus 1       Bus 2   CKT         P" // &
                        "                    Q")
  call branch_io%branch_write("")
!
  call bus_io%bus_header("")
  call bus_io%bus_header("   Bus Voltages and Phase Angles")
  call branch_io%branch_header("")
  call bus_io%bus_header(" Number      Phase Angle      Voltage Magnitude")
  call bus_io%bus_write("")
!
!  Close out calculation
!
  call bus_io%bus_destroy()
  call branch_io%branch_destroy()
  call factory%destroy()
  call y_map%destroy()
  call j_map%destroy()
  call v_map%destroy()
  call grid%destroy()
  call crs%finalize()
  call config%finalize()
  CALL comm%finalize()
  CALL gridpack_finalize_math()
  CALL gridpack_finalize_parallel()
END PROGRAM powerflow
