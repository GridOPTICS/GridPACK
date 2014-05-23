! Wrappers for Fortran bus component functions
!
! Return size of matrix block on the diagonal contributed by a bus component
!
logical(C_BOOL) function p_bus_matrix_diag_size(c_network, c_idx, &
       c_isize, c_jsize) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  integer(C_INT), intent(out) :: c_isize, c_jsize
  integer network, idx, isize, jsize
  logical bus_matrix_diag_size
  network = c_network
  idx = c_idx
  p_bus_matrix_diag_size = bus_matrix_diag_size(network, idx, isize, jsize)
  c_isize = isize
  c_jsize = jsize
  return
end function p_bus_matrix_diag_size
!
! Return value of matrix block on the diagonal contributed by a bus component.
! Values are returned in row-major order
!
logical(C_BOOL) function p_bus_matrix_diag_values(c_network, c_idx, &
       c_values) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  complex(C_DOUBLE_COMPLEX), intent(out) :: c_values(*)
  integer network, idx
  logical bus_matrix_diag_values
  network = c_network
  idx = c_idx
  p_bus_matrix_diag_values = bus_matrix_diag_values(network, idx, c_values)
  return
end function p_bus_matrix_diag_values
!
! Return size of off-diagonal matrix block contributed by bus component. The
! values are for the forward direction
!
logical(C_BOOL) function p_bus_matrix_forward_size(c_network, c_idx, &
       c_isize, c_jsize) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  integer(C_INT), intent(out) :: c_isize, c_jsize
  integer network, idx, isize, jsize
  logical bus_matrix_forward_size
  network = c_network
  idx = c_idx
  p_bus_matrix_forward_size = bus_matrix_forward_size(network, idx, isize, jsize)
  c_isize = isize
  c_jsize = jsize
  return
end function p_bus_matrix_forward_size
!
! Return size of off-diagonal matrix block contributed by bus component. The
! values are for the reverse direction
!
logical(C_BOOL) function p_bus_matrix_reverse_size(c_network, c_idx, &
       c_isize, c_jsize) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  integer(C_INT), intent(out) :: c_isize, c_jsize
  integer network, idx, isize, jsize
  logical bus_matrix_reverse_size
  network = c_network
  idx = c_idx
  p_bus_matrix_reverse_size = bus_matrix_reverse_size(network, idx, isize, jsize)
  c_isize = isize
  c_jsize = jsize
  return
end function p_bus_matrix_reverse_size
!
! Return the values of an off-diagonal matrix block contributed by a bus
! component. The values are for the forward direction and are returned in
! row-major order.
!
logical(C_BOOL) function p_bus_matrix_forward_values(c_network, c_idx, &
       c_values) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  complex(C_DOUBLE_COMPLEX), intent(out) :: c_values(*)
  integer network, idx
  logical bus_matrix_forward_values
  network = c_network
  idx = c_idx
  p_bus_matrix_forward_values = bus_matrix_forward_values(network, idx, c_values)
  return
end function p_bus_matrix_forward_values
!
! Return the values of an off-diagonal matrix block contributed by a bus
! component. The values are for the reverse direction and are returned in
! row-major order.
!
logical(C_BOOL) function p_bus_matrix_reverse_values(c_network, c_idx, &
       c_values) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  complex(C_DOUBLE_COMPLEX), intent(out) :: c_values(*)
  integer network, idx
  logical bus_matrix_reverse_values
  network = c_network
  idx = c_idx
  p_bus_matrix_reverse_values = bus_matrix_reverse_values(network, idx, c_values)
  return
end function p_bus_matrix_reverse_values
!
! Return the size of vector block contributed by bus component
!
logical(C_BOOL) function p_bus_vector_size(c_network, c_idx, c_size) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  integer(C_INT), intent(out) :: c_size
  integer network, idx, nsize
  logical bus_vector_size
  network = c_network
  idx = c_idx
  p_bus_vector_size = bus_vector_size(network, idx, nsize)
  c_size = nsize
  return
end function p_bus_vector_size
!
! Set the values in the bus component based on values in a vector or matrix
!
subroutine p_bus_set_values(c_network, c_idx, c_values) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  complex(C_DOUBLE_COMPLEX), intent(in) :: c_values(*)
  integer network, idx
  network = c_network
  idx = c_idx
  call bus_set_values(network,idx,c_values)
  return
end subroutine p_bus_set_values
!
! Return the values of a vector block
!
logical(C_BOOL) function p_bus_vector_values(c_network, c_idx, &
       c_values) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  complex(C_DOUBLE_COMPLEX), intent(out) :: c_values(*)
  integer network, idx
  logical bus_vector_values
  network = c_network
  idx = c_idx
  p_bus_vector_values = bus_vector_values(network, idx, c_values)
  return
end function p_bus_vector_values
!
! Load data from data collection into corresponding bus component.
!
subroutine p_bus_load(c_network, c_idx) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  integer network, idx
  network = c_network
  idx = c_idx
  call bus_load(network, idx)
  return
end subroutine p_bus_load
!
! Return the size of the buffer needed for data exchanges. Note that this must
! be the same size for all buses even if all buses do not need to exchange the
! parameters. Thus, the buffer must be big enough to exchange all variables
! that a bus might need, even if individual buses don't need all the variables
!
integer(C_INT) function p_bus_get_xc_buf_size(c_network, c_idx) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  integer network, idx
  integer bus_get_xc_buf_size
  network = c_network
  idx = c_idx
  p_bus_get_xc_buf_size = bus_get_xc_buf_size(network, idx)
  return 
end function p_bus_get_xc_buf_size
!
! Set an internal variable that can be used to control the behavior of the
! component. This function doesn't need to be implemented, but if needed, it can
! be used to change the behavior of the network in different phases of the
! calculation. For example, if a different matrix needs to be generated at
! different times, the mode of the calculation can be changed to get different
! values from the matrix interface functions.
!
subroutine p_bus_set_mode(c_network, c_idx, c_mode) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx, c_mode
  integer network, idx, mode
  network = c_network
  idx = c_idx
  mode = c_mode
  call bus_set_mode(network, idx, mode)
  return
end subroutine p_bus_set_mode
!
! Copy a string for output into a buffer. The behavior of this method can be
! altered by inputting different values for the signal string.
!
logical(C_BOOL) function p_bus_serial_write(c_network, c_idx, c_string, &
       c_bufsize, c_signal) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx, c_bufsize
  character(kind=C_CHAR), intent(out) :: c_string(*)
  character(kind=C_CHAR), intent(in) :: c_signal(*)
  integer network, idx, bufsize
  logical bus_serial_write
  network = c_network
  idx = c_idx
  bufsize = c_bufsize
  p_bus_serial_write = bus_serial_write(network,idx,c_string,bufsize,c_signal)
end function p_bus_serial_write
!
! Wrappers for Fortran branch component functions
!
! Return size of matrix block on the diagonal contributed by a branch component
!
logical(C_BOOL) function p_branch_matrix_diag_size(c_network, c_idx, &
       c_isize, c_jsize) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  integer(C_INT), intent(out) :: c_isize, c_jsize
  integer network, idx, isize, jsize
  logical branch_matrix_diag_size
  network = c_network
  idx = c_idx
  p_branch_matrix_diag_size = branch_matrix_diag_size(network, idx, isize, jsize)
  c_isize = isize
  c_jsize = jsize
  return
end function p_branch_matrix_diag_size
!
! Return value of matrix block on the diagonal contributed by a branch component.
! Values are returned in row-major order
!
logical(C_BOOL) function p_branch_matrix_diag_values(c_network, c_idx, &
       c_values) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  complex(C_DOUBLE_COMPLEX), intent(out) :: c_values(*)
  integer network, idx, nsize, i
  double complex, allocatable :: values(:)
  logical branch_matrix_diag_values
  network = c_network
  idx = c_idx
  p_branch_matrix_diag_values = branch_matrix_diag_values(network, idx, c_values)
  return
end function p_branch_matrix_diag_values
!
! Return size of off-diagonal matrix block contributed by branch component. The
! values are for the forward direction
!
logical(C_BOOL) function p_branch_matrix_forward_size(c_network, c_idx, &
       c_isize, c_jsize) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  integer(C_INT), intent(out) :: c_isize, c_jsize
  integer network, idx, isize, jsize
  logical branch_matrix_forward_size
  network = c_network
  idx = c_idx
  p_branch_matrix_forward_size = branch_matrix_forward_size(network, idx, isize, jsize)
  c_isize = isize
  c_jsize = jsize
  return
end function p_branch_matrix_forward_size
!
! Return size of off-diagonal matrix block contributed by branch component. The
! values are for the reverse direction
!
logical(C_BOOL) function p_branch_matrix_reverse_size(c_network, c_idx, &
       c_isize, c_jsize) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  integer(C_INT), intent(out) :: c_isize, c_jsize
  integer network, idx, isize, jsize
  logical branch_matrix_reverse_size
  network = c_network
  idx = c_idx
  p_branch_matrix_reverse_size = branch_matrix_reverse_size(network, idx, isize, jsize)
  c_isize = isize
  c_jsize = jsize
  return
end function p_branch_matrix_reverse_size
!
! Return the values of an off-diagonal matrix block contributed by a branch
! component. The values are for the forward direction and are returned in
! row-major order.
!
logical(C_BOOL) function p_branch_matrix_forward_values(c_network, c_idx, &
       c_values) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  complex(C_DOUBLE_COMPLEX), intent(out) :: c_values(*)
  integer network, idx
  logical branch_matrix_forward_values
  network = c_network
  idx = c_idx
  p_branch_matrix_forward_values = branch_matrix_forward_values(network, idx, c_values)
  return
end function p_branch_matrix_forward_values
!
! Return the values of an off-diagonal matrix block contributed by a branch
! component. The values are for the reverse direction and are returned in
! row-major order.
!
logical(C_BOOL) function p_branch_matrix_reverse_values(c_network, c_idx, &
       c_values) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  complex(C_DOUBLE_COMPLEX), intent(out) :: c_values(*)
  integer network, idx
  logical branch_matrix_reverse_values
  network = c_network
  idx = c_idx
  p_branch_matrix_reverse_values = branch_matrix_reverse_values(network, idx, c_values)
  return
end function p_branch_matrix_reverse_values
!
! Return the size of vector block contributed by branch component
!
logical(C_BOOL) function p_branch_vector_size(c_network, c_idx, c_size) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  integer(C_INT), intent(out) :: c_size
  integer network, idx, nsize
  logical branch_vector_size
  network = c_network
  idx = c_idx
  p_branch_vector_size = branch_vector_size(network, idx, nsize)
  c_size = nsize
  return
end function p_branch_vector_size
!
! Set the values in the branch component based on values in a vector or matrix
!
subroutine p_branch_set_values(c_network, c_idx, c_values) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  complex(C_DOUBLE_COMPLEX), intent(in) :: c_values(*)
  integer network, idx
  network = c_network
  idx = c_idx
  call branch_set_values(network,idx,c_values)
  return
end subroutine p_branch_set_values
!
! Return the values of a vector block
!
logical(C_BOOL) function p_branch_vector_values(c_network, c_idx, &
       c_values) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  complex(C_DOUBLE_COMPLEX), intent(out) :: c_values(*)
  integer network, idx
  logical branch_vector_values
  network = c_network
  idx = c_idx
  p_branch_vector_values = branch_vector_values(network, idx, c_values)
  return
end function p_branch_vector_values
!
! Load data from data collection into corresponding branch component.
!
subroutine p_branch_load(c_network, c_idx) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  integer network, idx
  network = c_network
  idx = c_idx
  call branch_load(network, idx)
  return
end subroutine p_branch_load
!
! Return the size of the buffer needed for data exchanges. Note that this must
! be the same size for all branches even if all branches do not need to exchange the
! parameters. Thus, the buffer must be big enough to exchange all variables
! that a branch might need, even if individual branches don't need all the variables
!
integer(C_INT) function p_branch_get_xc_buf_size(c_network, c_idx) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx
  integer network, idx
  integer branch_get_xc_buf_size
  network = c_network
  idx = c_idx
  p_branch_get_xc_buf_size = branch_get_xc_buf_size(network, idx)
  return 
end function p_branch_get_xc_buf_size
!
! Set an internal variable that can be used to control the behavior of the
! component. This function doesn't need to be implemented, but if needed, it can
! be used to change the behavior of the network in different phases of the
! calculation. For example, if a different matrix needs to be generated at
! different times, the mode of the calculation can be changed to get different
! values from the matrix interface functions.
!
subroutine p_branch_set_mode(c_network, c_idx, c_mode) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx, c_mode
  integer network, idx, mode
  network = c_network
  idx = c_idx
  mode = c_mode
  call branch_set_mode(network, idx, mode)
  return
end subroutine p_branch_set_mode
!
! Copy a string for output into a buffer. The behavior of this method can be
! altered by inputting different values for the signal string.
!
logical(C_BOOL) function p_branch_serial_write(c_network, c_idx, c_string, &
       c_bufsize, c_signal) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), intent(in) :: c_network, c_idx, c_bufsize
  character(kind=C_CHAR), intent(out) :: c_string(*)
  character(kind=C_CHAR), intent(in) :: c_signal(*)
  integer network, idx, bufsize
  logical branch_serial_write
  network = c_network
  idx = c_idx
  bufsize = c_bufsize
  p_branch_serial_write = branch_serial_write(network,idx,c_string,bufsize,c_signal)
end function p_branch_serial_write
