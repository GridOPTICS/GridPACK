!
!  Fortran bus network component
!
module gridpack_bus_component
  use, intrinsic :: iso_c_binding
  use gridpack_data_collection
  implicit none
!
!  Define bus component type
!
  private
  type, public :: bus_component
    type (C_PTR) :: p_bus
    class(data_collection), pointer :: p_data
    integer :: p_idx
    integer :: p_mode
    contains
!
!  Matrix-vector interface calls
!
    procedure::matrix_diag_size
    procedure::matrix_diag_values
    procedure::matrix_forward_size
    procedure::matrix_reverse_size
    procedure::matrix_forward_values
    procedure::matrix_reverse_values
    procedure::vector_size
    procedure::vector_values
    procedure::set_values
!
!  Might not need wrappers for these on the Fortran side
!
!    procedure::set_mat_vec_index
!    procedure::get_mat_vec_index
!    procedure::set_mat_vec_indices
!    procedure::get_mat_vec_indices
!
!  Generalized matrix-vector interface calls
!
    procedure::matrix_num_rows
    procedure::matrix_num_cols
    procedure::matrix_set_row_index
    procedure::matrix_set_col_index
    procedure::matrix_get_row_index
    procedure::matrix_get_col_index
    procedure::matrix_num_values
    procedure::matrix_get_values
    procedure::vector_num_elements
    procedure::vector_set_element_index
    procedure::vector_get_element_indices
    procedure::vector_get_element_values
    procedure::vector_set_element_values
!
!  Base component calls
!
    procedure::load
!    procedure::get_xc_buf_size
!    procedure::set_xc_buf
    procedure::set_mode
    procedure::serial_write
!
!  Base bus calls
!
#if 0
    procedure::add_branch
    procedure::add_bus
    procedure::get_neighbor_branches
    procedure::get_neighbor_buses
    procedure::clear_branches
    procedure::clear_buses
    procedure::set_reference_bus
    procedure::get_reference_bus
    procedure::set_original_index
    procedure::get_original_index
    procedure::set_global_index
    procedure::get_global_index
#endif
  end type
!
  type, public :: bus_wrapper
    class(bus_component), pointer :: bus
  end type
!  type(bus_component), target, allocatable :: p_buses(:)
!  type(data_collection), target, allocatable :: p_data(:)
!
!  Interface declaration to C calls
!
  interface
#if 0
!
! Add a branch to the list of branches that a bus is connected to
! @param bus GridPACK bus object
! @param idx index of branch that branch is connected to
!
    subroutine bus_add_branch(bus, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
      integer(C_INT), value, intent(in) :: idx
    end subroutine bus_add_branch
!
! Add a bus to the list of branches that a bus is connected to via a branch
! @param bus GridPACK bus object
! @param idx index of bus that branch is connected to via a branch
!
    subroutine bus_add_bus(bus, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
      integer(C_INT), value, intent(in) :: idx
    end subroutine bus_add_bus
!
! Get the number of neighboring branches/buses
! @param bus GridPACK bus object
! @return number of attached branches/buses
!
    integer(C_INT) function bus_get_num_neighbors(bus) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end function bus_get_num_neighbors
!
! Get indices of branches that are attached to bus
! @param bus GridPACK bus object
! @param nghbrs indices of branches that are attached to bus
!
    subroutine bus_get_neighbor_branches(bus, nghbrs) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
      integer(C_INT), intent(out) :: nghbrs(*)
    end subroutine bus_get_neighbor_branches
!
! Get indices of buses that are attached to bus via a branch
! @param bus GridPACK bus object
! @param nghbrs indices of buses that are attached to bus via a branch
!
    subroutine bus_get_neighbor_buses(bus, nghbrs) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
      integer(C_INT), intent(out) :: nghbrs(*)
    end subroutine bus_get_neighbor_buses
!
! Clear all branches attached to bus
! @param bus GridPACK bus object
!
    subroutine bus_clear_branches(bus) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end subroutine bus_clear_branches
!
! Clear all buses attached to bus
! @param bus GridPACK bus object
!
    subroutine bus_clear_buses(bus) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end subroutine bus_clear_buses
!
! Set reference bus status
! @param bus GridPACK bus object
! @param status true if bus is reference bus
!
    subroutine bus_set_reference_bus(bus, status) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
      type(C_BOOL), value, intent(in) :: status
    end subroutine bus_set_reference_bus
!
! Get reference bus status
! @param bus GridPACK bus object
! @return true if bus is reference bus
!
    logical(C_BOOL) function bus_get_reference_bus(bus) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end function bus_get_reference_bus
!
! Set original index (from input file)
! @param bus GridPACK bus object
! @param idx original index from network
!
    subroutine bus_set_original_index(bus, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end subroutine bus_set_original_index
!
! Get original index (from input file)
! @param bus GridPACK bus object
! @return original index from network
!
    integer(C_INT) function bus_get_original_index(bus) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end function bus_get_original_index
!
! Set global index
! @param bus GridPACK bus object
! @param idx global index from network
!
    subroutine bus_set_global_index(bus, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end subroutine bus_set_global_index
!
! Get global index
! @param bus GridPACK bus object
! @return global index from network
!
    integer(C_INT) function bus_get_global_index(bus) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end function bus_get_global_index
#endif
  end interface
  contains
!
! Return size of matrix block on the diagonal contributed by component
! @param bus GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function matrix_diag_size(bus, isize, jsize)
    implicit none
    class(bus_component), intent(in) :: bus
    integer isize, jsize
    isize = 1
    jsize = 1
    matrix_diag_size = .true.
    return
  end function matrix_diag_size
!
! Return the values for a diagonal matrix block. The values are returned in
! row-major order.
! @param bus GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function matrix_diag_values(bus, values)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(out) :: values(*)
    values(1) = cmplx(1.0d00,0.0d00)
    matrix_diag_values = .true.
    return
  end function matrix_diag_values
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the forward direction
! @param bus GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function matrix_forward_size(bus, isize, jsize)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, intent(out) :: isize, jsize
    isize = 1
    jsize = 1
    matrix_forward_size = .true.
    return
  end function matrix_forward_size
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the reverse direction
! @param bus GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function matrix_reverse_size(bus, isize, jsize)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_component), intent(in) :: bus
    integer, intent(out) :: isize, jsize 
    isize = 1
    jsize = 1
    matrix_reverse_size = .true.
    return
  end function matrix_reverse_size
!
! Return the values of off-diagonal matrix block. The values are for the forward
! direction and are returned in row-major order.
! @param bus GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function matrix_forward_values(bus, values)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(out) :: values(*)
    values(1) = cmplx(1.0d00,0.0d00)
    matrix_forward_values = .true.
    return
  end function matrix_forward_values
!
! Return the values of off-diagonal matrix block. The values are for the reverse
! direction and are returned in row-major order.
! @param bus GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function matrix_reverse_values(bus, values)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(out) :: values(*)
    values(1) = cmplx(1.0d00,0.0d00)
    matrix_reverse_values = .true.
    return
  end function matrix_reverse_values
!
! Return the size of vector block contributed by component
! @param bus GridPACK bus object
! @param isize number of vector elements
! @return false if network component does not contribute vector element
!
  logical function vector_size(bus, isize)
    implicit none
    class(bus_component), intent(in) :: bus
    integer(C_INT), intent(out) :: isize
    isize = 1
    vector_size = .true.
  end function vector_size
!
! Return the values of the vector block
! @param bus GridPACK bus object
! @param values array that contains vector values
! @return false if network component does not contribute vector element
! 
  logical function vector_values(bus, values)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(out) :: values(*)
    values(1) = cmplx(1.0d00,0.0d00)
    vector_values = .true.
  end function vector_values
!
! Set values in network component based on values in a vector or matrix
! @param bus GridPACK bus object
! @param values array that contains vector values
! 
  subroutine set_values(bus, values)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(in) :: values(*)
  end subroutine set_values
!
! Return the number of rows in matrix from component
! @param bus GridPACK bus object
! @return number of rows from component
!
  integer function matrix_num_rows(bus)
    implicit none
    class(bus_component), intent(in) :: bus
  end function matrix_num_rows
!
! Return the number of columns in matrix from component
! @param bus GridPACK bus object
! @return number of columns from component
!
  integer function matrix_num_cols(bus)
    implicit none
    class(bus_component), intent(in) :: bus
  end function matrix_num_cols
!
! Set row indices corresponding to the rows contributed by this component
! @param bus GridPACK bus object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @param idx matrix index of row irow
!
  subroutine matrix_set_row_index(bus, irow, idx)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, value, intent(in) :: irow, idx
  end subroutine matrix_set_row_index
!
! Set column indices corresponding to the columns contributed by this component
! @param bus GridPACK bus object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @param idx matrix index of column icol
!
  subroutine matrix_set_col_index(bus, icol, idx)
    implicit none
    class(bus_component), intent(in) :: bus
    integer(C_INT), value, intent(in) :: icol, idx
  end subroutine matrix_set_col_index
!
! Get the row index for the rows contributed by this component
! @param bus GridPACK bus object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @return matrix index of row irow
!
  integer function matrix_get_row_index(bus, irow)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, value, intent(in) :: irow
  end function matrix_get_row_index
!
! Get the col index for the columns contributed by this component
! @param bus GridPACK bus object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @return matrix index of column icol
!
  integer function matrix_get_col_index(bus, icol)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, value, intent(in) :: icol
  end function matrix_get_col_index
!
! Return the number of matrix values contributed by this component
! @param bus GridPACK bus object
! @return number of matrix values
!
  integer function matrix_num_values(bus)
    implicit none
    class(bus_component), intent(in) :: bus
  end function matrix_num_values
!
! Get a list of matrix values contributed by this components along with their
! matrix indices
! @param bus GridPACK bus object
! @param values list of matrix element values
! @param rows row indices for the matrix elements
! @param cols column indices for the matrix elements
!
  subroutine matrix_get_values(bus, values, rows, cols)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(out) :: values(*)
    integer, intent(out) :: rows(*), cols(*)
  end subroutine matrix_get_values
!
! Return the number of elements in vector coming from component
! @param bus GridPACK bus object
! @return number of elements contributed from component
!
  integer function vector_num_elements(bus)
    implicit none
    class(bus_component), intent(in) :: bus
  end function vector_num_elements
!
! Set indices corresponding to the elements contributed by this component
! @param bus GridPACK bus object
! @param ielem index of element contributed by this component (e.g. if component
! contributes 3 elements then ielem is between 0 and 2)
! @param idx vector index of element ielem
!
  subroutine vector_set_element_index(bus, ielem, idx)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, value, intent(in) :: ielem, idx
  end subroutine vector_set_element_index
!
! Get a list of vector indices from the component
! @param bus GridPACK bus object
! @param idx list of indices the component maps onto
!
  subroutine vector_get_element_indices(bus, idx)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, intent(out) :: idx(*)
  end subroutine vector_get_element_indices
!
! Get a list of vector values contributed by this component and their indices
! @param bus GridPACK bus object
! @param values list of vector element values
! @param idx indices of the vector elements
!
  subroutine vector_get_element_values(bus, values, idx)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(out) :: values(*)
    integer, intent(out) :: idx(*)
  end subroutine vector_get_element_values
!
! Transfer vector values to component
! @param bus GridPACK bus object
! @param values list of vector element values
!
  subroutine vector_set_element_values(bus, values)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(out) :: values(*)
  end subroutine vector_set_element_values
!
! Load data from DataCollection object into corresponding component.
! @param bus GridPACK bus object
! @param data DataCollection object associated with component
!
  subroutine load(bus, data)
    implicit none
    class(bus_component), intent(in) :: bus
    integer data
  end subroutine load
!
! Return the size of the buffer needed for the data exchanges. Note that this
! must be the same size for all buses even if all buses do not require the same
! parameters. Thus, the buffer must be big enough to exchange all variables a
! bus might need, even if individual buses do use all the variables
! @param bus GridPACK bus object
! @return size of buffer in bytes
!
!  integer function get_xc_buf_size(bus)
!    implicit none
!    class(bus_component), intent(in) :: bus
!  end function get_xc_buf_size
!
! Assign the location of the data exchange buffer. These buffers are allocated
! and deallocated by the network
! @param bus GridPACK bus object
! @return size of buffer in bytes
!
!  subroutine set_xc_buf(bus, buf)
!    implicit none
!    class(bus_component), intent(in) :: bus
!    type(C_PTR), value, intent(in) :: buf
!  end subroutine set_xc_buf
!
! Set an internal variable that can be used to control the behavior of the
! component. This function doesn't need to be implemented, but if needed it can
! be used to change the behavior of the network in different phases of the
! calculation. For example, if a different matrix needs to be generated at
! different times, the mode can be change to get different values from the
! matrix-vector interface functions.
! @param bus GridPACK bus object
! @param mode integer indicating which mode should be used
!
  subroutine set_mode(bus, mode)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, value, intent(in) :: mode
  end subroutine set_mode
!
! Copy a string for output into buffer. The behavior of this method can be
! altered by inputting different values for the signal string
! @param bus GridPACK bus object
! @param string buffer containing string to be written to output
! @param bufsize size of string buffer in bytes
! @param signal string to control behavior of routine (e.g. what
! properties to write
! @return true if component is writing a contribution, false otherwise
!
  logical function serial_write(bus, string, bufsize, signal)
    implicit none
    class(bus_component), intent(in) :: bus
    character, intent(out) :: string(*)
    integer, value, intent(in) :: bufsize
    character, intent(in) :: signal(*)
  end function serial_write
!
! Allocate Fortran bus objects in internal array and assign internal index to
! that object
! @param bus_num number of buses on this processor
!
#if 0
  subroutine bus_allocate(bus_num)
    implicit none
    integer bus_num
    integer i
    allocate(p_buses(bus_num))
    allocate(p_data(bus_num))
    do i = 1, bus_num
      p_buses(i)%p_idx = i
      p_data(i)%p_idx = i
      p_buses(i)%p_data => p_data(i)
    end do
    return
  end subroutine bus_allocate 
#endif
  type(C_PTR) function bus_allocate() bind(c)
    implicit none
    type(bus_wrapper), pointer :: wbus
    class(bus_component), pointer :: bus
    allocate(bus)
    allocate(wbus)
    wbus%bus => bus
    bus_allocate = C_LOC(wbus)
    return
  end function bus_allocate
!
! Deallocate Fortran bus objects
!
#if 0
  subroutine bus_deallocate()
    implicit none
    if (allocated(p_buses)) deallocate(p_buses)
    return
  end subroutine bus_deallocate 
#endif
  subroutine bus_deallocate(ptr) bind(c)
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: wbus
    call C_F_POINTER(ptr,wbus)
    if (associated(wbus)) then
      deallocate(wbus%bus)
      deallocate(wbus)
    endif
    return
  end subroutine bus_deallocate
!
! Return size of matrix block on the diagonal contributed by bus
! @param ptr C pointer to GridPACK bus object
! @param isize,jsize number of rows and columns in block
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function bus_matrix_diag_size(ptr, isize, jsize) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer(C_INT), intent(out) :: isize, jsize
    logical f_ret
    integer f_isize, f_jsize
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%matrix_diag_size(f_isize, f_jsize)
    isize = f_isize
    jsize = f_isize
    bus_matrix_diag_size = f_ret
    return
  end function bus_matrix_diag_size
!
! Return the values for a diagonal matrix block. The values are returned in
! row-major order.
! @param ptr C pointer to GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function bus_matrix_diag_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    logical f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%matrix_diag_values(values)
    bus_matrix_diag_values = f_ret
    return
  end function bus_matrix_diag_values
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the forward direction
! @param ptr C pointer to GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function bus_matrix_forward_size(ptr, isize, jsize) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer(C_INT), intent(out) :: isize, jsize
    logical f_ret
    integer f_isize, f_jsize
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%matrix_forward_size(f_isize, f_jsize)
    isize = f_isize
    jsize = f_jsize
    bus_matrix_forward_size = f_ret
    return
  end function bus_matrix_forward_size
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the reverse direction
! @param ptr C pointer to GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function bus_matrix_reverse_size(ptr, isize, jsize) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer(C_INT), intent(out) :: isize, jsize 
    logical f_ret
    integer f_isize, f_jsize
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%matrix_reverse_size(f_isize, f_jsize)
    isize = f_isize
    jsize = f_jsize
    bus_matrix_reverse_size = f_ret
    return
  end function bus_matrix_reverse_size
!
! Return the values of off-diagonal matrix block. The values are for the forward
! direction and are returned in row-major order.
! @param ptr C pointer to GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function bus_matrix_forward_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    logical f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%matrix_forward_values(values)
    bus_matrix_forward_values = f_ret
    return
  end function bus_matrix_forward_values
!
! Return the values of off-diagonal matrix block. The values are for the reverse
! direction and are returned in row-major order.
! @param ptr C pointer to GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function bus_matrix_reverse_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    logical f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%matrix_reverse_values(values)
    bus_matrix_reverse_values = f_ret
    return
  end function bus_matrix_reverse_values
!
! Return the size of vector block contributed by component
! @param ptr C pointer to GridPACK bus object
! @param isize number of vector elements
! @return false if network component does not contribute vector element
!
  logical(C_BOOL) function bus_vector_size(ptr, isize) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer(C_INT), intent(out) :: isize
    integer f_isize
    logical f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%vector_size(f_isize)
    bus_vector_size = f_ret
  end function bus_vector_size
!
! Return the values of the vector block
! @param ptr C pointer to GridPACK bus object
! @param values array that contains vector values
! @return false if network component does not contribute vector element
! 
  logical(C_BOOL) function bus_vector_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    logical f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%vector_values(values)
    bus_vector_values = f_ret
  end function bus_vector_values
!
! Set values in network component based on values in a vector or matrix
! @param ptr C pointer to GridPACK bus object
! @param values array that contains vector values
! 
  subroutine bus_set_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(in) :: values(*)
    call C_F_POINTER(ptr,bus)
    call bus%bus%set_values(values)
  end subroutine bus_set_values
!
! Return the number of rows in matrix from component
! @param ptr C pointer to GridPACK bus object
! @return number of rows from component
!
  integer(C_INT) function bus_matrix_num_rows(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%matrix_num_rows()
    bus_matrix_num_rows = f_ret
  end function bus_matrix_num_rows
!
! Return the number of columns in matrix from component
! @param ptr C pointer to GridPACK bus object
! @return number of columns from component
!
  integer(C_INT) function bus_matrix_num_cols(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%matrix_num_cols()
    bus_matrix_num_cols = f_ret
  end function bus_matrix_num_cols
!
! Set row indices corresponding to the rows contributed by this component
! @param ptr C pointer to GridPACK bus object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @param idx matrix index of row irow
!
  subroutine bus_matrix_set_row_index(ptr, irow, r_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer(C_INT), value, intent(in) :: irow, r_idx
    integer f_row, f_idx
    call C_F_POINTER(ptr,bus)
    f_row = irow
    f_idx = r_idx
    call bus%bus%matrix_set_row_index(f_row, f_idx)
    return
  end subroutine bus_matrix_set_row_index
!
! Set column indices corresponding to the columns contributed by this component
! @param ptr C pointer to GridPACK bus object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @param idx matrix index of column icol
!
  subroutine bus_matrix_set_col_index(ptr, icol, c_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer(C_INT), value, intent(in) :: icol, c_idx
    integer f_col, f_idx
    call C_F_POINTER(ptr,bus)
    f_col = icol
    f_idx = c_idx
    call bus%bus%matrix_set_col_index(f_col, f_idx)
    return
  end subroutine bus_matrix_set_col_index
!
! Get the row index for the rows contributed by this component
! @param ptr C pointer to GridPACK bus object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @return matrix index of row irow
!
  integer(C_INT) function bus_matrix_get_row_index(ptr, irow) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer(C_INT), value, intent(in) :: irow
    integer f_row, f_ret
    call C_F_POINTER(ptr,bus)
    f_row = irow
    f_ret = bus%bus%matrix_get_row_index(f_row)
    bus_matrix_get_row_index = f_ret
    return
  end function bus_matrix_get_row_index
!
! Get the col index for the columns contributed by this component
! @param ptr C pointer to GridPACK bus object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @return matrix index of column icol
!
  integer(C_INT) function bus_matrix_get_col_index(ptr, icol) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer(C_INT), value, intent(in) :: icol
    integer f_col, f_ret
    call C_F_POINTER(ptr,bus)
    f_col = icol
    f_ret = bus%bus%matrix_get_col_index(f_col)
    bus_matrix_get_col_index = f_ret
  end function bus_matrix_get_col_index
!
! Return the number of matrix values contributed by this component
! @param ptr C pointer to GridPACK bus object
! @return number of matrix values
!
  integer(C_INT) function bus_matrix_num_values(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%matrix_num_values()
    bus_matrix_num_values = f_ret
  end function bus_matrix_num_values
!
! Get a list of matrix values contributed by this components along with their
! matrix indices
! @param ptr C pointer to GridPACK bus object
! @param values list of matrix element values
! @param rows row indices for the matrix elements
! @param cols column indices for the matrix elements
!
  subroutine bus_matrix_get_values(ptr, values, rows, cols) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    integer(C_INT), intent(out) :: rows(*), cols(*)
    integer nval, i
    integer, allocatable :: f_rows(:), f_cols(:)
    call C_F_POINTER(ptr,bus)
    nval = bus%bus%matrix_num_values()
    allocate(f_rows(nval))
    allocate(f_cols(nval))
    call bus%bus%matrix_get_values(values,f_rows,f_cols)
    do i = 1, nval
      rows(i) = f_rows(i)
      cols(i) = f_cols(i)
    end do
    deallocate(f_rows)
    deallocate(f_cols)
    return
  end subroutine bus_matrix_get_values
!
! Return the number of elements in vector coming from component
! @param ptr C pointer to GridPACK bus object
! @return number of elements contributed from component
!
  integer(C_INT) function bus_vector_num_elements(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%vector_num_elements()
    bus_vector_num_elements = f_ret
    return
  end function bus_vector_num_elements
!
! Set indices corresponding to the elements contributed by this component
! @param ptr C pointer to GridPACK bus object
! @param ielem index of element contributed by this component (e.g. if component
! contributes 3 elements then ielem is between 0 and 2)
! @param idx vector index of element ielem
!
  subroutine bus_vector_set_element_index(ptr, ielem, e_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer(C_INT), value, intent(in) :: ielem, e_idx
    integer f_elem, f_idx
    call C_F_POINTER(ptr,bus)
    f_elem = ielem
    f_idx = e_idx
    call bus%bus%vector_set_element_index(f_elem, f_idx)
    return
  end subroutine bus_vector_set_element_index
!
! Get a list of vector indices from the component
! @param ptr C pointer to GridPACK bus object
! @param idx list of indices the component maps onto
!
  subroutine bus_vector_get_element_indices(ptr, e_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer(C_INT), intent(out) :: e_idx(*)
    integer, allocatable :: f_idx(:)
    integer nval, i
    call C_F_POINTER(ptr,bus)
    nval = bus%bus%vector_num_elements()
    allocate(f_idx(nval))
    call bus%bus%vector_get_element_indices(f_idx)
    do i = 1, nval
      e_idx(i) = f_idx(i)
    end do
    deallocate(f_idx)
    return
  end subroutine bus_vector_get_element_indices
!
! Get a list of vector values contributed by this component and their indices
! @param ptr C pointer to GridPACK bus object
! @param values list of vector element values
! @param idx indices of the vector elements
!
  subroutine bus_vector_get_element_values(ptr, values, e_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    integer(C_INT), intent(out) :: e_idx(*)
    integer, allocatable :: f_idx(:)
    integer nval, i
    call C_F_POINTER(ptr,bus)
    nval = bus%bus%vector_num_elements()
    allocate(f_idx(nval))
    call bus%bus%vector_get_element_values(values, f_idx)
    do i = 1, nval
      e_idx(i) = f_idx(i)
    end do
    deallocate(f_idx)
  end subroutine bus_vector_get_element_values
!
! Transfer vector values to component
! @param ptr C pointer to GridPACK bus object
! @param values list of vector element values
!
  subroutine bus_vector_set_element_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    call C_F_POINTER(ptr,bus)
    call bus%bus%vector_set_element_values(values)
    return
  end subroutine bus_vector_set_element_values
!
! Load data from DataCollection object into corresponding component.
! @param ptr C pointer to GridPACK bus object
! @param data DataCollection object associated with component
!
  subroutine bus_load(ptr, data) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer(C_INT) data
    integer f_data
    call C_F_POINTER(ptr,bus)
    f_data = data
    call bus%bus%load(f_data)
    return
  end subroutine bus_load
!
! Return the size of the buffer needed for the data exchanges. Note that this
! must be the same size for all buses even if all buses do not require the same
! parameters. Thus, the buffer must be big enough to exchange all variables a
! bus might need, even if individual buses do use all the variables
! @param ptr C pointer to GridPACK bus object
! @return size of buffer in bytes
!
!  integer(C_INT) function bus_get_xc_buf_size(ptr) bind(c)
!    use, intrinsic :: iso_c_binding
!    implicit none
!    type(C_PTR), value, intent(in) :: ptr
!    class(bus_wrapper), pointer :: bus
!    integer f_ret
!    call C_F_POINTER(ptr,bus)
!    f_ret = bus%bus%get_xc_buf_size()
!    bus_get_xc_buf_size = f_ret
!    return
!  end function bus_get_xc_buf_size
!
! Assign the location of the data exchange buffer. These buffers are allocated
! and deallocated by the network
! @param ptr C pointer to GridPACK bus object
! @return size of buffer in bytes
!
!  subroutine bus_set_xc_buf(ptr, buf) bind(c)
!    use, intrinsic :: iso_c_binding
!    implicit none
!    type(C_PTR), value, intent(in) :: ptr
!    class(bus_wrapper), pointer :: bus
!    call C_F_POINTER(ptr,bus)
!    call bus%bus%set_xc_buf(buf)
!    return
!  end subroutine bus_set_xc_buf
!
! Set an internal variable that can be used to control the behavior of the
! component. This function doesn't need to be implemented, but if needed it can
! be used to change the behavior of the network in different phases of the
! calculation. For example, if a different matrix needs to be generated at
! different times, the mode can be change to get different values from the
! matrix-vector interface functions.
! @param ptr C pointer to GridPACK bus object
! @param mode integer indicating which mode should be used
!
  subroutine bus_set_mode(ptr, mode) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    integer(C_INT), value, intent(in) :: mode
    integer f_mode
    call C_F_POINTER(ptr,bus)
    f_mode = mode
    call bus%bus%set_mode(f_mode)
    return
  end subroutine bus_set_mode
!
! Copy a string for output into buffer. The behavior of this method can be
! altered by inputting different values for the signal string
! @param ptr C pointer to GridPACK bus object
! @param string buffer containing string to be written to output
! @param bufsize size of string buffer in bytes
! @param signal string to control behavior of routine (e.g. what
! properties to write
! @return true if component is writing a contribution, false otherwise
!
  logical(C_BOOL) function bus_serial_write(ptr, string, bufsize, &
      signal) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    class(bus_wrapper), pointer :: bus
    character(C_CHAR), intent(out) :: string(*)
    integer(C_INT), value, intent(in) :: bufsize
    character(C_CHAR), intent(in) :: signal(*)
    logical f_ret
    integer f_bufsize
    call C_F_POINTER(ptr,bus)
    f_bufsize = bufsize
    f_ret = bus%bus%serial_write(string,f_bufsize,signal)
    bus_serial_write = f_ret
    return
  end function bus_serial_write
end module gridpack_bus_component
