! ----------------------------------------------------------------
! file: bus_component_f.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created September 30, 2014 by Bruce Palmer
! ----------------------------------------------------------------
!
!  Fortran bus network component
!
module gridpack_component
  use, intrinsic :: iso_c_binding
  use gridpack_data_collection
  implicit none
!
!  Define bus component type
!
!  private
!
  type, abstract, public :: bus_component
    type (C_PTR) :: p_bus
    type (C_PTR) :: c_this
    contains
!
!  Matrix-vector interface calls
!
    procedure::bus_matrix_diag_size => base_bus_matrix_diag_size
    procedure::bus_matrix_diag_values => base_bus_matrix_diag_values
    procedure::bus_matrix_forward_size => base_bus_matrix_forward_size
    procedure::bus_matrix_reverse_size => base_bus_matrix_reverse_size
    procedure::bus_matrix_forward_values => base_bus_matrix_forward_values
    procedure::bus_matrix_reverse_values => base_bus_matrix_reverse_values
    procedure::bus_vector_size => base_bus_vector_size
    procedure::bus_vector_values => base_bus_vector_values
    procedure::bus_set_values => base_bus_set_values
    procedure, non_overridable::bus_get_mat_vec_index
!
!  Generalized matrix-vector interface calls
!
    procedure::bus_matrix_num_rows => base_bus_matrix_num_rows
    procedure::bus_matrix_num_cols => base_bus_matrix_num_cols
    procedure::bus_matrix_set_row_index => base_bus_matrix_set_row_index
    procedure::bus_matrix_set_col_index => base_bus_matrix_set_col_index
    procedure::bus_matrix_get_row_index => base_bus_matrix_get_row_index
    procedure::bus_matrix_get_col_index => base_bus_matrix_get_col_index
    procedure::bus_matrix_num_values => base_bus_matrix_num_values
    procedure::bus_matrix_get_values => base_bus_matrix_get_values
    procedure::bus_vector_num_elements => base_bus_vector_num_elements
    procedure::bus_vector_set_element_index => base_bus_vector_set_element_index
    procedure::bus_vector_get_element_indices => base_bus_vector_get_element_indices
    procedure::bus_vector_get_element_values => base_bus_vector_get_element_values
    procedure::bus_vector_set_element_values => base_bus_vector_set_element_values
!
!  Base component calls
!
    procedure::bus_load => base_bus_load
    procedure::bus_set_mode => base_bus_set_mode
    procedure::bus_serial_write => base_bus_serial_write
    procedure(i_bus_get_xc_buf_size), deferred::bus_get_xc_buf_size
    procedure(i_bus_get_xc_buf), deferred::bus_get_xc_buf
!
!  Base bus calls
!
    procedure, non_overridable::bus_get_num_neighbors
    procedure, non_overridable::bus_clear_branches
    procedure, non_overridable::bus_clear_buses
    procedure, non_overridable::bus_set_reference_bus
    procedure, non_overridable::bus_get_reference_bus
    procedure, non_overridable::bus_get_original_index
    procedure, non_overridable::bus_get_global_index
    procedure, non_overridable::bus_compare
  end type
!
!  Define branch component type
!
  type, abstract, public :: branch_component
    type (C_PTR) :: p_branch
    type (C_PTR) :: c_this
    contains
!
!  Matrix-vector interface calls
!
    procedure::branch_matrix_diag_size => base_branch_matrix_diag_size
    procedure::branch_matrix_diag_values => base_branch_matrix_diag_values
    procedure::branch_matrix_forward_size => base_branch_matrix_forward_size
    procedure::branch_matrix_reverse_size => base_branch_matrix_reverse_size
    procedure::branch_matrix_forward_values => base_branch_matrix_forward_values
    procedure::branch_matrix_reverse_values => base_branch_matrix_reverse_values
    procedure::branch_vector_size => base_branch_vector_size
    procedure::branch_vector_values => base_branch_vector_values
    procedure::branch_set_values => base_branch_set_values
    procedure, non_overridable::branch_get_mat_vec_indices
!
!  Generalized matrix-vector interface calls
!
    procedure::branch_matrix_num_rows => base_branch_matrix_num_rows
    procedure::branch_matrix_num_cols => base_branch_matrix_num_cols
    procedure::branch_matrix_set_row_index => base_branch_matrix_set_row_index
    procedure::branch_matrix_set_col_index => base_branch_matrix_set_col_index
    procedure::branch_matrix_get_row_index => base_branch_matrix_get_row_index
    procedure::branch_matrix_get_col_index => base_branch_matrix_get_col_index
    procedure::branch_matrix_num_values => base_branch_matrix_num_values
    procedure::branch_matrix_get_values => base_branch_matrix_get_values
    procedure::branch_vector_num_elements => base_branch_vector_num_elements
    procedure::branch_vector_set_element_index => base_branch_vector_set_element_index
    procedure::branch_vector_get_element_indices => base_branch_vector_get_element_indices
    procedure::branch_vector_get_element_values => base_branch_vector_get_element_values
    procedure::branch_vector_set_element_values => base_branch_vector_set_element_values
!
!  Base component calls
!
    procedure::branch_load => base_branch_load
    procedure(i_branch_get_xc_buf_size), deferred::branch_get_xc_buf_size
    procedure(i_branch_get_xc_buf), deferred::branch_get_xc_buf
    procedure::branch_set_mode => base_branch_set_mode
    procedure::branch_serial_write => base_branch_serial_write
!
!  Base branch calls
!
    procedure, non_overridable::branch_clear_buses
    procedure, non_overridable::branch_get_bus1_original_index
    procedure, non_overridable::branch_get_bus2_original_index
    procedure, non_overridable::branch_get_bus1_global_index
    procedure, non_overridable::branch_get_bus2_global_index
    procedure, non_overridable::branch_compare
  end type
!
  type, public :: bus_wrapper
    class(bus_component), pointer :: bus
  end type
!
  type, public :: branch_wrapper
    class(branch_component), pointer :: branch
  end type
!
!  Abstract interface for buffer calls
!
  abstract interface
!
! Return the size of the buffer needed for the data exchanges. Note that this
! must be the same size for all buses even if all buses do not require the same
! parameters. Thus, the buffer must be big enough to exchange all variables a
! bus might need, even if individual buses do use all the variables
! @param bus GridPACK bus object
! @return size of buffer in bytes
!
    integer function i_bus_get_xc_buf_size(bus)
      use, intrinsic :: iso_c_binding
      import bus_component
      implicit none
      class(bus_component), intent(in) :: bus
    end function i_bus_get_xc_buf_size
!
! Return the location of the data exchange buffer.
! @param bus GridPACK bus object
! @param pointer to exchange buffer
!
    subroutine i_bus_get_xc_buf(bus, buf)
      use, intrinsic :: iso_c_binding
      import bus_component
      implicit none
      class(bus_component), target, intent(in) :: bus
      type(C_PTR), intent(out) :: buf
    end subroutine i_bus_get_xc_buf
!
! Return the size of the buffer needed for the data exchanges. Note that this
! must be the same size for all branches even if all branches do not require the
! same parameters. Thus, the buffer must be big enough to exchange all variables
! a branch might need, even if individual branches do use all the variables
! @param branch GridPACK branch object
! @return size of buffer in bytes
!
    integer function i_branch_get_xc_buf_size(branch)
      use, intrinsic :: iso_c_binding
      import branch_component
      implicit none
      class(branch_component), intent(in) :: branch
    end function i_branch_get_xc_buf_size
!
! Return the location of the data exchange buffer.
! @param branch GridPACK branch object
! @param pointer to exchange buffer
!
    subroutine i_branch_get_xc_buf(branch, buf)
      use, intrinsic :: iso_c_binding
      import branch_component
      implicit none
      class(branch_component), target, intent(in) :: branch
      type(C_PTR), intent(out) :: buf
    end subroutine i_branch_get_xc_buf
  end interface
!
!  Interface declaration to C calls
!
  interface
!
! Get the matrix index for component, based on location of
! component in network
! @param bus GridPACK bus object
! @param idx matrix index of bus
!
    subroutine p_bus_get_mat_vec_index(bus,idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
      integer(C_INT), intent(out) :: idx
    end subroutine p_bus_get_mat_vec_index
!
! Get the number of neighboring branches/buses
! @param bus GridPACK bus object
! @return number of attached branches/buses
!
    integer(C_INT) function p_bus_get_num_neighbors(bus) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end function p_bus_get_num_neighbors
!
! Get get pointer to branch that is attached to calling bus
! @param bus GridPACK bus object
! @param idx index of neighboring branch (value is between 0 and number of
! neighbors -1)
! @return pointer to branch wrapper
!
    type(C_PTR) function p_bus_get_neighbor_branch(bus, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
      integer(C_INT), value, intent(in) :: idx
    end function p_bus_get_neighbor_branch
!
! Get get pointer to bus that is attached to calling bus via a branch
! @param bus GridPACK bus object
! @param idx index of neighboring bus (value is between 0 and number of
! neighbors -1)
! @return pointer to bus wrapper
!
    type(C_PTR) function p_bus_get_neighbor_bus(bus, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
      integer(C_INT), value, intent(in) :: idx
    end function p_bus_get_neighbor_bus
!
! Clear all pointers to neighboring branches
! @param bus GridPACK bus object
!
    subroutine p_bus_clear_branches(bus) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end subroutine p_bus_clear_branches
!
! Clear all pointers to neighboring buses
! @param bus GridPACK bus object
!
    subroutine p_bus_clear_buses(bus) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end subroutine p_bus_clear_buses
!
! Set reference bus status
! @param bus GridPACK bus object
! @param status true if bus is reference bus
!
    subroutine p_bus_set_reference_bus(bus, status) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
      logical(C_BOOL), value, intent(in) :: status
    end subroutine p_bus_set_reference_bus
!
! Get reference bus status
! @param bus GridPACK bus object
! @return true if bus is reference bus
!
    logical(C_BOOL) function p_bus_get_reference_bus(bus) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end function p_bus_get_reference_bus
!
! Get original index (from input file)
! @param bus GridPACK bus object
! @return original index from network
!
    integer(C_INT) function p_bus_get_original_index(bus) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end function p_bus_get_original_index
!
! Get global index
! @param bus GridPACK bus object
! @return global index from network
!
    integer(C_INT) function p_bus_get_global_index(bus) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end function p_bus_get_global_index
#if 0
!
! Set original index (from input file)
! @param bus GridPACK bus object
! @param idx original index from network
!
    subroutine p_bus_set_original_index(bus, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end subroutine p_bus_set_original_index
!
! Set global index
! @param bus GridPACK bus object
! @param idx global index from network
!
    subroutine p_bus_set_global_index(bus, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: bus
    end subroutine p_bus_set_global_index
#endif
!
! Get the matrix indices for component, based on location of
! component in network
! @param bus GridPACK bus object
! @param idx matrix row index of branch
! @param jdx matrix column index of branch
!
    subroutine p_branch_get_mat_vec_indices(branch,idx,jdx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: branch
      integer(C_INT), intent(out) :: idx, jdx
    end subroutine p_branch_get_mat_vec_indices
!
! Get pointer to bus that is attached to "from" end of branch
! @param branch GridPACK branch object
! @return pointer to bus wrapper
!
    type(C_PTR) function p_branch_get_bus1(branch) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: branch
    end function p_branch_get_bus1
!
! Get pointer to bus that is attached to "to" end of branch
! @param branch GridPACK branch object
! @return pointer to bus wrapper
!
    type(C_PTR) function p_branch_get_bus2(branch) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: branch
    end function p_branch_get_bus2
!
! Clear all pointers to buses at each end of branch
! @param branch GridPACK branch object
!
    subroutine p_branch_clear_buses(branch) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: branch
    end subroutine p_branch_clear_buses
!
! Get original index of "from" bus
! @param branch GridPACK branch object
! @return original index from network
!
    integer(C_INT) function p_branch_get_bus1_original_index(branch) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: branch
    end function p_branch_get_bus1_original_index
!
! Get original index of "to" bus
! @param branch GridPACK branch object
! @return original index from network
!
    integer(C_INT) function p_branch_get_bus2_original_index(branch) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: branch
    end function p_branch_get_bus2_original_index
!
! Get global index of "from" bus
! @param branch GridPACK branch object
! @return global index from network
!
    integer(C_INT) function p_branch_get_bus1_global_index(branch) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: branch
    end function p_branch_get_bus1_global_index
!
! Get global index of "to" bus
! @param branch GridPACK branch object
! @return global index from network
!
    integer(C_INT) function p_branch_get_bus2_global_index(branch) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: branch
    end function p_branch_get_bus2_global_index
  end interface
  contains
!
! Return size of matrix block on the diagonal contributed by component
! @param bus GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function base_bus_matrix_diag_size(bus, isize, jsize)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, intent(out) :: isize, jsize
    base_bus_matrix_diag_size = .false.
    return
  end function base_bus_matrix_diag_size
!
! Return the values for a diagonal matrix block. The values are returned in
! row-major order.
! @param bus GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function base_bus_matrix_diag_values(bus, values)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(out) :: values(*)
    base_bus_matrix_diag_values = .false.
    return
  end function base_bus_matrix_diag_values
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the forward direction
! @param bus GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function base_bus_matrix_forward_size(bus, isize, jsize)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, intent(out) :: isize, jsize
    base_bus_matrix_forward_size = .false.
    return
  end function base_bus_matrix_forward_size
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the reverse direction
! @param bus GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function base_bus_matrix_reverse_size(bus, isize, jsize)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_component), intent(in) :: bus
    integer, intent(out) :: isize, jsize 
    base_bus_matrix_reverse_size = .false.
    return
  end function base_bus_matrix_reverse_size
!
! Return the values of off-diagonal matrix block. The values are for the forward
! direction and are returned in row-major order.
! @param bus GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function base_bus_matrix_forward_values(bus, values)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(out) :: values(*)
    base_bus_matrix_forward_values = .false.
    return
  end function base_bus_matrix_forward_values
!
! Return the values of off-diagonal matrix block. The values are for the reverse
! direction and are returned in row-major order.
! @param bus GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function base_bus_matrix_reverse_values(bus, values)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(out) :: values(*)
    base_bus_matrix_reverse_values = .false.
    return
  end function base_bus_matrix_reverse_values
!
! Return the size of vector block contributed by component
! @param bus GridPACK bus object
! @param isize number of vector elements
! @return false if network component does not contribute vector element
!
  logical function base_bus_vector_size(bus, isize)
    implicit none
    class(bus_component), intent(in) :: bus
    integer(C_INT), intent(out) :: isize
    base_bus_vector_size = .false.
  end function base_bus_vector_size
!
! Return the values of the vector block
! @param bus GridPACK bus object
! @param values array that contains vector values
! @return false if network component does not contribute vector element
! 
  logical function base_bus_vector_values(bus, values)
    implicit none
    class(bus_component), intent(inout) :: bus
    double complex, intent(out) :: values(*)
    base_bus_vector_values = .false.
  end function base_bus_vector_values
!
! Set values in network component based on values in a vector or matrix
! @param bus GridPACK bus object
! @param values array that contains vector values
! 
  subroutine base_bus_set_values(bus, values)
    implicit none
    class(bus_component) :: bus
    double complex, intent(in) :: values(*)
  end subroutine base_bus_set_values
!
! Return the number of rows in matrix from component
! @param bus GridPACK bus object
! @return number of rows from component
!
  integer function base_bus_matrix_num_rows(bus)
    implicit none
    class(bus_component), intent(in) :: bus
  end function base_bus_matrix_num_rows
!
! Return the number of columns in matrix from component
! @param bus GridPACK bus object
! @return number of columns from component
!
  integer function base_bus_matrix_num_cols(bus)
    implicit none
    class(bus_component), intent(in) :: bus
    base_bus_matrix_num_cols = 0
  end function base_bus_matrix_num_cols
!
! Set row indices corresponding to the rows contributed by this component
! @param bus GridPACK bus object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @param idx matrix index of row irow
!
  subroutine base_bus_matrix_set_row_index(bus, irow, idx)
    implicit none
    class(bus_component) :: bus
    integer, value, intent(in) :: irow, idx
  end subroutine base_bus_matrix_set_row_index
!
! Set column indices corresponding to the columns contributed by this component
! @param bus GridPACK bus object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @param idx matrix index of column icol
!
  subroutine base_bus_matrix_set_col_index(bus, icol, idx)
    implicit none
    class(bus_component) :: bus
    integer(C_INT), value, intent(in) :: icol, idx
  end subroutine base_bus_matrix_set_col_index
!
! Get the row index for the rows contributed by this component
! @param bus GridPACK bus object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @return matrix index of row irow
!
  integer function base_bus_matrix_get_row_index(bus, irow)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, value, intent(in) :: irow
    base_bus_matrix_get_row_index = -1
  end function base_bus_matrix_get_row_index
!
! Get the col index for the columns contributed by this component
! @param bus GridPACK bus object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @return matrix index of column icol
!
  integer function base_bus_matrix_get_col_index(bus, icol)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, value, intent(in) :: icol
    base_bus_matrix_get_col_index = -1
  end function base_bus_matrix_get_col_index
!
! Return the number of matrix values contributed by this component
! @param bus GridPACK bus object
! @return number of matrix values
!
  integer function base_bus_matrix_num_values(bus)
    implicit none
    class(bus_component), intent(in) :: bus
    base_bus_matrix_num_values = 0
  end function base_bus_matrix_num_values
!
! Get a list of matrix values contributed by this components along with their
! matrix indices
! @param bus GridPACK bus object
! @param values list of matrix element values
! @param rows row indices for the matrix elements
! @param cols column indices for the matrix elements
!
  subroutine base_bus_matrix_get_values(bus, values, rows, cols)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(out) :: values(*)
    integer, intent(out) :: rows(*), cols(*)
  end subroutine base_bus_matrix_get_values
!
! Return the number of elements in vector coming from component
! @param bus GridPACK bus object
! @return number of elements contributed from component
!
  integer function base_bus_vector_num_elements(bus)
    implicit none
    class(bus_component), intent(in) :: bus
    base_bus_vector_num_elements = 0
  end function base_bus_vector_num_elements
!
! Set indices corresponding to the elements contributed by this component
! @param bus GridPACK bus object
! @param ielem index of element contributed by this component (e.g. if component
! contributes 3 elements then ielem is between 0 and 2)
! @param idx vector index of element ielem
!
  subroutine base_bus_vector_set_element_index(bus, ielem, idx)
    implicit none
    class(bus_component) :: bus
    integer, value, intent(in) :: ielem, idx
  end subroutine base_bus_vector_set_element_index
!
! Get a list of vector indices from the component
! @param bus GridPACK bus object
! @param idx list of indices the component maps onto
!
  subroutine base_bus_vector_get_element_indices(bus, idx)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, intent(out) :: idx(*)
  end subroutine base_bus_vector_get_element_indices
!
! Get a list of vector values contributed by this component and their indices
! @param bus GridPACK bus object
! @param values list of vector element values
! @param idx indices of the vector elements
!
  subroutine base_bus_vector_get_element_values(bus, values, idx)
    implicit none
    class(bus_component), intent(in) :: bus
    double complex, intent(out) :: values(*)
    integer, intent(out) :: idx(*)
  end subroutine base_bus_vector_get_element_values
!
! Transfer vector values to component
! @param bus GridPACK bus object
! @param values list of vector element values
!
  subroutine base_bus_vector_set_element_values(bus, values)
    implicit none
    class(bus_component) :: bus
    double complex, intent(out) :: values(*)
  end subroutine base_bus_vector_set_element_values
!
! Load data from DataCollection object into corresponding component.
! @param bus GridPACK bus object
! @param data DataCollection object associated with component
!
  subroutine base_bus_load(bus, data)
    implicit none
    class(bus_component) :: bus
    class(data_collection), intent(in) :: data
  end subroutine base_bus_load
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
  subroutine base_bus_set_mode(bus, mode)
    implicit none
    class(bus_component) :: bus
    integer, value, intent(in) :: mode
  end subroutine base_bus_set_mode
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
  logical function base_bus_serial_write(bus, string, bufsize, signal)
    implicit none
    class(bus_component) :: bus
    character(len=*), intent(inout) :: string
    integer, value, intent(in) :: bufsize
    character(len=*), intent(in) :: signal
  end function base_bus_serial_write
!
! Get the matrix index for component, based on location of
! component in network
! @param bus GridPACK bus object
! @param idx matrix index of bus
!
  subroutine bus_get_mat_vec_index(bus,idx)
    implicit none
    class(bus_component), intent(in) :: bus
    integer, intent(out) :: idx
    integer(C_INT) c_idx
    call p_bus_get_mat_vec_index(bus%c_this, c_idx)
    idx = c_idx
    return
  end subroutine bus_get_mat_vec_index
!
! Get the number of neighboring branches/buses
! @param bus GridPACK bus object
! @return number of attached branches/buses
!
  integer function bus_get_num_neighbors(bus)
    implicit none
    class(bus_component), intent(in) :: bus
    integer(C_INT) c_ret
    c_ret = p_bus_get_num_neighbors(bus%c_this)
    bus_get_num_neighbors = c_ret
    return
  end function bus_get_num_neighbors
!
! Clear all pointers to neighboring branches
! @param bus GridPACK bus object
!
  subroutine bus_clear_branches(bus)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_component), value, intent(in) :: bus
    call p_bus_clear_branches(bus%c_this)
    return
  end subroutine bus_clear_branches
!
! Clear all pointers to neighboring buses
! @param bus GridPACK bus object
!
  subroutine bus_clear_buses(bus)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_component), intent(in) :: bus
    call p_bus_clear_buses(bus%c_this)
    return
  end subroutine bus_clear_buses
!
! Set reference bus status
! @param bus GridPACK bus object
! @param status true if bus is reference bus
!
  subroutine bus_set_reference_bus(bus, status)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_component), value, intent(in) :: bus
    logical, value, intent(in) :: status
    logical(C_BOOL) c_status
    c_status = status
    call p_bus_set_reference_bus(bus%c_this,c_status)
    return
  end subroutine bus_set_reference_bus
!
! Get reference bus status
! @param bus GridPACK bus object
! @return true if bus is reference bus
!
  logical function bus_get_reference_bus(bus)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_component), intent(in) :: bus
    logical(C_BOOL) c_ret
    c_ret = p_bus_get_reference_bus(bus%c_this)
    bus_get_reference_bus = c_ret
    return
  end function bus_get_reference_bus
!
! Get original index (from input file)
! @param bus GridPACK bus object
! @return original index from network
!
  integer function bus_get_original_index(bus)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_component), value, intent(in) :: bus
    integer(C_INT) c_ret
    c_ret = p_bus_get_original_index(bus%c_this)
    bus_get_original_index = c_ret
    return
  end function bus_get_original_index
!
! Get global index
! @param bus GridPACK bus object
! @return global index from network
!
  integer function bus_get_global_index(bus)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_component), value, intent(in) :: bus
    integer(C_INT) c_ret
    c_ret = p_bus_get_global_index(bus%c_this)
    bus_get_global_index = c_ret
    return
  end function bus_get_global_index
!
! Compare two Fortran bus objects to see if they point to the same C object
! @param bus GridPACK bus object
! @param other second GridPACK bus object
! @return true if both buses point to the same C object, false otherwise
!
  logical function bus_compare(bus,other)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_component), value, intent(in) :: bus
    class(bus_component), value, intent(in) :: other
    integer idx1, idx2
    idx1 = bus%bus_get_global_index()
    idx2 = other%bus_get_global_index()
    if (idx1.eq.idx2) then
      bus_compare = .true.
    else
      bus_compare = .false.
    endif
    return
  end function bus_compare
#if 0
!
! Return size of matrix block on the diagonal contributed by bus
! @param ptr C pointer to GridPACK bus object
! @param isize,jsize number of rows and columns in block
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function p_bus_matrix_diag_size(ptr, isize, jsize) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer(C_INT), intent(out) :: isize, jsize
    logical f_ret
    integer f_isize, f_jsize
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_matrix_diag_size(f_isize, f_jsize)
    isize = f_isize
    jsize = f_isize
    p_bus_matrix_diag_size = f_ret
    return
  end function p_bus_matrix_diag_size
!
! Return the values for a diagonal matrix block. The values are returned in
! row-major order.
! @param ptr C pointer to GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function p_bus_matrix_diag_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    logical f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_matrix_diag_values(values)
    p_bus_matrix_diag_values = f_ret
    return
  end function p_bus_matrix_diag_values
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the forward direction
! @param ptr C pointer to GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function p_bus_matrix_forward_size(ptr, isize, jsize) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer(C_INT), intent(out) :: isize, jsize
    logical f_ret
    integer f_isize, f_jsize
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_matrix_forward_size(f_isize, f_jsize)
    isize = f_isize
    jsize = f_jsize
    p_bus_matrix_forward_size = f_ret
    return
  end function p_bus_matrix_forward_size
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the reverse direction
! @param ptr C pointer to GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function p_bus_matrix_reverse_size(ptr, isize, jsize) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer(C_INT), intent(out) :: isize, jsize 
    logical f_ret
    integer f_isize, f_jsize
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_matrix_reverse_size(f_isize, f_jsize)
    isize = f_isize
    jsize = f_jsize
    p_bus_matrix_reverse_size = f_ret
    return
  end function p_bus_matrix_reverse_size
!
! Return the values of off-diagonal matrix block. The values are for the forward
! direction and are returned in row-major order.
! @param ptr C pointer to GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function p_bus_matrix_forward_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    logical f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_matrix_forward_values(values)
    p_bus_matrix_forward_values = f_ret
    return
  end function p_bus_matrix_forward_values
!
! Return the values of off-diagonal matrix block. The values are for the reverse
! direction and are returned in row-major order.
! @param ptr C pointer to GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function p_bus_matrix_reverse_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    logical f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_matrix_reverse_values(values)
    p_bus_matrix_reverse_values = f_ret
    return
  end function p_bus_matrix_reverse_values
!
! Return the size of vector block contributed by component
! @param ptr C pointer to GridPACK bus object
! @param isize number of vector elements
! @return false if network component does not contribute vector element
!
  logical(C_BOOL) function p_bus_vector_size(ptr, isize) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer(C_INT), intent(out) :: isize
    integer f_isize
    logical f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_vector_size(f_isize)
    p_bus_vector_size = f_ret
  end function p_bus_vector_size
!
! Return the values of the vector block
! @param ptr C pointer to GridPACK bus object
! @param values array that contains vector values
! @return false if network component does not contribute vector element
! 
  logical(C_BOOL) function p_bus_vector_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    logical f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_vector_values(values)
    p_bus_vector_values = f_ret
  end function p_bus_vector_values
!
! Set values in network component based on values in a vector or matrix
! @param ptr C pointer to GridPACK bus object
! @param values array that contains vector values
! 
  subroutine p_bus_set_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(in) :: values(*)
    call C_F_POINTER(ptr,bus)
    call bus%bus%bus_set_values(values)
  end subroutine p_bus_set_values
!
! Return the number of rows in matrix from component
! @param ptr C pointer to GridPACK bus object
! @return number of rows from component
!
  integer(C_INT) function p_bus_matrix_num_rows(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_matrix_num_rows()
    p_bus_matrix_num_rows = f_ret
  end function p_bus_matrix_num_rows
!
! Return the number of columns in matrix from component
! @param ptr C pointer to GridPACK bus object
! @return number of columns from component
!
  integer(C_INT) function p_bus_matrix_num_cols(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_matrix_num_cols()
    p_bus_matrix_num_cols = f_ret
  end function p_bus_matrix_num_cols
!
! Set row indices corresponding to the rows contributed by this component
! @param ptr C pointer to GridPACK bus object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @param idx matrix index of row irow
!
  subroutine p_bus_matrix_set_row_index(ptr, irow, r_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer(C_INT), value, intent(in) :: irow, r_idx
    integer f_row, f_idx
    call C_F_POINTER(ptr,bus)
    f_row = irow
    f_idx = r_idx
    call bus%bus%bus_matrix_set_row_index(f_row, f_idx)
    return
  end subroutine p_bus_matrix_set_row_index
!
! Set column indices corresponding to the columns contributed by this component
! @param ptr C pointer to GridPACK bus object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @param idx matrix index of column icol
!
  subroutine p_bus_matrix_set_col_index(ptr, icol, c_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer(C_INT), value, intent(in) :: icol, c_idx
    integer f_col, f_idx
    call C_F_POINTER(ptr,bus)
    f_col = icol
    f_idx = c_idx
    call bus%bus%bus_matrix_set_col_index(f_col, f_idx)
    return
  end subroutine p_bus_matrix_set_col_index
!
! Get the row index for the rows contributed by this component
! @param ptr C pointer to GridPACK bus object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @return matrix index of row irow
!
  integer(C_INT) function p_bus_matrix_get_row_index(ptr, irow) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer(C_INT), value, intent(in) :: irow
    integer f_row, f_ret
    call C_F_POINTER(ptr,bus)
    f_row = irow
    f_ret = bus%bus%bus_matrix_get_row_index(f_row)
    p_bus_matrix_get_row_index = f_ret
    return
  end function p_bus_matrix_get_row_index
!
! Get the col index for the columns contributed by this component
! @param ptr C pointer to GridPACK bus object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @return matrix index of column icol
!
  integer(C_INT) function p_bus_matrix_get_col_index(ptr, icol) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer(C_INT), value, intent(in) :: icol
    integer f_col, f_ret
    call C_F_POINTER(ptr,bus)
    f_col = icol
    f_ret = bus%bus%bus_matrix_get_col_index(f_col)
    p_bus_matrix_get_col_index = f_ret
  end function p_bus_matrix_get_col_index
!
! Return the number of matrix values contributed by this component
! @param ptr C pointer to GridPACK bus object
! @return number of matrix values
!
  integer(C_INT) function p_bus_matrix_num_values(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_matrix_num_values()
    p_bus_matrix_num_values = f_ret
  end function p_bus_matrix_num_values
!
! Get a list of matrix values contributed by this components along with their
! matrix indices
! @param ptr C pointer to GridPACK bus object
! @param values list of matrix element values
! @param rows row indices for the matrix elements
! @param cols column indices for the matrix elements
!
  subroutine p_bus_matrix_get_values(ptr, values, rows, cols) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    integer(C_INT), intent(out) :: rows(*), cols(*)
    integer nval, i
    integer, allocatable :: f_rows(:), f_cols(:)
    call C_F_POINTER(ptr,bus)
    nval = bus%bus%bus_matrix_num_values()
    allocate(f_rows(nval))
    allocate(f_cols(nval))
    call bus%bus%bus_matrix_get_values(values,f_rows,f_cols)
    do i = 1, nval
      rows(i) = f_rows(i)
      cols(i) = f_cols(i)
    end do
    deallocate(f_rows)
    deallocate(f_cols)
    return
  end subroutine p_bus_matrix_get_values
!
! Return the number of elements in vector coming from component
! @param ptr C pointer to GridPACK bus object
! @return number of elements contributed from component
!
  integer(C_INT) function p_bus_vector_num_elements(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_vector_num_elements()
    p_bus_vector_num_elements = f_ret
    return
  end function p_bus_vector_num_elements
!
! Set indices corresponding to the elements contributed by this component
! @param ptr C pointer to GridPACK bus object
! @param ielem index of element contributed by this component (e.g. if component
! contributes 3 elements then ielem is between 0 and 2)
! @param idx vector index of element ielem
!
  subroutine p_bus_vector_set_element_index(ptr, ielem, e_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer(C_INT), value, intent(in) :: ielem, e_idx
    integer f_elem, f_idx
    call C_F_POINTER(ptr,bus)
    f_elem = ielem
    f_idx = e_idx
    call bus%bus%bus_vector_set_element_index(f_elem, f_idx)
    return
  end subroutine p_bus_vector_set_element_index
!
! Get a list of vector indices from the component
! @param ptr C pointer to GridPACK bus object
! @param idx list of indices the component maps onto
!
  subroutine p_bus_vector_get_element_indices(ptr, e_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer(C_INT), intent(out) :: e_idx(*)
    integer, allocatable :: f_idx(:)
    integer nval, i
    call C_F_POINTER(ptr,bus)
    nval = bus%bus%bus_vector_num_elements()
    allocate(f_idx(nval))
    call bus%bus%bus_vector_get_element_indices(f_idx)
    do i = 1, nval
      e_idx(i) = f_idx(i)
    end do
    deallocate(f_idx)
    return
  end subroutine p_bus_vector_get_element_indices
!
! Get a list of vector values contributed by this component and their indices
! @param ptr C pointer to GridPACK bus object
! @param values list of vector element values
! @param idx indices of the vector elements
!
  subroutine p_bus_vector_get_element_values(ptr, values, e_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    integer(C_INT), intent(out) :: e_idx(*)
    integer, allocatable :: f_idx(:)
    integer nval, i
    call C_F_POINTER(ptr,bus)
    nval = bus%bus%bus_vector_num_elements()
    allocate(f_idx(nval))
    call bus%bus%bus_vector_get_element_values(values, f_idx)
    do i = 1, nval
      e_idx(i) = f_idx(i)
    end do
    deallocate(f_idx)
  end subroutine p_bus_vector_get_element_values
!
! Transfer vector values to component
! @param ptr C pointer to GridPACK bus object
! @param values list of vector element values
!
  subroutine p_bus_vector_set_element_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    call C_F_POINTER(ptr,bus)
    call bus%bus%bus_vector_set_element_values(values)
    return
  end subroutine p_bus_vector_set_element_values
!
! Load data from DataCollection object into corresponding component.
! @param ptr C pointer to GridPACK bus object
! @param data DataCollection object associated with component
!
  subroutine p_bus_load(ptr, data_ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(C_PTR), value, intent(in) :: data_ptr
    type(bus_wrapper), pointer :: bus
    type(data_collection), pointer :: data
    call C_F_POINTER(ptr,bus)
    call C_F_POINTER(data_ptr,data)
    call bus%bus%bus_load(data)
    return
  end subroutine p_bus_load
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
  subroutine p_bus_set_mode(ptr, mode) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    integer(C_INT), value, intent(in) :: mode
    integer f_mode
    call C_F_POINTER(ptr,bus)
    f_mode = mode
    call bus%bus%bus_set_mode(f_mode)
    return
  end subroutine p_bus_set_mode
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
  logical(C_BOOL) function p_bus_serial_write(ptr, string, bufsize, &
      signal) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(bus_wrapper), pointer :: bus
    character(C_CHAR), intent(out) :: string(*)
    integer(C_INT), value, intent(in) :: bufsize
    character(C_CHAR), intent(in) :: signal(*)
    logical f_ret
    integer f_bufsize
    call C_F_POINTER(ptr,bus)
    f_bufsize = bufsize
    f_ret = bus%bus%bus_serial_write(string,f_bufsize,signal)
    p_bus_serial_write = f_ret
    return
  end function p_bus_serial_write
#endif
!
! Return size of matrix block on the diagonal contributed by component
! @param branch GridPACK branch object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function base_branch_matrix_diag_size(branch, isize, jsize)
    implicit none
    class(branch_component), intent(in) :: branch
    integer isize, jsize
    base_branch_matrix_diag_size = .false.
    return
  end function base_branch_matrix_diag_size
!
! Return the values for a diagonal matrix block. The values are returned in
! row-major order.
! @param branch GridPACK branch object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function base_branch_matrix_diag_values(branch, values)
    implicit none
    class(branch_component), intent(in) :: branch
    double complex, intent(out) :: values(*)
    base_branch_matrix_diag_values = .false.
    return
  end function base_branch_matrix_diag_values
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the forward direction
! @param branch GridPACK branch object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function base_branch_matrix_forward_size(branch, isize, jsize)
    implicit none
    class(branch_component), intent(in) :: branch
    integer, intent(out) :: isize, jsize
    base_branch_matrix_forward_size = .false.
    return
  end function base_branch_matrix_forward_size
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the reverse direction
! @param branch GridPACK branch object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function base_branch_matrix_reverse_size(branch, isize, jsize)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_component), intent(in) :: branch
    integer, intent(out) :: isize, jsize 
    base_branch_matrix_reverse_size = .false.
    return
  end function base_branch_matrix_reverse_size
!
! Return the values of off-diagonal matrix block. The values are for the forward
! direction and are returned in row-major order.
! @param branch GridPACK branch object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function base_branch_matrix_forward_values(branch, values)
    implicit none
    class(branch_component), intent(in) :: branch
    double complex, intent(out) :: values(*)
    base_branch_matrix_forward_values = .false.
    return
  end function base_branch_matrix_forward_values
!
! Return the values of off-diagonal matrix block. The values are for the reverse
! direction and are returned in row-major order.
! @param branch GridPACK branch object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function base_branch_matrix_reverse_values(branch, values)
    implicit none
    class(branch_component), intent(in) :: branch
    double complex, intent(out) :: values(*)
    base_branch_matrix_reverse_values = .false.
    return
  end function base_branch_matrix_reverse_values
!
! Return the size of vector block contributed by component
! @param branch GridPACK branch object
! @param isize number of vector elements
! @return false if network component does not contribute vector element
!
  logical function base_branch_vector_size(branch, isize)
    implicit none
    class(branch_component), intent(in) :: branch
    integer(C_INT), intent(out) :: isize
    base_branch_vector_size = .false.
  end function base_branch_vector_size
!
! Return the values of the vector block
! @param branch GridPACK branch object
! @param values array that contains vector values
! @return false if network component does not contribute vector element
! 
  logical function base_branch_vector_values(branch, values)
    implicit none
    class(branch_component), intent(in) :: branch
    double complex, intent(out) :: values(*)
    base_branch_vector_values = .false.
  end function base_branch_vector_values
!
! Set values in network component based on values in a vector or matrix
! @param branch GridPACK branch object
! @param values array that contains vector values
! 
  subroutine base_branch_set_values(branch, values)
    implicit none
    class(branch_component) :: branch
    double complex, intent(in) :: values(*)
  end subroutine base_branch_set_values
!
! Return the number of rows in matrix from component
! @param branch GridPACK branch object
! @return number of rows from component
!
  integer function base_branch_matrix_num_rows(branch)
    implicit none
    class(branch_component), intent(in) :: branch
  end function base_branch_matrix_num_rows
!
! Return the number of columns in matrix from component
! @param branch GridPACK branch object
! @return number of columns from component
!
  integer function base_branch_matrix_num_cols(branch)
    implicit none
    class(branch_component), intent(in) :: branch
  end function base_branch_matrix_num_cols
!
! Set row indices corresponding to the rows contributed by this component
! @param branch GridPACK branch object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @param idx matrix index of row irow
!
  subroutine base_branch_matrix_set_row_index(branch, irow, idx)
    implicit none
    class(branch_component) :: branch
    integer, value, intent(in) :: irow, idx
  end subroutine base_branch_matrix_set_row_index
!
! Set column indices corresponding to the columns contributed by this component
! @param branch GridPACK branch object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @param idx matrix index of column icol
!
  subroutine base_branch_matrix_set_col_index(branch, icol, idx)
    implicit none
    class(branch_component) :: branch
    integer(C_INT), value, intent(in) :: icol, idx
  end subroutine base_branch_matrix_set_col_index
!
! Get the row index for the rows contributed by this component
! @param branch GridPACK branch object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @return matrix index of row irow
!
  integer function base_branch_matrix_get_row_index(branch, irow)
    implicit none
    class(branch_component), intent(in) :: branch
    integer, value, intent(in) :: irow
    base_branch_matrix_get_row_index = -1
  end function base_branch_matrix_get_row_index
!
! Get the col index for the columns contributed by this component
! @param branch GridPACK branch object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @return matrix index of column icol
!
  integer function base_branch_matrix_get_col_index(branch, icol)
    implicit none
    class(branch_component), intent(in) :: branch
    integer, value, intent(in) :: icol
    base_branch_matrix_get_col_index = -1
  end function base_branch_matrix_get_col_index
!
! Return the number of matrix values contributed by this component
! @param branch GridPACK branch object
! @return number of matrix values
!
  integer function base_branch_matrix_num_values(branch)
    implicit none
    class(branch_component), intent(in) :: branch
    base_branch_matrix_num_values = 0
  end function base_branch_matrix_num_values
!
! Get a list of matrix values contributed by this components along with their
! matrix indices
! @param branch GridPACK branch object
! @param values list of matrix element values
! @param rows row indices for the matrix elements
! @param cols column indices for the matrix elements
!
  subroutine base_branch_matrix_get_values(branch, values, rows, cols)
    implicit none
    class(branch_component), intent(in) :: branch
    double complex, intent(out) :: values(*)
    integer, intent(out) :: rows(*), cols(*)
  end subroutine base_branch_matrix_get_values
!
! Return the number of elements in vector coming from component
! @param branch GridPACK branch object
! @return number of elements contributed from component
!
  integer function base_branch_vector_num_elements(branch)
    implicit none
    class(branch_component), intent(in) :: branch
    base_branch_vector_num_elements = 0
  end function base_branch_vector_num_elements
!
! Set indices corresponding to the elements contributed by this component
! @param branch GridPACK branch object
! @param ielem index of element contributed by this component (e.g. if component
! contributes 3 elements then ielem is between 0 and 2)
! @param idx vector index of element ielem
!
  subroutine base_branch_vector_set_element_index(branch, ielem, idx)
    implicit none
    class(branch_component) :: branch
    integer, value, intent(in) :: ielem, idx
  end subroutine base_branch_vector_set_element_index
!
! Get a list of vector indices from the component
! @param branch GridPACK branch object
! @param idx list of indices the component maps onto
!
  subroutine base_branch_vector_get_element_indices(branch, idx)
    implicit none
    class(branch_component), intent(in) :: branch
    integer, intent(out) :: idx(*)
  end subroutine base_branch_vector_get_element_indices
!
! Get a list of vector values contributed by this component and their indices
! @param branch GridPACK branch object
! @param values list of vector element values
! @param idx indices of the vector elements
!
  subroutine base_branch_vector_get_element_values(branch, values, idx)
    implicit none
    class(branch_component), intent(in) :: branch
    double complex, intent(out) :: values(*)
    integer, intent(out) :: idx(*)
  end subroutine base_branch_vector_get_element_values
!
! Transfer vector values to component
! @param branch GridPACK branch object
! @param values list of vector element values
!
  subroutine base_branch_vector_set_element_values(branch, values)
    implicit none
    class(branch_component) :: branch
    double complex, intent(out) :: values(*)
  end subroutine base_branch_vector_set_element_values
!
! Load data from DataCollection object into corresponding component.
! @param branch GridPACK branch object
! @param data DataCollection object associated with component
!
  subroutine base_branch_load(branch, data)
    implicit none
    class(branch_component) :: branch
    class(data_collection), intent(in) :: data
  end subroutine base_branch_load
!
! Set an internal variable that can be used to control the behavior of the
! component. This function doesn't need to be implemented, but if needed it can
! be used to change the behavior of the network in different phases of the
! calculation. For example, if a different matrix needs to be generated at
! different times, the mode can be change to get different values from the
! matrix-vector interface functions.
! @param branch GridPACK branch object
! @param mode integer indicating which mode should be used
!
  subroutine base_branch_set_mode(branch, mode)
    implicit none
    class(branch_component) :: branch
    integer, value, intent(in) :: mode
  end subroutine base_branch_set_mode
!
! Copy a string for output into buffer. The behavior of this method can be
! altered by inputting different values for the signal string
! @param branch GridPACK branch object
! @param string buffer containing string to be written to output
! @param bufsize size of string buffer in bytes
! @param signal string to control behavior of routine (e.g. what
! properties to write
! @return true if component is writing a contribution, false otherwise
!
  logical function base_branch_serial_write(branch, string, bufsize, signal)
    implicit none
    class(branch_component) :: branch
    character(len=*), intent(inout) :: string
    integer, value, intent(in) :: bufsize
    character(len=*), intent(in) :: signal
  end function base_branch_serial_write
!
! Get the matrix indices for component, based on location of
! component in network
! @param bus GridPACK bus object
! @param idx matrix row index of branch
! @param jdx matrix column index of branch
!
  subroutine branch_get_mat_vec_indices(branch, idx, jdx)
    implicit none
    class(branch_component) :: branch
    integer, intent(out) :: idx, jdx
    integer(C_INT) c_idx, c_jdx
    call p_branch_get_mat_vec_indices(branch%c_this, c_idx, c_jdx)
    idx = c_idx
    jdx = c_jdx
    return
  end subroutine branch_get_mat_vec_indices
!
! Clear all pointers to buses at each end of branch
! @param branch GridPACK branch object
!
  subroutine branch_clear_buses(branch)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_component), value, intent(in) :: branch
    call p_branch_clear_buses(branch%c_this)
    return
  end subroutine branch_clear_buses
!
! Get original index of "from" bus
! @param branch GridPACK branch object
! @return original index from network
!
  integer function branch_get_bus1_original_index(branch)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_component), value, intent(in) :: branch
    integer(C_INT) c_idx
    c_idx = p_branch_get_bus1_original_index(branch%c_this)
    branch_get_bus1_original_index = c_idx
    return
  end function branch_get_bus1_original_index
!
! Get original index of "to" bus
! @param branch GridPACK branch object
! @return original index from network
!
  integer function branch_get_bus2_original_index(branch)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_component), value, intent(in) :: branch
    integer(C_INT) c_idx
    c_idx = p_branch_get_bus2_original_index(branch%c_this)
    branch_get_bus2_original_index = c_idx
    return
  end function branch_get_bus2_original_index
!
! Get global index of "from" bus
! @param branch GridPACK branch object
! @return global index from network
!
  integer function branch_get_bus1_global_index(branch)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_component), value, intent(in) :: branch
    integer(C_INT) c_idx
    c_idx = p_branch_get_bus1_global_index(branch%c_this)
    branch_get_bus1_global_index = c_idx
    return
  end function branch_get_bus1_global_index
!
! Get global index of "to" bus
! @param branch GridPACK branch object
! @return global index from network
!
  integer function branch_get_bus2_global_index(branch)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_component), value, intent(in) :: branch
    integer(C_INT) c_idx
    c_idx = p_branch_get_bus2_global_index(branch%c_this)
    branch_get_bus2_global_index = c_idx
    return
  end function branch_get_bus2_global_index
!
! Compare two Fortran branch objects to see if they point to the same C object
! @param branch GridPACK branch object
! @param other second GridPACK branch object
! @return true if both branches point to the same C object, false otherwise
!
  logical function branch_compare(branch,other)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_component), value, intent(in) :: branch
    class(branch_component), value, intent(in) :: other
    integer idx1, idx2, o1, o2
    idx1 = branch%branch_get_bus1_global_index()
    idx2 = branch%branch_get_bus2_global_index()
    o1 = other%branch_get_bus1_global_index()
    o2 = other%branch_get_bus2_global_index()
    if (idx1.eq.o1.and.idx2.eq.o2) then
      branch_compare = .true.
    else
      branch_compare = .false.
    endif
    return
  end function branch_compare
end module gridpack_component
