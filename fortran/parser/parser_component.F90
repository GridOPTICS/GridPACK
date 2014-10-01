! ----------------------------------------------------------------
! file: parser_component.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created September 24, 2014 by Bruce Palmer
! ----------------------------------------------------------------
!
!  Fortran application component
!
module application_components
  use gridpack_component
  use gridpack_data_collection
  use, intrinsic :: iso_c_binding
  implicit none
!
!  Define bus and branch component type and exchange buffers
!
  private
!
!  Create application-specific types for data exchanges. Must use iso_c_binding
!  types for data declarations
!
  type, bind(c), public :: bus_xc_data
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
  end type

  type, bind(c), public :: branch_xc_data
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
  end type

  type, extends(bus_component), public :: application_bus
!
!  Application specific data elements go here
!
!
!  Required data elements are defined here
!
    type(bus_xc_data) :: xc_buf
    type(C_PTR) :: xc_ptr
    contains
!
!  Add user-defined function here
!
!  Matrix-vector interface calls
!
    procedure::bus_matrix_diag_size
    procedure::bus_matrix_diag_values
    procedure::bus_matrix_forward_size
    procedure::bus_matrix_reverse_size
    procedure::bus_matrix_forward_values
    procedure::bus_matrix_reverse_values
    procedure::bus_vector_size
    procedure::bus_vector_values
    procedure::bus_set_values
!
!  Generalized matrix-vector interface calls
!
    procedure::bus_matrix_num_rows
    procedure::bus_matrix_num_cols
    procedure::bus_matrix_set_row_index
    procedure::bus_matrix_set_col_index
    procedure::bus_matrix_get_row_index
    procedure::bus_matrix_get_col_index
    procedure::bus_matrix_num_values
    procedure::bus_matrix_get_values
    procedure::bus_vector_num_elements
    procedure::bus_vector_set_element_index
    procedure::bus_vector_get_element_indices
    procedure::bus_vector_get_element_values
    procedure::bus_vector_set_element_values
!
!  Base component calls
!
    procedure::bus_load
    procedure::bus_set_mode
    procedure::bus_serial_write
!
    procedure::bus_get_neighbor_branch
    procedure::bus_get_neighbor_bus
    procedure::bus_get_xc_buf_size
    procedure::bus_get_xc_buf
  end type
!
!  Define branch component type
!
  type, extends(branch_component), public :: application_branch
!
!  Application specific data elements go here
!
!
!  Required data elements are defined here
!
    type(branch_xc_data) :: xc_buf
    type(C_PTR) :: xc_ptr
    contains
!
!  Add user-defined function here
!
!  Matrix-vector interface calls
!
    procedure::branch_matrix_diag_size
    procedure::branch_matrix_diag_values
    procedure::branch_matrix_forward_size
    procedure::branch_matrix_reverse_size
    procedure::branch_matrix_forward_values
    procedure::branch_matrix_reverse_values
    procedure::branch_vector_size
    procedure::branch_vector_values
    procedure::branch_set_values
!
!  Generalized matrix-vector interface calls
!
    procedure::branch_matrix_num_rows
    procedure::branch_matrix_num_cols
    procedure::branch_matrix_set_row_index
    procedure::branch_matrix_set_col_index
    procedure::branch_matrix_get_row_index
    procedure::branch_matrix_get_col_index
    procedure::branch_matrix_num_values
    procedure::branch_matrix_get_values
    procedure::branch_vector_num_elements
    procedure::branch_vector_set_element_index
    procedure::branch_vector_get_element_indices
    procedure::branch_vector_get_element_values
    procedure::branch_vector_set_element_values
!
!  Base component calls
!
    procedure::branch_load
    procedure::branch_set_mode
    procedure::branch_serial_write
!
    procedure::branch_get_bus1
    procedure::branch_get_bus2
    procedure::branch_get_xc_buf_size
    procedure::branch_get_xc_buf
  end type
!
  type, public :: application_bus_wrapper
    type(application_bus), pointer :: bus
  end type
!
  type, public :: application_branch_wrapper
    type(application_branch), pointer :: branch
  end type
  public :: bus_cast
  public :: branch_cast
!
  contains
!
! Return size of matrix block on the diagonal contributed by component
! @param bus GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function bus_matrix_diag_size(bus, isize, jsize)
    implicit none
    class(application_bus), intent(in) :: bus
    integer, intent(out) :: isize, jsize
    bus_matrix_diag_size = .false.
    return
  end function bus_matrix_diag_size
!
! Return the values for a diagonal matrix block. The values are returned in
! row-major order.
! @param bus GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function bus_matrix_diag_values(bus, values)
    implicit none
    class(application_bus), intent(in) :: bus
    double complex, intent(out) :: values(*)
    bus_matrix_diag_values = .false.
    return
  end function bus_matrix_diag_values
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the forward direction
! @param bus GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function bus_matrix_forward_size(bus, isize, jsize)
    implicit none
    class(application_bus), intent(in) :: bus
    integer, intent(out) :: isize, jsize
    bus_matrix_forward_size = .false.
    return
  end function bus_matrix_forward_size
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the reverse direction
! @param bus GridPACK bus object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function bus_matrix_reverse_size(bus, isize, jsize)
    use, intrinsic :: iso_c_binding
    implicit none
    class(application_bus), intent(in) :: bus
    integer, intent(out) :: isize, jsize 
    bus_matrix_reverse_size = .false.
    return
  end function bus_matrix_reverse_size
!
! Return the values of off-diagonal matrix block. The values are for the forward
! direction and are returned in row-major order.
! @param bus GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function bus_matrix_forward_values(bus, values)
    implicit none
    class(application_bus), intent(in) :: bus
    double complex, intent(out) :: values(*)
    bus_matrix_forward_values = .false.
    return
  end function bus_matrix_forward_values
!
! Return the values of off-diagonal matrix block. The values are for the reverse
! direction and are returned in row-major order.
! @param bus GridPACK bus object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function bus_matrix_reverse_values(bus, values)
    implicit none
    class(application_bus), intent(in) :: bus
    double complex, intent(out) :: values(*)
    bus_matrix_reverse_values = .false.
    return
  end function bus_matrix_reverse_values
!
! Return the size of vector block contributed by component
! @param bus GridPACK bus object
! @param isize number of vector elements
! @return false if network component does not contribute vector element
!
  logical function bus_vector_size(bus, isize)
    implicit none
    class(application_bus), intent(in) :: bus
    integer(C_INT), intent(out) :: isize
    bus_vector_size = .false.
  end function bus_vector_size
!
! Return the values of the vector block
! @param bus GridPACK bus object
! @param values array that contains vector values
! @return false if network component does not contribute vector element
! 
  logical function bus_vector_values(bus, values)
    implicit none
    class(application_bus), intent(inout) :: bus
    double complex, intent(out) :: values(*)
    bus_vector_values = .false.
  end function bus_vector_values
!
! Set values in network component based on values in a vector or matrix
! @param bus GridPACK bus object
! @param values array that contains vector values
! 
  subroutine bus_set_values(bus, values)
    implicit none
    class(application_bus) :: bus
    double complex, intent(in) :: values(*)
  end subroutine bus_set_values
!
! Return the number of rows in matrix from component
! @param bus GridPACK bus object
! @return number of rows from component
!
  integer function bus_matrix_num_rows(bus)
    implicit none
    class(application_bus), intent(in) :: bus
    bus_matrix_num_rows = 0
  end function bus_matrix_num_rows
!
! Return the number of columns in matrix from component
! @param bus GridPACK bus object
! @return number of columns from component
!
  integer function bus_matrix_num_cols(bus)
    implicit none
    class(application_bus), intent(in) :: bus
    bus_matrix_num_cols = 0
  end function bus_matrix_num_cols
!
! Set row indices corresponding to the rows contributed by this component
! @param bus GridPACK bus object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @param idx matrix index of row irow
!
  subroutine bus_matrix_set_row_index(bus, irow, idx)
    implicit none
    class(application_bus) :: bus
    integer, value, intent(in) :: irow, idx
  end subroutine bus_matrix_set_row_index
!
! Set column indices corresponding to the columns contributed by this component
! @param bus GridPACK bus object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @param idx matrix index of column icol
!
  subroutine bus_matrix_set_col_index(bus, icol, idx)
    implicit none
    class(application_bus) :: bus
    integer(C_INT), value, intent(in) :: icol, idx
  end subroutine bus_matrix_set_col_index
!
! Get the row index for the rows contributed by this component
! @param bus GridPACK bus object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @return matrix index of row irow
!
  integer function bus_matrix_get_row_index(bus, irow)
    implicit none
    class(application_bus), intent(in) :: bus
    integer, value, intent(in) :: irow
    bus_matrix_get_row_index = -1
  end function bus_matrix_get_row_index
!
! Get the col index for the columns contributed by this component
! @param bus GridPACK bus object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @return matrix index of column icol
!
  integer function bus_matrix_get_col_index(bus, icol)
    implicit none
    class(application_bus), intent(in) :: bus
    integer, value, intent(in) :: icol
    bus_matrix_get_col_index = -1
  end function bus_matrix_get_col_index
!
! Return the number of matrix values contributed by this component
! @param bus GridPACK bus object
! @return number of matrix values
!
  integer function bus_matrix_num_values(bus)
    implicit none
    class(application_bus), intent(in) :: bus
    bus_matrix_num_values = 0
  end function bus_matrix_num_values
!
! Get a list of matrix values contributed by this components along with their
! matrix indices
! @param bus GridPACK bus object
! @param values list of matrix element values
! @param rows row indices for the matrix elements
! @param cols column indices for the matrix elements
!
  subroutine bus_matrix_get_values(bus, values, rows, cols)
    implicit none
    class(application_bus), intent(in) :: bus
    double complex, intent(out) :: values(*)
    integer, intent(out) :: rows(*), cols(*)
  end subroutine bus_matrix_get_values
!
! Return the number of elements in vector coming from component
! @param bus GridPACK bus object
! @return number of elements contributed from component
!
  integer function bus_vector_num_elements(bus)
    implicit none
    class(application_bus), intent(in) :: bus
    bus_vector_num_elements = 0
  end function bus_vector_num_elements
!
! Set indices corresponding to the elements contributed by this component
! @param bus GridPACK bus object
! @param ielem index of element contributed by this component (e.g. if component
! contributes 3 elements then ielem is between 0 and 2)
! @param idx vector index of element ielem
!
  subroutine bus_vector_set_element_index(bus, ielem, idx)
    implicit none
    class(application_bus) :: bus
    integer, value, intent(in) :: ielem, idx
  end subroutine bus_vector_set_element_index
!
! Get a list of vector indices from the component
! @param bus GridPACK bus object
! @param idx list of indices the component maps onto
!
  subroutine bus_vector_get_element_indices(bus, idx)
    implicit none
    class(application_bus), intent(in) :: bus
    integer, intent(out) :: idx(*)
  end subroutine bus_vector_get_element_indices
!
! Get a list of vector values contributed by this component and their indices
! @param bus GridPACK bus object
! @param values list of vector element values
! @param idx indices of the vector elements
!
  subroutine bus_vector_get_element_values(bus, values, idx)
    implicit none
    class(application_bus), intent(in) :: bus
    double complex, intent(out) :: values(*)
    integer, intent(out) :: idx(*)
  end subroutine bus_vector_get_element_values
!
! Transfer vector values to component
! @param bus GridPACK bus object
! @param values list of vector element values
!
  subroutine bus_vector_set_element_values(bus, values)
    implicit none
    class(application_bus) :: bus
    double complex, intent(out) :: values(*)
  end subroutine bus_vector_set_element_values
!
! Load data from DataCollection object into corresponding component.
! @param bus GridPACK bus object
! @param data DataCollection object associated with component
!
  subroutine bus_load(bus, data)
    implicit none
    class(application_bus) :: bus
    class(data_collection), intent(in) :: data
  end subroutine bus_load
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
  subroutine bus_set_mode(bus, mode)
    implicit none
    class(application_bus) :: bus
    integer, value, intent(in) :: mode
  end subroutine bus_set_mode
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
  logical function bus_serial_write(bus, string, bufsize, signal)
    implicit none
    class(application_bus) :: bus
    character(len=*), intent(inout) :: string
    integer, value, intent(in) :: bufsize
    character(len=*), intent(in) :: signal
  end function bus_serial_write
!
! Return size of matrix block on the diagonal contributed by component
! @param branch GridPACK branch object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function branch_matrix_diag_size(branch, isize, jsize)
    implicit none
    class(application_branch), intent(in) :: branch
    integer isize, jsize
    branch_matrix_diag_size = .false.
    return
  end function branch_matrix_diag_size
!
! Return the values for a diagonal matrix block. The values are returned in
! row-major order.
! @param branch GridPACK branch object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function branch_matrix_diag_values(branch, values)
    implicit none
    class(application_branch), intent(in) :: branch
    double complex, intent(out) :: values(*)
    branch_matrix_diag_values = .false.
    return
  end function branch_matrix_diag_values
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the forward direction
! @param branch GridPACK branch object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function branch_matrix_forward_size(branch, isize, jsize)
    implicit none
    class(application_branch), intent(in) :: branch
    integer, intent(out) :: isize, jsize
    branch_matrix_forward_size = .false.
    return
  end function branch_matrix_forward_size
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the reverse direction
! @param branch GridPACK branch object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical function branch_matrix_reverse_size(branch, isize, jsize)
    use, intrinsic :: iso_c_binding
    implicit none
    class(application_branch), intent(in) :: branch
    integer, intent(out) :: isize, jsize 
    branch_matrix_reverse_size = .false.
    return
  end function branch_matrix_reverse_size
!
! Return the values of off-diagonal matrix block. The values are for the forward
! direction and are returned in row-major order.
! @param branch GridPACK branch object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function branch_matrix_forward_values(branch, values)
    implicit none
    class(application_branch), intent(in) :: branch
    double complex, intent(out) :: values(*)
    branch_matrix_forward_values = .false.
    return
  end function branch_matrix_forward_values
!
! Return the values of off-diagonal matrix block. The values are for the reverse
! direction and are returned in row-major order.
! @param branch GridPACK branch object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical function branch_matrix_reverse_values(branch, values)
    implicit none
    class(application_branch), intent(in) :: branch
    double complex, intent(out) :: values(*)
    branch_matrix_reverse_values = .false.
    return
  end function branch_matrix_reverse_values
!
! Return the size of vector block contributed by component
! @param branch GridPACK branch object
! @param isize number of vector elements
! @return false if network component does not contribute vector element
!
  logical function branch_vector_size(branch, isize)
    implicit none
    class(application_branch), intent(in) :: branch
    integer(C_INT), intent(out) :: isize
    branch_vector_size = .false.
  end function branch_vector_size
!
! Return the values of the vector block
! @param branch GridPACK branch object
! @param values array that contains vector values
! @return false if network component does not contribute vector element
! 
  logical function branch_vector_values(branch, values)
    implicit none
    class(application_branch), intent(in) :: branch
    double complex, intent(out) :: values(*)
    branch_vector_values = .false.
  end function branch_vector_values
!
! Set values in network component based on values in a vector or matrix
! @param branch GridPACK branch object
! @param values array that contains vector values
! 
  subroutine branch_set_values(branch, values)
    implicit none
    class(application_branch) :: branch
    double complex, intent(in) :: values(*)
  end subroutine branch_set_values
!
! Return the number of rows in matrix from component
! @param branch GridPACK branch object
! @return number of rows from component
!
  integer function branch_matrix_num_rows(branch)
    implicit none
    class(application_branch), intent(in) :: branch
    branch_matrix_num_rows = 0
  end function branch_matrix_num_rows
!
! Return the number of columns in matrix from component
! @param branch GridPACK branch object
! @return number of columns from component
!
  integer function branch_matrix_num_cols(branch)
    implicit none
    class(application_branch), intent(in) :: branch
    branch_matrix_num_cols = 0
  end function branch_matrix_num_cols
!
! Set row indices corresponding to the rows contributed by this component
! @param branch GridPACK branch object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @param idx matrix index of row irow
!
  subroutine branch_matrix_set_row_index(branch, irow, idx)
    implicit none
    class(application_branch) :: branch
    integer, value, intent(in) :: irow, idx
  end subroutine branch_matrix_set_row_index
!
! Set column indices corresponding to the columns contributed by this component
! @param branch GridPACK branch object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @param idx matrix index of column icol
!
  subroutine branch_matrix_set_col_index(branch, icol, idx)
    implicit none
    class(application_branch) :: branch
    integer(C_INT), value, intent(in) :: icol, idx
  end subroutine branch_matrix_set_col_index
!
! Get the row index for the rows contributed by this component
! @param branch GridPACK branch object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @return matrix index of row irow
!
  integer function branch_matrix_get_row_index(branch, irow)
    implicit none
    class(application_branch), intent(in) :: branch
    integer, value, intent(in) :: irow
    branch_matrix_get_row_index = -1
    return
  end function branch_matrix_get_row_index
!
! Get the col index for the columns contributed by this component
! @param branch GridPACK branch object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @return matrix index of column icol
!
  integer function branch_matrix_get_col_index(branch, icol)
    implicit none
    class(application_branch), intent(in) :: branch
    integer, value, intent(in) :: icol
    branch_matrix_get_col_index = -1
    return
  end function branch_matrix_get_col_index
!
! Return the number of matrix values contributed by this component
! @param branch GridPACK branch object
! @return number of matrix values
!
  integer function branch_matrix_num_values(branch)
    implicit none
    class(application_branch), intent(in) :: branch
    branch_matrix_num_values = 0
  end function branch_matrix_num_values
!
! Get a list of matrix values contributed by this components along with their
! matrix indices
! @param branch GridPACK branch object
! @param values list of matrix element values
! @param rows row indices for the matrix elements
! @param cols column indices for the matrix elements
!
  subroutine branch_matrix_get_values(branch, values, rows, cols)
    implicit none
    class(application_branch), intent(in) :: branch
    double complex, intent(out) :: values(*)
    integer, intent(out) :: rows(*), cols(*)
  end subroutine branch_matrix_get_values
!
! Return the number of elements in vector coming from component
! @param branch GridPACK branch object
! @return number of elements contributed from component
!
  integer function branch_vector_num_elements(branch)
    implicit none
    class(application_branch), intent(in) :: branch
    branch_vector_num_elements = 0
  end function branch_vector_num_elements
!
! Set indices corresponding to the elements contributed by this component
! @param branch GridPACK branch object
! @param ielem index of element contributed by this component (e.g. if component
! contributes 3 elements then ielem is between 0 and 2)
! @param idx vector index of element ielem
!
  subroutine branch_vector_set_element_index(branch, ielem, idx)
    implicit none
    class(application_branch) :: branch
    integer, value, intent(in) :: ielem, idx
  end subroutine branch_vector_set_element_index
!
! Get a list of vector indices from the component
! @param branch GridPACK branch object
! @param idx list of indices the component maps onto
!
  subroutine branch_vector_get_element_indices(branch, idx)
    implicit none
    class(application_branch), intent(in) :: branch
    integer, intent(out) :: idx(*)
  end subroutine branch_vector_get_element_indices
!
! Get a list of vector values contributed by this component and their indices
! @param branch GridPACK branch object
! @param values list of vector element values
! @param idx indices of the vector elements
!
  subroutine branch_vector_get_element_values(branch, values, idx)
    implicit none
    class(application_branch), intent(in) :: branch
    double complex, intent(out) :: values(*)
    integer, intent(out) :: idx(*)
  end subroutine branch_vector_get_element_values
!
! Transfer vector values to component
! @param branch GridPACK branch object
! @param values list of vector element values
!
  subroutine branch_vector_set_element_values(branch, values)
    implicit none
    class(application_branch) :: branch
    double complex, intent(out) :: values(*)
  end subroutine branch_vector_set_element_values
!
! Load data from DataCollection object into corresponding component.
! @param branch GridPACK branch object
! @param data DataCollection object associated with component
!
  subroutine branch_load(branch, data)
    implicit none
    class(application_branch) :: branch
    class(data_collection), intent(in) :: data
  end subroutine branch_load
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
  subroutine branch_set_mode(branch, mode)
    implicit none
    class(application_branch) :: branch
    integer, value, intent(in) :: mode
  end subroutine branch_set_mode
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
  logical function branch_serial_write(branch, string, bufsize, signal)
    implicit none
    class(application_branch) :: branch
    character(len=*), intent(inout) :: string
    integer, value, intent(in) :: bufsize
    character(len=*), intent(in) :: signal
  end function branch_serial_write
!
!  DO NOT REMOVE THIS INCLUDE FILE. IT CONTAINS FUNCTIONS THAT MUST
!  BE INCLUDED IN THIS FILE BUT SHOULD NOT BE MODIFIED BY THE
!  APPLICATION DEVELOPER
!
  include 'component_inc.F90'
end module application_components
