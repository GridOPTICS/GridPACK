! ----------------------------------------------------------------
! file: pf_components.f90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
!
!  Fortran application component
!
module application_components
  use gridpack_component
  use gridpack_data_collection
  use, intrinsic :: iso_c_binding
  implicit none
!#define LARGE_MATRIX
!
!  Declare global variables
!
  integer, public :: JACOBIAN, YBUS, RHS, STATE, SCAL
  parameter(JACOBIAN=1, YBUS=2, RHS=3, STATE=4, SCAL=5)
!
!  Define bus and branch component type and exchange buffers
!
  private
!
!  Create application-specific types for data exchanges. Must use iso_c_binding
!  types for data declarations
!
  type, bind(c), public :: bus_xc_data
    real(C_DOUBLE) v_ang
    real(C_DOUBLE) v_mag
  end type
!
!  This is not used but leave an integer in it so that it doesn't mess up the
!  rest of the branch functions
!
  type, bind(c), public :: branch_xc_data
    integer(C_INT) dummy
  end type

  type, extends(bus_component), public :: application_bus
!
!  Application specific data elements go here
!
!
!  Y-matrix parameters
!
    double precision p_shunt_gs
    double precision p_shunt_bs
    logical p_shunt
    logical p_load
    integer p_mode
    logical p_isolated
    double precision p_ybusr, p_ybusi
!
!  Powerflow parameters
!
    double precision p_v, p_a, p_P0, p_Q0
    double precision p_theta, p_angle, p_voltage
    double precision, allocatable :: p_pg(:), p_qg(:)
    double precision, allocatable :: p_qmax(:), p_qmin(:), p_vs(:)
    integer, allocatable :: p_gstatus(:)
    integer p_ngen
    double precision p_pl, p_ql, p_sbase, p_pinj, p_qinj
    logical p_ispv, p_saveispv
    character(32), allocatable :: p_gid(:)
!
!  Required data elements are defined here
!
    type(bus_xc_data) :: xc_buf
    type(C_PTR) :: xc_ptr
    contains
!
!  Add user-defined functions here
!
    procedure :: bus_set_y_matrix !set_ybus
    procedure :: bus_get_y_matrix !get_ybus
    procedure :: bus_is_isolated
    procedure :: bus_set_isolated
    procedure :: bus_reset_voltage
    procedure :: bus_get_voltage
    procedure :: bus_get_complex_voltage
    procedure :: bus_get_phase
    procedure :: bus_is_pv
    procedure :: bus_set_is_pv
    procedure :: bus_reset_is_pv
    procedure :: bus_set_sbus
    procedure :: bus_chk_q_lim
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
!  Functions that are defined in this module but should not be modified by
!  application developer
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
!  Y-matrix parameters
!
    double precision, allocatable :: p_reactance(:)
    double precision, allocatable :: p_resistance(:)
    double precision, allocatable :: p_tap_ratio(:)
    double precision, allocatable :: p_phase_shift(:)
    double precision, allocatable :: p_charging(:)
    double precision, allocatable :: p_shunt_admt_g1(:)
    double precision, allocatable :: p_shunt_admt_b1(:)
    double precision, allocatable :: p_shunt_admt_g2(:)
    double precision, allocatable :: p_shunt_admt_b2(:)
    double precision, allocatable :: p_rate_a(:)
    logical, allocatable :: p_xform(:), p_shunt(:)
    integer p_mode
    double precision p_ybusr_frwd, p_ybusi_frwd
    double precision p_ybusr_rvrs, p_ybusi_rvrs
    logical, allocatable :: p_branch_status(:), p_switched(:)
    character(2), allocatable :: p_tag(:)
    double precision p_theta, p_sbase
    integer p_elems
    logical p_isolated, p_active
!
!  Required data elements are defined here
!
    type(branch_xc_data) :: xc_buf
    type(C_PTR) :: xc_ptr
    contains
!
!  Add user-defined function here
!
    procedure::branch_set_y_matrix
    procedure::branch_get_forward_y_matrix
    procedure::branch_get_reverse_y_matrix
    procedure::branch_get_admittance
    procedure::branch_get_transformer
    procedure::branch_get_shunt
    procedure::branch_get_line_element
    procedure::branch_get_line_tags
    procedure::branch_get_jacobian
    procedure::branch_get_pq
    procedure::branch_get_complex_power
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
!  Functions that are defined in this module but should not be modified by
!  application developer
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
    isize = 1
    jsize = 1
    bus_matrix_diag_size = .true.
    if (bus%p_mode.eq.JACOBIAN) then
      if (.not.bus%bus_is_isolated()) then
#ifdef LARGE_MATRIX
        isize = 2
        jsize = 2
        bus_matrix_diag_size = .true.
#else
        if (bus%bus_get_reference_bus()) then
          bus_matrix_diag_size = .false.
          return
        else if (bus%p_ispv) then
          isize = 1
          jsize = 1
        else
          isize = 2
          jsize = 2
        endif
        bus_matrix_diag_size = .true.
        return
#endif
      else
        bus_matrix_diag_size = .false.
      endif
    else if (bus%p_mode.eq.YBUS) then
      if (.not.bus%bus_is_isolated()) then
        bus_matrix_diag_size = .true.
        isize = 1
        jsize = 1
      else
        bus_matrix_diag_size = .false.
      endif
      return
    endif
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
    if (bus%p_mode.eq.JACOBIAN) then
      if (.not.bus%bus_is_isolated()) then
#ifdef LARGE_MATRIX
        if (.not.bus%bus_get_reference_bus()) then
          values(1) = -bus%p_qinj - bus%p_ybusi * bus%p_v * bus%p_v
          values(2) = bus%p_pinj - bus%p_ybusr * bus%p_v * bus%p_v
          values(3) = bus%p_pinj/bus%p_v + bus%p_ybusr *  bus%p_v
          values(4) = bus%p_qinj/bus%p_v - bus%p_ybusi *  bus%p_v
! Fix up elements if bus is PV bus
          if (p_ispv)
            values(2) = dcmplx(0.0d00,0.0d00)
            values(3) = dcmplx(0.0d00,0.0d00)
            values(4) = dcmplx(1.0d00,0.0d00)
          endif
        else
          values(1) = dcmplx(1.0d00,0.0d00)
          values(2) = dcmplx(0.0d00,0.0d00)
          values(3) = dcmplx(0.0d00,0.0d00)
          values(4) = dcmplx(1.0d00,0.0d00)
        endif
        bus_matrix_diag_values = .true.
        return
#else
        if (.not.bus%bus_get_reference_bus().and.(.not.bus%p_ispv)) then
          values(1) = -bus%p_qinj - bus%p_ybusi * bus%p_v * bus%p_v
          values(2) = bus%p_pinj - bus%p_ybusr * bus%p_v * bus%p_v
          values(3) = bus%p_pinj/bus%p_v + bus%p_ybusr *  bus%p_v
          values(4) = bus%p_qinj/bus%p_v - bus%p_ybusi *  bus%p_v
          bus_matrix_diag_values = .true.
        else if (.not.bus%bus_get_reference_bus().and.bus%p_ispv) then
          values(1) = -bus%p_qinj - bus%p_ybusi * bus%p_v * bus%p_v
          bus_matrix_diag_values = .true.
        else
          bus_matrix_diag_values = .false.
        endif
        return
#endif
      endif
    else if (bus%p_mode.eq.YBUS) then
      if (.not.bus%bus_is_isolated()) then
        bus_matrix_diag_values = .true.
        values(1) = dcmplx(bus%p_ybusr,bus%p_ybusi)
      else
        bus_matrix_diag_values = .false.
      endif
      return
    endif
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
    if (bus%p_mode.eq.RHS.or.bus%p_mode.eq.STATE) then
      if (.not.bus%bus_is_isolated()) then
#ifdef LARGE_MATRIX
        bus_vector_size = .true.
        isize = 2
#else
        if (bus%bus_get_reference_bus()) then
          bus_vector_size = .false.
        else if (bus%p_ispv) then
          bus_vector_size = .true.
          isize = 1
        else
          bus_vector_size = .true.
          isize = 2
        endif
#endif
        return
      else
        bus_vector_size = .false.
        return
      endif
    else if (bus%p_mode.eq.SCAL) then
      bus_vector_size = .true.
      isize = 1
    else
      bus_vector_size = .true.
      isize = 2
    endif
    return
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
    double precision retr, reti
    type(application_branch), pointer :: branch
    integer nbranch, i
    double precision pp, qq, p, q
    bus_vector_values = .true.
    if (bus%p_mode.eq.SCAL) then
      retr = bus%p_v *cos(bus%p_a)
      reti = bus%p_v *sin(bus%p_a)
      values(1) = dcmplx(retr,reti)
    else if (bus%p_mode.eq.STATE) then
      values(1) = dcmplx(bus%p_v,0.0d00)
      values(2) = dcmplx(bus%p_a,0.0d00)
    else if (bus%p_mode.eq.RHS) then
      if (.not.bus%bus_is_isolated()) then
        if (.not.bus%bus_get_reference_bus()) then
          nbranch = bus%bus_get_num_neighbors()
          pp = 0.0d00
          qq = 0.0d00
          do i = 1, nbranch
            branch => bus%bus_get_neighbor_branch(i)
            call branch%branch_get_pq(bus,p,q)
            pp = pp + p
            qq = qq + q
          end do
! Also add bus i's own Pi, Qi
          pp = pp + bus%p_v*bus%p_v*bus%p_ybusr
          qq = qq + bus%p_v*bus%p_v*(-bus%p_ybusi)
          bus%p_pinj = pp
          bus%p_qinj = qq
          pp = pp - bus%p_p0
          qq = qq - bus%p_q0
          values(1) = dcmplx(pp,0.0d00)
#ifdef LARGE_MATRIX
          if (.not.bus%p_ispv) then
            values(2) = dcmplx(qq,0.0d00)
          else
            values(2) = dcmplx(0.0d00,0.0d00)
          endif
#else
          if (.not.bus%p_ispv) then
            values(2) = dcmplx(qq,0.0d00)
          endif
#endif
          return
        else
#ifdef LARGE_MATRIX
          values(1) = dcmplx(0.0d00,0.0d00)
          values(2) = dcmplx(0.0d00,0.0d00)
#else
          bus_vector_values = .false.
#endif
          return
        endif
      else
        bus_vector_values = .false.
        return
      endif
      return
    end if
  end function bus_vector_values
!
! Check QLIM
! @return false if no violations exist, true otherwise
! @param bus GridPACK bus object
!
  logical function bus_chk_q_lim(bus)
    implicit none
    class(application_bus), intent(inout) :: bus
    type(application_branch), pointer :: branch
    double precision qmin, qmax
    double precision pp, qq, p, q
    integer i, nsize
    if (bus%p_ispv) then
      qmin = 0.0d00
      qmax = 0.0d00
      do i = 1, bus%p_ngen
        if (bus%p_gstatus(i).eq.1) then
          qmin = qmin + bus%p_qmin(i)
          qmax = qmax + bus%p_qmax(i)
        endif
      end do
      write(6,'(a,i4,a,f16.8,a,f16.8)') ' PV Check: Gen ', &
        bus%bus_get_original_index(),' p_ql = ',bus%p_ql,' QMAX = ',qmax
      nsize = bus%bus_get_num_neighbors()
      pp = 0.0d00
      qq = 0.0d00
      do i = 1, nsize
        branch => bus%bus_get_neighbor_branch(i)
        call branch%branch_get_pq(bus,p,q)
        pp = pp + p
        qq = qq + q
      end do
      write(6,'(a,i4,3(a,f16.8))') ' Gen ',bus%bus_get_original_index(), &
        ': Q = ',-qq,', p_QL = ',bus%p_ql,', Q+p_Q0 = ', &
        -qq+bus%p_ql,', QMAX = ',qmax
      if (-qq+bus%p_ql.gt.qmax) then
        write(6,'(a,i4,2(a,f16.8))') ' Gen ',bus%bus_get_original_index(), &
          ' exceeds the QMAX limit ',-qq+bus%p_ql,' vs ',qmax
        bus%p_ql = bus%p_ql+qmax
        bus%p_ispv = .false.
        do i = 1, bus%p_ngen
          bus%p_gstatus(i) = 0
        end do
        bus_chk_q_lim = .true.
        return
      else if (-qq+bus%p_ql.lt.qmin) then
        write(6,'(a,i4,2(a,f16.8))') ' Gen ',bus%bus_get_original_index(), &
          ' exceeds the QMIN limit ',-qq+bus%p_ql,' vs ',qmin
        bus%p_ql = bus%p_ql+qmin
        bus%p_ispv = .false.
        do i = 1, bus%p_ngen
          bus%p_gstatus(i) = 0
        end do
        bus_chk_q_lim = .true.
        return
      else
        bus_chk_q_lim = .false.
        return
      endif
    else
      write(6,'(a,i4,a,f16.8)') ' PQ Check: bus: ', &
        bus%bus_get_original_index(),' p_ql = ',bus%p_ql
      bus_chk_q_lim = .false.
      return
    endif
  end function bus_chk_q_lim
!
! Set values in network component based on values in a vector or matrix
! @param bus GridPACK bus object
! @param values array that contains vector values
! 
  subroutine bus_set_values(bus, values)
    implicit none
    class(application_bus) :: bus
    double complex, intent(in) :: values(*)
    double precision vt, at
    vt = bus%p_v
    at = bus%p_a
    bus%p_a = bus%p_a - real(values(1))
#ifdef LARGE_MATRIX
    bus%p_v = bus%p_v - real(values(2))
#else
    if (.not.bus%bus_is_pv()) then
      bus%p_v = bus%p_v - real(values(2))
    endif
#endif
    bus%xc_buf%v_ang = bus%p_a
    bus%xc_buf%v_mag = bus%p_v
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
    double precision pi, pg, qg, vs, qmin, qmax, sbase
    integer itype, i, ngen, gstatus, icnt
    logical ok, lgen
    character(32) id
!
!  Default values
!
    bus%p_shunt_gs = 0.0d00
    bus%p_shunt_bs = 0.0d00
    bus%p_v = 0.0d00
    bus%p_a = 0.0d00
    bus%p_theta = 0.0d00
    bus%p_angle = 0.0d00
    bus%p_voltage = 0.0d00
    bus%p_pl = 0.0d00
    bus%p_ql = 0.0d00
    bus%p_sbase = 0.0d00
    bus%p_mode = YBUS
    bus%p_isolated = .false.
!
!  Load Y-matrix parameters
!
    ok = data%get_int_value('BUS_TYPE',itype)
    ok = data%get_double_value('CASE_SBASE',sbase)
    bus%p_shunt = .true.
    bus%p_shunt = bus%p_shunt.and.data%get_double_value('BUS_SHUNT_GL', &
        bus%p_shunt_gs)
    bus%p_shunt = bus%p_shunt.and.data%get_double_value('BUS_SHUNT_BL', &
        bus%p_shunt_bs)
    bus%p_shunt_gs = bus%p_shunt_gs/sbase
    bus%p_shunt_bs = bus%p_shunt_bs/sbase
    bus%p_sbase = sbase
!
    if (itype.eq.3) then
      call bus%bus_set_reference_bus(.true.)
    else if (itype.eq.4) then
      bus%p_isolated = .true.
    endif
!
!  Load powerflow parameters
!
    ok = data%get_double_value('CASE_SBASE',bus%p_sbase)
    ok = data%get_double_value('BUS_VOLTAGE_ANG',bus%p_angle)
    ok = data%get_double_value('BUS_VOLTAGE_MAG',bus%p_voltage)
    bus%p_v = bus%p_voltage
    pi = 4.0d00*atan(1.0d00)
    bus%p_angle = bus%p_angle*pi/180.0d00
    bus%p_a = bus%p_angle
!
!  If BUS_TYPE = 2 and gstatus is 1, then bus is PV bus
!
    bus%p_ispv = .false.
!
!  add p_pg, p_qg, p_pl, p_ql, p_sbase
!
    bus%p_load = .true.
    bus%p_load = bus%p_load.and.data%get_double_value('LOAD_PL',bus%p_pl)
    bus%p_load = bus%p_load.and.data%get_double_value('LOAD_QL',bus%p_ql)
!
    bus%p_ngen = 0
    if (data%get_int_value('GENERATOR_NUMBER',ngen)) then
      bus%p_ngen = ngen
      if (ngen.gt.0) then
        allocate(bus%p_pg(ngen))
        allocate(bus%p_qg(ngen))
        allocate(bus%p_gstatus(ngen))
        allocate(bus%p_qmin(ngen))
        allocate(bus%p_qmax(ngen))
        allocate(bus%p_gid(ngen))
        icnt = 0
        do i = 1, ngen
          lgen = .true.
          lgen = lgen.and.data%get_double_indexed_value('GENERATOR_PG',pg,i)
          lgen = lgen.and.data%get_double_indexed_value('GENERATOR_QG',qg,i)
          lgen = lgen.and.data%get_double_indexed_value('GENERATOR_VS',vs,i)
          lgen = lgen.and.data%get_int_indexed_value('GENERATOR_STAT',gstatus,i)
          lgen = lgen.and.data%get_double_indexed_value('GENERATOR_QMIN',qmin,i)
          lgen = lgen.and.data%get_double_indexed_value('GENERATOR_QMAX',qmax,i)
          if (lgen) then
            icnt = icnt + 1
            bus%p_pg(icnt) = pg 
            bus%p_qg(icnt) = qg 
            bus%p_gstatus(icnt) = gstatus 
            bus%p_qmin(icnt) = qmin 
            bus%p_qmax(icnt) = qmax 
            if (gstatus.eq.1) then
              bus%p_v = vs
              if (itype.eq.2) bus%p_ispv = .true.
            endif
            id = '-1'
            ok = data%get_string_indexed_value('GENERATOR_ID',id,i)
            bus%p_gid(icnt) = id
          endif
        end do
        bus%p_ngen = icnt
      endif
    endif
    bus%p_saveisPV = bus%p_ispv
    bus%xc_buf%v_ang = bus%p_a
    bus%xc_buf%v_mag = bus%p_v
    return
  end subroutine bus_load
!
! Set values of the Y-matrix. These can be used in subseqent calculations
! @param bus GridPACK bus object
!
  subroutine bus_set_y_matrix(bus)
    implicit none
    class(application_bus), intent(inout) :: bus
    type(application_branch), pointer :: branch
    double complex ybus
    integer i,nnghbr
    ybus = dcmplx(0.0d00,0.0d00)
    nnghbr = bus%bus_get_num_neighbors()
    do i = 1, nnghbr
      branch => bus%bus_get_neighbor_branch(i)
      ybus = ybus - branch%branch_get_admittance()
      ybus = ybus - branch%branch_get_transformer(bus)
      ybus = ybus + branch%branch_get_shunt(bus)
    end do
    if (bus%p_shunt) then
      ybus = ybus + dcmplx(bus%p_shunt_gs,bus%p_shunt_bs)
    endif
    bus%p_ybusr = real(ybus)
    bus%p_ybusi = dimag(ybus)
    return
  end subroutine bus_set_y_matrix
!
! Get values of the Y-matrix.
! @param bus GridPACK bus object
! @return value of y-matrix for this bus
!
   double complex function bus_get_y_matrix(bus)
     implicit none
     class(application_bus), intent(in) :: bus
     bus_get_y_matrix = dcmplx(bus%p_ybusr,bus%p_ybusi)
     return
   end function bus_get_y_matrix
!
! Reset the voltage and phase angle to initial values
! @param bus GridPACK bus object
!
   subroutine bus_reset_voltage(bus)
     implicit none
     class(application_bus), intent(inout) :: bus
     bus%p_v = bus%p_voltage
     bus%p_a = bus%p_angle
     bus%xc_buf%v_mag = bus%p_v
     bus%xc_buf%v_ang = bus%p_a
     return
   end subroutine bus_reset_voltage
!
! Return the voltage magnitude on this bus
! @param bus GridPACK bus object
!
   double precision function bus_get_voltage(bus)
     implicit none
     class(application_bus), intent(in) :: bus
     bus_get_voltage = bus%xc_buf%v_mag
     return
   end function bus_get_voltage
!
! Return the value of the phase angle on this bus
! @param bus GridPACK bus object
!
   double precision function bus_get_phase(bus)
     implicit none
     class(application_bus), intent(in) :: bus
     bus_get_phase = bus%xc_buf%v_ang
     return
   end function bus_get_phase
!
! Return whether or not bus is a PV bus (V held fixed in powerflow equations)
! @param bus GridPACK bus object
! @return true if bus is PV bus
!
  logical function bus_is_pv(bus)
     implicit none
     class(application_bus), intent(in) :: bus
     bus_is_pv = bus%p_ispv
     return
  end function bus_is_pv
!
! Set is PV status
! @param bus GridPACK bus object
! @param status PV status (1 if bus is PV)
!
  subroutine bus_set_is_pv(bus,status)
    implicit none
    class(application_bus), intent(inout) :: bus
    integer, intent(in) :: status
    bus%p_saveispv = bus%p_ispv
    if (status.eq.0) then
      bus%p_ispv = .false.
    else
      bus%p_ispv = .true.
    endif
    bus%p_v = bus%p_voltage
    return
  end subroutine bus_set_is_pv
!
! Reset is PV status
! @param bus GridPACK bus object
! @param status PV status (1 if bus is PV)
!
  subroutine bus_reset_is_pv(bus)
    implicit none
    class(application_bus), intent(inout) :: bus
    bus%p_ispv = bus%p_saveispv
    return
  end subroutine bus_reset_is_pv
!
! Set S-bus
! @param bus GridPACK bus object
!
  subroutine bus_set_sbus(bus)
    implicit none
    class(application_bus), intent(inout) :: bus
    integer i, ngen
    double precision pg, qg
    logical usegen
    double complex sbus
    usegen = .false.
    ngen = bus%p_ngen
    do i = 1, ngen
      if (bus%p_gstatus(i).eq.1) then
        pg = pg + bus%p_pg(i)
        qg = qg + bus%p_qg(i)
        usegen = .true.
      endif
    end do
    if (usegen) then
      sbus = dcmplx((pg-bus%p_pl)/bus%p_sbase,(qg-bus%p_ql)/bus%p_sbase)
      bus%p_p0 = real(sbus)
      bus%p_q0 = dimag(sbus)
    else
      sbus = dcmplx(-bus%p_pl/bus%p_sbase,-bus%p_ql/bus%p_sbase)
      bus%p_p0 = real(sbus)
      bus%p_q0 = dimag(sbus)
    endif
    return
  end subroutine bus_set_sbus
!
! Return the complex voltage on this bus
!
  double complex function bus_get_complex_voltage(bus)
    implicit none
    class(application_bus), intent(inout) :: bus
    double complex ret
    bus%p_a = bus%xc_buf%v_ang
    bus%p_v = bus%xc_buf%v_mag
    ret = dcmplx(cos(bus%p_a),sin(bus%p_a))
    ret = bus%p_v*ret
    bus_get_complex_voltage = ret
    return
  end function bus_get_complex_voltage
!
! Return whether or not a bus is isolated
! @param bus GridPACK bus object
! @return true if bus is isolated
!
  logical function bus_is_isolated(bus)
    implicit none
    class(application_bus), intent(in) :: bus
    bus_is_isolated = bus%p_isolated
    return
  end function bus_is_isolated
!
! Change isolated status of bus
! @param bus GridPACK bus object
!
  subroutine bus_set_isolated(bus, status)
    implicit none
    class(application_bus), intent(inout) :: bus
    logical, intent(in) :: status
    bus%p_isolated = status
    return
  end subroutine bus_set_isolated
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
    bus%p_mode = mode
    return
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
    double precision pi, angle
    double complex v(2)
    integer nnghbr
    logical ok
    if (len(trim(signal)).eq.0) then
      pi = 4.0d00*atan(1.0d00) 
      angle = bus%p_a*180.0d00/pi
      write(string,'(a,i6,a,f12.6,a,f12.6)') '     ', &
        bus%bus_get_original_index(),'      ',angle,'      ',bus%p_v
      bus_serial_write = .true.
      return
    else if (trim(signal).eq.'pq') then
      ok = bus%bus_vector_values(v)
      nnghbr = bus%bus_get_num_neighbors()
      write(string,'(a,i6,a,f12.6,a,f12.6,a,i2)') '     ', &
        bus%bus_get_original_index(),'      ',real(v(1)),'      ', &
        real(v(2)),'      ',nnghbr
      bus_serial_write = .true.
      return
    endif
    bus_serial_write = .false.
    return
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
    logical ok, bus1pv, bus2pv
    class(application_bus), pointer :: bus1, bus2
    if (branch%p_mode.eq.JACOBIAN) then
      bus1 => branch%branch_get_bus1()
      bus2 => branch%branch_get_bus2()
      ok = .not.bus1%bus_get_reference_bus()
      ok = ok.and.(.not.bus2%bus_get_reference_bus())
      ok = ok.and.(.not.bus1%bus_is_isolated())
      ok = ok.and.(.not.bus2%bus_is_isolated())
      ok = ok.and.branch%p_active
      if (ok) then
#ifdef LARGE_MATRIX
        isize = 2
        jsize = 2
        branch_matrix_forward_size = .true.
#else
        bus1PV = bus1%bus_is_pv()
        bus2PV = bus2%bus_is_pv()
        if (bus1PV.and.bus2PV) then
          isize = 1
          jsize = 1
          branch_matrix_forward_size = .true.
        else if (bus1PV) then
          isize = 1
          jsize = 2
          branch_matrix_forward_size = .true.
        else if (bus2PV) then
          isize = 2
          jsize = 1
          branch_matrix_forward_size = .true.
        else
          isize = 2
          jsize = 2
          branch_matrix_forward_size = .true.
        endif
#endif
      else
        branch_matrix_forward_size = .false.
      endif
      return
    else if (branch%p_mode.eq.YBUS) then
      bus1 => branch%branch_get_bus1()
      bus2 => branch%branch_get_bus1()
      ok = .not.bus1%bus_is_isolated()
      ok = ok.and.(.not.bus2%bus_is_isolated())
      if (branch%p_active.and.ok) then
        isize = 1
        jsize = 1
        branch_matrix_forward_size = .true.
      else
        branch_matrix_forward_size = .false.
      endif
      return
    endif
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
    logical ok, bus1pv, bus2pv
    class(application_bus), pointer :: bus1, bus2
    if (branch%p_mode.eq.JACOBIAN) then
      bus1 => branch%branch_get_bus1()
      bus2 => branch%branch_get_bus2()
      ok = .not.bus1%bus_get_reference_bus()
      ok = ok.and.(.not.bus2%bus_get_reference_bus())
      ok = ok.and.(.not.bus1%bus_is_isolated())
      ok = ok.and.(.not.bus2%bus_is_isolated())
      ok = ok.and.branch%p_active
      if (ok) then
#ifdef LARGE_MATRIX
        isize = 2
        jsize = 2
        branch_matrix_reverse_size = .true.
#else
        bus1PV = bus1%bus_is_pv()
        bus2PV = bus2%bus_is_pv()
        if (bus1PV.and.bus2PV) then
          isize = 1
          jsize = 1
          branch_matrix_reverse_size = .true.
        else if (bus1PV) then
          isize = 2
          jsize = 1
          branch_matrix_reverse_size = .true.
        else if (bus2PV) then
          isize = 1
          jsize = 2
          branch_matrix_reverse_size = .true.
        else
          isize = 2
          jsize = 2
          branch_matrix_reverse_size = .true.
        endif
#endif
      else
        branch_matrix_reverse_size = .false.
      endif
      return
    else if (branch%p_mode.eq.YBUS) then
      bus1 => branch%branch_get_bus1()
      bus2 => branch%branch_get_bus1()
      ok = .not.bus1%bus_is_isolated()
      ok = ok.and.(.not.bus2%bus_is_isolated())
      if (branch%p_active.and.ok) then
        isize = 1
        jsize = 1
        branch_matrix_reverse_size = .true.
      else
        branch_matrix_reverse_size = .false.
      endif
      return
    endif
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
    class(application_bus), pointer :: bus1, bus2
    double precision t11, t12, t21, t22
    double precision cs, sn, ybusr, ybusi, v1, v2
    logical ok, bus1PV, bus2PV
    if (branch%p_mode.eq.JACOBIAN) then
      bus1 => branch%branch_get_bus1()
      bus2 => branch%branch_get_bus2()
      ok = .not.bus1%bus_get_reference_bus()
      ok = ok.and.(.not.bus2%bus_get_reference_bus())
      ok = ok.and.(.not.bus1%bus_is_isolated())
      ok = ok.and.(.not.bus2%bus_is_isolated())
      ok = ok.and.branch%p_active
      if (ok) then
        cs = cos(branch%p_theta)
        sn = sin(branch%p_theta)
        bus1PV = bus1%bus_is_pv()
        bus2PV = bus2%bus_is_pv()
        ybusr = branch%p_ybusr_frwd
        ybusi = branch%p_ybusi_frwd
        v1 = bus1%bus_get_voltage()
        v2 = bus2%bus_get_voltage()
#ifdef LARGE_MATRIX
        values(1) = ybusr*sn - ybusi*cs
        values(2) = ybusr*cs + ybusi*sn
        values(3) = ybusr*cs + ybusi*sn
        values(4) = ybusr*sn - ybusi*cs
        values(1) = v1*v2*values(1)
        values(2) = -v1*v2*values(2)
        values(3) = v1*values(3)
        values(4) = v1*values(4)
!
!  Fix up matrix if one or both buses at the end of the branch is a PV bus
!
        if (bus1PV.and.bus2PV) then
          values(2) = dcmplx(0.0d00,0.0d00)
          values(3) = dcmplx(0.0d00,0.0d00)
          values(4) = dcmplx(0.0d00,0.0d00)
        else if (bus1PV) then
          values(2) = dcmplx(0.0d00,0.0d00)
          values(4) = dcmplx(0.0d00,0.0d00)
        else if (bus2PV) then
          values(3) = dcmplx(0.0d00,0.0d00)
          values(4) = dcmplx(0.0d00,0.0d00)
        endif
#else
        if (bus1PV.and.bus2PV) then
          values(1) = ybusr*sn - ybusi*cs
          values(1) = v1*v2*values(1)
        else if (bus1PV) then
          values(1) = ybusr*sn - ybusi*cs
          values(2) = ybusr*cs + ybusi*sn
          values(1) = v1*v2*values(1)
          values(2) = v1*values(2)
        else if (bus2PV) then
          values(1) = ybusr*sn - ybusi*cs
          values(2) = ybusr*cs + ybusi*sn
          values(1) = v1*v2*values(1)
          values(2) = -v1*v2*values(2)
        else
          values(1) = ybusr*sn - ybusi*cs
          values(2) = ybusr*cs + ybusi*sn
          values(3) = ybusr*cs + ybusi*sn
          values(4) = ybusr*sn - ybusi*cs
          values(1) = v1*v2*values(1)
          values(2) = -v1*v2*values(2)
          values(3) = v1*values(3)
          values(4) = v1*values(4)
        endif
#endif
        branch_matrix_forward_values = .true.
        return
      else        
        branch_matrix_forward_values = .false.
        return
      endif
    else if (branch%p_mode.eq.YBUS) then
      bus1 => branch%branch_get_bus1()
      bus2 => branch%branch_get_bus1()
      ok = .not.bus1%bus_is_isolated()
      ok = ok.and.(.not.bus2%bus_is_isolated())
      ok = ok.and.branch%p_active
      if (ok) then
        values(1) = dcmplx(branch%p_ybusr_frwd,branch%p_ybusi_frwd)
        branch_matrix_forward_values = .true.
        return
      else
        branch_matrix_forward_values = .false.
        return
      endif
    endif
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
    class(application_bus), pointer :: bus1, bus2
    double precision t11, t12, t21, t22
    double precision cs, sn, ybusr, ybusi, v1, v2
    logical ok, bus1PV, bus2PV
    if (branch%p_mode.eq.JACOBIAN) then
      bus1 => branch%branch_get_bus1()
      bus2 => branch%branch_get_bus2()
      ok = .not.bus1%bus_get_reference_bus()
      ok = ok.and.(.not.bus2%bus_get_reference_bus())
      ok = ok.and.(.not.bus1%bus_is_isolated())
      ok = ok.and.(.not.bus2%bus_is_isolated())
      ok = ok.and.branch%p_active
      if (ok) then
        cs = cos(-branch%p_theta)
        sn = sin(-branch%p_theta)
        bus1PV = bus1%bus_is_pv()
        bus2PV = bus2%bus_is_pv()
        ybusr = branch%p_ybusr_rvrs
        ybusi = branch%p_ybusi_rvrs
        v1 = bus1%bus_get_voltage()
        v2 = bus2%bus_get_voltage()
#ifdef LARGE_MATRIX
        values(1) = ybusr*sn - ybusi*cs
        values(2) = ybusr*cs + ybusi*sn
        values(3) = ybusr*cs + ybusi*sn
        values(4) = ybusr*sn - ybusi*cs
        values(1) = v1*v2*values(1)
        values(2) = -v1*v2*values(2)
        values(3) = v2*values(3)
        values(4) = v2*values(4)
!
!  Fix up matrix if one or both buses at the end of the branch is a PV bus
!
        if (bus1PV.and.bus2PV) then
          values(2) = dcmplx(0.0d00,0.0d00)
          values(3) = dcmplx(0.0d00,0.0d00)
          values(4) = dcmplx(0.0d00,0.0d00)
        else if (bus1PV) then
          values(3) = dcmplx(0.0d00,0.0d00)
          values(4) = dcmplx(0.0d00,0.0d00)
        else if (bus2PV) then
          values(2) = dcmplx(0.0d00,0.0d00)
          values(4) = dcmplx(0.0d00,0.0d00)
        endif
#else
        if (bus1PV.and.bus2PV) then
          values(1) = ybusr*sn - ybusi*cs
          values(1) = v1*v2*values(1)
        else if (bus1PV) then
          values(1) = ybusr*sn - ybusi*cs
          values(2) = ybusr*cs + ybusi*sn
          values(1) = v1*v2*values(1)
          values(2) = -v1*v2*values(2)
        else if (bus2PV) then
          values(1) = ybusr*sn - ybusi*cs
          values(2) = ybusr*cs + ybusi*sn
          values(1) = v1*v2*values(1)
          values(2) = v2*values(2)
        else
          values(1) = ybusr*sn - ybusi*cs
          values(2) = ybusr*cs + ybusi*sn
          values(3) = ybusr*cs + ybusi*sn
          values(4) = ybusr*sn - ybusi*cs
          values(1) = v1*v2*values(1)
          values(2) = -v1*v2*values(2)
          values(3) = v2*values(3)
          values(4) = v2*values(4)
        endif
#endif
        branch_matrix_reverse_values = .true.
        return
      else        
        branch_matrix_reverse_values = .false.
        return
      endif
    else if (branch%p_mode.eq.YBUS) then
      bus1 => branch%branch_get_bus1()
      bus2 => branch%branch_get_bus1()
      ok = .not.bus1%bus_is_isolated()
      ok = ok.and.(.not.bus2%bus_is_isolated())
      ok = ok.and.branch%p_active
      if (ok) then
        values(1) = dcmplx(branch%p_ybusr_rvrs,branch%p_ybusi_rvrs)
        branch_matrix_reverse_values = .true.
        return
      else
        branch_matrix_reverse_values = .false.
        return
      endif
    endif
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
! Calculate contributions to the admittance matrix from the branches
! @param branch GridPACK branch object
!
  subroutine branch_set_y_matrix(branch)
    implicit none
    class(application_branch), intent(inout) :: branch
    class(application_bus), pointer :: bus1, bus2
    integer i, nelems
    double complex ret,a
    branch%p_ybusr_frwd = 0.0d00
    branch%p_ybusi_frwd = 0.0d00
    branch%p_ybusr_rvrs = 0.0d00
    branch%p_ybusi_rvrs = 0.0d00
    nelems = branch%p_elems
    do i = 1, nelems
      ret = dcmplx(branch%p_resistance(i),branch%p_reactance(i))
      ret = dcmplx(-1.0d00,0.0d00)/ret
      a = dcmplx(cos(branch%p_phase_shift(i)),sin(branch%p_phase_shift(i)))
      a = a*dcmplx(branch%p_tap_ratio(i),0.0d00)
      if (branch%p_switched(i)) a = conjg(a)
      if (branch%p_branch_status(i)) then
        if (branch%p_xform(i)) then
          branch%p_ybusr_frwd = branch%p_ybusr_frwd + real(ret/conjg(a))
          branch%p_ybusi_frwd = branch%p_ybusi_frwd + dimag(ret/conjg(a))
          branch%p_ybusr_rvrs = branch%p_ybusr_rvrs + real(ret/a)
          branch%p_ybusi_rvrs = branch%p_ybusi_rvrs + dimag(ret/a)
        else
          branch%p_ybusr_frwd = branch%p_ybusr_frwd + real(ret)
          branch%p_ybusi_frwd = branch%p_ybusi_frwd + dimag(ret)
          branch%p_ybusr_rvrs = branch%p_ybusr_rvrs + real(ret)
          branch%p_ybusi_rvrs = branch%p_ybusi_rvrs + dimag(ret)
        endif
      endif
    end do
    bus1 => branch%branch_get_bus1()
    bus2 => branch%branch_get_bus2()
    branch%p_theta = bus1%bus_get_phase() - bus2%bus_get_phase()
    return
  end subroutine branch_set_y_matrix
!
! Get values of y-matrix. These can be used in subsequent calculations
! @param branch GridPACK branch object
! @return complex forward y-matrix element 
!
  double complex function branch_get_forward_y_matrix(branch)
    implicit none
    class(application_branch), intent(in) :: branch
    branch_get_forward_y_matrix = dcmplx(branch%p_ybusr_frwd,branch%p_ybusi_frwd)
    return
  end function branch_get_forward_y_matrix
!
! Get values of y-matrix. These can be used in subsequent calculations
! @param branch GridPACK branch object
! @return complex reverse y-matrix element 
!
  double complex function branch_get_reverse_y_matrix(branch)
    implicit none
    class(application_branch), intent(in) :: branch
    branch_get_reverse_y_matrix = dcmplx(branch%p_ybusr_rvrs,branch%p_ybusi_rvrs)
    return
  end function branch_get_reverse_y_matrix
!
! Evaluate the complex admittance for the branch
! @param branch GridPACK branch object
! @return complex admittance from branch 
!
  double complex function branch_get_admittance(branch)
    implicit none
    class(application_branch), intent(in) :: branch
    integer i
    double complex ret, tmp
    ret = dcmplx(0.0d00,0.0d00)
    do i = 1, branch%p_elems
      tmp = dcmplx(branch%p_resistance(i),branch%p_reactance(i))
      if (.not.branch%p_xform(i).and.branch%p_branch_status(i)) then
        tmp = -1.0/tmp
      else
        tmp = dcmplx(0.0d00,0.0d00)
      endif
      ret = ret + tmp
    end do
    branch_get_admittance = ret
    return
  end function branch_get_admittance
!
! Evaluate the contribution from transformers for the branch
! @param branch GridPACK branch object
! @param bus pointer to bus at one or the other end of the branch
! @return transformer contribution to Y-matrix from branch 
!
  double complex function branch_get_transformer(branch, bus)
    implicit none
    class(application_branch), intent(in) :: branch
    class(application_bus), intent(in) :: bus
    class(application_bus), pointer :: bus1, bus2
    double complex ret, tmp, tmpb, a
    integer i
    ret = dcmplx(0.0d00,0.0d00)
    do i = 1, branch%p_elems
      tmp = dcmplx(branch%p_resistance(i),branch%p_reactance(i))
      tmpb = dcmplx(0.0d00,0.5d00*branch%p_charging(i))
      if (branch%p_xform(i).and.branch%p_branch_status(i)) then
        tmp = dcmplx(-1.0d00,0.0d00)/tmp
        tmp = tmp - tmpb
        a = dcmplx(cos(branch%p_phase_shift(i)),sin(branch%p_phase_shift(i)))
        a = a*branch%p_tap_ratio(i)
        bus1 => branch%branch_get_bus1()
        bus2 => branch%branch_get_bus2()
        if ((.not.branch%p_switched(i).and.bus%bus_compare(bus1)).or. &
            (branch%p_switched(i).and.bus%bus_compare(bus2))) then
          tmp = tmp/(conjg(a)*a)
        endif
      else
        tmp = dcmplx(0.0d00,0.0d00)
      endif
      ret = ret + tmp
    end do
    branch_get_transformer = ret
    return
  end function branch_get_transformer
!
! Evaluate the contribution from shunts for the branch
! @param branch GridPACK branch object
! @param bus pointer to bus at one or the other end of the branch
! @return shunt contribution to Y-matrix from branch 
!
  double complex function branch_get_shunt(branch, bus)
    implicit none
    class(application_branch), intent(in) :: branch
    class(application_bus), intent(in) :: bus
    class(application_bus), pointer :: bus1, bus2
    double precision retr, reti, tmpr, tmpi
    integer i
    retr = 0.0d00
    reti = 0.0d00
    do i = 1, branch%p_elems
      if (branch%p_shunt(i).and.branch%p_branch_status(i)) then
        tmpr = 0.0d00
        tmpi = 0.0d00
        if (.not.branch%p_xform(i)) then
          tmpi = 0.5d00*branch%p_charging(i)
          tmpr = 0.0d00
        endif
        bus1 => branch%branch_get_bus1()
        bus2 => branch%branch_get_bus2()
        if (bus%bus_compare(bus1)) then
          tmpr = tmpr + branch%p_shunt_admt_g1(i)
          tmpi = tmpi + branch%p_shunt_admt_b1(i)
        else if (bus%bus_compare(bus2)) then
          tmpr = tmpr + branch%p_shunt_admt_g2(i)
          tmpi = tmpi + branch%p_shunt_admt_b2(i)
        endif
      else
        tmpr = 0.0d00
        tmpi = 0.0d00
      endif
      retr = retr + tmpr
      reti = reti + tmpi
    end do
    branch_get_shunt = dcmplx(retr,reti)
    return
  end function branch_get_shunt
!
! Return contribution to Y-matrix from a specific transmission element
! @param branch GridPACK branch object
! @param tag character string for transmission element
! @param yii contribution from "from" bus
! @param yij contribution from line element
!
  subroutine branch_get_line_element(branch, tag, yii, yij)
    implicit none
    class(application_branch), intent(in) :: branch
    character(len=*), intent(in) :: tag
    double complex, intent(out) :: yii, yij
    double complex zero, flow, y, aij, bij
    integer i, idx
    idx = -1
    zero = dcmplx(0.0d00,0.0d00)
    flow = zero
    yii = zero
    yij = zero
    do i = 1, branch%p_elems
      if (trim(tag).eq.trim(branch%p_tag(i))) then
        idx = i
        exit
      endif
    end do
    if (idx.gt.0) then
      y = dcmplx(branch%p_resistance(idx),branch%p_reactance(idx))
      bij = dcmplx(0.0d00,branch%p_charging(idx))
      if (y.ne.zero) y = dcmplx(-1.0d00,0.0d00)/y
      if (branch%p_xform(idx)) then
        aij = dcmplx(cos(branch%p_phase_shift(idx)),sin(branch%p_phase_shift(idx)))
        aij = aij*branch%p_tap_ratio(idx)
        if (aij.ne.zero) then
          yij = y/conjg(aij)
          yii = -(y-0.5d00*bij)
          yii = yii/(aij*conjg(aij))
        else
          yij = y/aij
          yii = -(y-0.5d00*bij)
        endif
      else
        yij = y
        yii = -(yij-0.5d00*bij)
      endif
    endif
  end subroutine branch_get_line_element
!
! Return vector of tags
! @param branch GridPACK branch object
! @param tags vector of character strings containing tags
! @param ilen length of character strings
!
  subroutine branch_get_line_tags(branch, tags, ilen)
    implicit none
    class(application_branch), intent(in) :: branch
    integer, intent(in) :: ilen
    character(len=ilen), intent(inout) :: tags(*)
    integer i
    do i = 1, branch%p_elems
      tags(i) = trim(branch%p_tag(i))
    end do
    return
  end subroutine branch_get_line_tags
!
!
!
! Return complex power for each line element
! @param branch GridPACK branch object
! @param tag line element identifier
! @return complex power associated with line elemenet labeled by tag
!
  double complex function branch_get_complex_power(branch, tag)
    implicit none
    class(application_branch), intent(in) :: branch
    character(len=*), intent(in) :: tag
    double complex vi, vj, yii, yij, s 
    class(application_bus), pointer :: bus1, bus2
    s = dcmplx(0.0d00,0.0d00)
    bus1 => branch%branch_get_bus1()
    vi = bus1%bus_get_complex_voltage()
    bus2 => branch%branch_get_bus2()
    vj = bus2%bus_get_complex_voltage()
    call branch%branch_get_line_element(tag,yii,yij)
    s = vi*conjg(yii*vi+yij*vj)*branch%p_sbase
    branch_get_complex_power = s
    return
  end function branch_get_complex_power
!
! Return the contribution to the Jacobian in the powerflow equations comming
! from a branch
! @param branch GridPACK branch object
! @param bus pointer to bus at one end of branch
! @param values an array of 4 real doubles that holds return matrix element  
!
  subroutine branch_get_jacobian(branch, bus, values)
    implicit none
    class(application_branch), intent(in) :: branch
    class(application_bus), pointer, intent(in) :: bus
    class(application_bus), pointer :: bus1, bus2
    double precision, intent(out) :: values(*)
    double precision v, cs, sn, ybusr, ybusi
    bus1 => branch%branch_get_bus1()
    bus2 => branch%branch_get_bus2()
    if (bus%bus_compare(bus1)) then
      v = bus2%bus_get_voltage()
      cs = cos(branch%p_theta)
      sn = sin(branch%p_theta)
      ybusr = branch%p_ybusr_frwd
      ybusi = branch%p_ybusi_frwd
    else if (bus%bus_compare(bus2)) then
      v = bus1%bus_get_voltage()
      cs = cos(-branch%p_theta)
      sn = sin(-branch%p_theta)
      ybusr = branch%p_ybusr_rvrs
      ybusi = branch%p_ybusi_rvrs
    else
      write(6,'(a)') 'No bus match in branch_get_jacobian call'
    endif
    values(1) = v*(ybusr*sn - ybusi*cs)
    values(2) = -v*(ybusr*cs + ybusi*sn)
    values(3) = ybusr*cs + ybusi*sn
    values(4) = ybusr*sn - ybusi*cs
    return
  end subroutine branch_get_jacobian
!
! Return contribution to constraints
! @param branch GridPACK branch object
! @param bus pointer to bus at one end of branch
! @param p real part of constraint
! @param q imaginary part of constraint
!
  subroutine branch_get_pq(branch,bus,p,q)
    class(application_branch), intent(inout) :: branch
    class(application_bus), intent(in) :: bus
    double precision, intent(out) :: p, q
    class(application_bus), pointer :: bus1, bus2
    double precision v1, v2, cs, sn, ybusr, ybusi
    bus1 => branch%branch_get_bus1()
    bus2 => branch%branch_get_bus2()
    v1 = bus1%bus_get_voltage()
    v2 = bus2%bus_get_voltage()
    branch%p_theta = bus1%bus_get_phase() - bus2%bus_get_phase()
    if (bus%bus_compare(bus1)) then
      cs = cos(branch%p_theta)
      sn = sin(branch%p_theta)
      ybusr = branch%p_ybusr_frwd
      ybusi = branch%p_ybusi_frwd
    else if (bus%bus_compare(bus2)) then
      cs = cos(-branch%p_theta)
      sn = sin(-branch%p_theta)
      ybusr = branch%p_ybusr_rvrs
      ybusi = branch%p_ybusi_rvrs
    else
      write(6,'(a)') 'No bus match in branch_get_pq call'
    endif
    p = v1*v2*(ybusr*cs + ybusi*sn)
    q = v1*v2*(ybusr*sn - ybusi*cs)
    return
  end subroutine branch_get_pq
!
! Load data from DataCollection object into corresponding component.
! @param branch GridPACK branch object
! @param data DataCollection object associated with component
!
  subroutine branch_load(branch, data)
    implicit none
    class(application_branch) :: branch
    class(data_collection), intent(in) :: data
    logical ok, lvar, xform, shunt, rate
    double precision rvar, pi
    character(32) svar
    integer ivar, idx, nelems
!
!  Initialize some values
!
    branch%p_elems = 0
    branch%p_theta = 0.0d00
    branch%p_sbase = 0.0d00
    branch%p_mode = YBUS
!
    pi = 4.0d00*atan(1.0d00)
    ok = data%get_int_value('BRANCH_NUM_ELEMENTS',branch%p_elems)
    if (.not.ok) branch%p_elems = 0
    ok = data%get_double_value('CASE_SBASE',branch%p_sbase)
    if (.not.ok) branch%p_sbase = 0.0d00
    ok = .true.
    nelems = branch%p_elems
    branch%p_active = .false.
    if (nelems.gt.0) then
      allocate(branch%p_reactance(nelems))
      allocate(branch%p_resistance(nelems))
      allocate(branch%p_phase_shift(nelems))
      allocate(branch%p_tap_ratio(nelems))
      allocate(branch%p_tag(nelems))
      allocate(branch%p_xform(nelems))
      allocate(branch%p_branch_status(nelems))
      allocate(branch%p_switched(nelems))
      allocate(branch%p_charging(nelems))
      allocate(branch%p_shunt_admt_g1(nelems))
      allocate(branch%p_shunt_admt_b1(nelems))
      allocate(branch%p_shunt_admt_g2(nelems))
      allocate(branch%p_shunt_admt_b2(nelems))
      allocate(branch%p_shunt(nelems))
      allocate(branch%p_rate_a(nelems))
      do idx = 1, nelems
        xform = .true.
        xform = xform.and.data%get_double_indexed_value('BRANCH_X',rvar,idx)
        branch%p_reactance(idx) = rvar
        xform = xform.and.data%get_double_indexed_value('BRANCH_R',rvar,idx)
        branch%p_resistance(idx) = rvar
        ok = ok.and.data%get_double_indexed_value('BRANCH_SHIFT',rvar,idx)
        rvar = -rvar*pi/180.0d00
        branch%p_phase_shift(idx) = rvar
        ok = ok.and.data%get_double_indexed_value('BRANCH_TAP',rvar,idx)
        branch%p_tap_ratio(idx) = rvar
        ok = ok.and.data%get_string_indexed_value('BRANCH_CKT',svar,idx)
        branch%p_tag(idx) = trim(svar)
        if (rvar.ne.0.0d00) then
          branch%p_xform(idx) = xform
        else
          branch%p_xform(idx) = .false.
        endif
        ivar = 1
        xform = xform.and.data%get_int_indexed_value('BRANCH_STATUS',ivar,idx)
        if (ivar.eq.1) then
          branch%p_branch_status(idx) = .true.
        else
          branch%p_branch_status(idx) = .false.
        endif
        if (ivar.eq.1) branch%p_active = .true.
        ok = data%get_logical_indexed_value('BRANCH_SWITCHED', lvar, idx)
        if (.not.ok) lvar = .false.
        branch%p_switched(idx) = lvar
        shunt = .true.
        shunt = shunt.and.data%get_double_indexed_value('BRANCH_B',rvar,idx)
        branch%p_charging(idx) = rvar
        shunt = shunt.and.data%get_double_indexed_value('BRANCH_SHUNT_ADMTTNC_G1',rvar,idx)
        branch%p_shunt_admt_g1(idx) = rvar
        shunt = shunt.and.data%get_double_indexed_value('BRANCH_SHUNT_ADMTTNC_B1',rvar,idx)
        branch%p_shunt_admt_b1(idx) = rvar
        shunt = shunt.and.data%get_double_indexed_value('BRANCH_SHUNT_ADMTTNC_G2',rvar,idx)
        branch%p_shunt_admt_g2(idx) = rvar
        shunt = shunt.and.data%get_double_indexed_value('BRANCH_SHUNT_ADMTTNC_B2',rvar,idx)
        branch%p_shunt_admt_b2(idx) = rvar
        branch%p_shunt(idx) = shunt
        rate = data%get_double_indexed_value('BRANCH_RATING_A',rvar,idx)
        branch%p_rate_a(idx) = rvar
      end do
    endif
    return
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
    branch%p_mode = mode
    return
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
    class(application_bus), pointer :: bus1, bus2
    character(128) buf
    character(12), allocatable :: tags(:)
    integer i, ilen, oldlen
    double complex s
    logical found
    double precision ss, p, q
    bus1 => branch%branch_get_bus1()
    bus2 => branch%branch_get_bus2()
    allocate(tags(branch%p_elems))
    call branch%branch_get_line_tags(tags,12)
    if (len(trim(signal)).ne.0.and.trim(signal).eq.'flow') then
      found = .false.
      ilen = 0
      oldlen = 0
      string(1:) = ''
      do i = 1, branch%p_elems
        s = branch%branch_get_complex_power(tags(i))
        p = real(s)
        q = dimag(s)
        if (.not.branch%p_branch_status(i)) p = 0.0d00
        if (.not.branch%p_branch_status(i)) q = 0.0d00
        ss = sqrt(p**2+q**2)
        if (ss.gt.branch%p_rate_a(i).and.branch%p_rate_a(i).ne.0.0d00) then
          write(buf,'(a,i6,a,i6,a,a2,a,f12.6,a,f12.6,a,f8.2,a,f8.2,a)') &
            '     ',bus1%bus_get_original_index(),'      ', &
            bus2%bus_get_original_index(),'        ',trim(tags(i)), &
            '  ',p,'         ',q,'     ',branch%p_rate_a(i),'     ', &
            branch%p_rate_a(i)*100.0d00,'%'
          ilen = ilen + len(trim(buf))
          if ((i.lt.branch%p_elems.and.ilen+1.lt.bufsize).or. &
              (ilen.lt.bufsize)) then
            string(oldlen+1:ilen) = trim(buf)
            if (i.lt.branch%p_elems) then
              ilen = ilen + 1
              string(ilen:ilen) = new_line('a')
            endif
            oldlen = ilen
          endif
        endif
      end do
      branch_serial_write = found
    else
      ilen = 0
      oldlen = 0
      string(1:) = ''
      do i = 1, branch%p_elems
        s = branch%branch_get_complex_power(tags(i))
        p = real(s)
        q = dimag(s)
        if (.not.branch%p_branch_status(i)) p = 0.0d00
        if (.not.branch%p_branch_status(i)) q = 0.0d00
        write(buf(1:),'(a,i6,a,i6,a,a,a,f12.6,a,f12.6)') '     ', &
          bus1%bus_get_original_index(),'      ', &
          bus2%bus_get_original_index(),'     ',trim(tags(i)),'   ', &
          p,'         ',q
        ilen = ilen + len(trim(buf))
        if ((i.lt.branch%p_elems.and.ilen+1.lt.bufsize).or.(ilen.lt.bufsize)) then
          string(oldlen+1:ilen) = trim(buf)
          if (i.lt.branch%p_elems) then
            ilen = ilen + 1
            string(ilen:ilen) = new_line('a')
          endif
          oldlen = ilen
        endif
      end do
      branch_serial_write = .true.
    endif
    deallocate(tags)
    return
  end function branch_serial_write
!
!  DO NOT REMOVE THIS INCLUDE FILE. IT CONTAINS FUNCTIONS THAT MUST
!  BE INCLUDED IN THIS FILE BUT SHOULD NOT BE MODIFIED BY THE
!  APPLICATION DEVELOPER
!
  include 'component_inc.F90'
end module application_components
