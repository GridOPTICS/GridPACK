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
  integer JACOBIAN, YBUS, RHS, STATE, SCAL
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
    integer p_mode
    logical p_isolated
    double precision p_ybusr, p_ybusi
!
!  Powerflow parameters
!
    double precision p_v, p_a, p_P0, p_Q0
    double precision p_angle, p_voltage
    double precision, allocatable :: p_pg(:), p_qg(:)
    double precision, allocatable :: p_qmax(:), p_qmin(:), p_vs(:)
    integer, allocatable :: p_gstatus(:)
    integer p_ngen
    double precision p_pl, p_ql, p_sbase, p_pinj, p_qinj
    logical p_ispv, p_saveispv
    character(32) p_gid(:)
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
    procedure :: bus_get_gen_status
    procedure :: bus_get_generators
    procedure :: bus_is_pv
    procedure :: bus_set_voltage
    procedure :: bus_set_phase
    procedure :: bus_set_gen_status
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
    if (p_mode.eq.JACOBIAN) then
      if (.not.bus%is_isolated()) then
#ifdef LARGE_MATRIX
        isize = 2
        jsize = 2
        bus_matrix_diag_size = .true.
#else
        if (bus%get_reference_bus()) then
          bus_matrix_diag_size = .false.
          return
        else if (p_ispv) then
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
    else if (p_mode.eq.YBUS) then
      if (.not.bus%is_isolated()) then
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
    if (p_mode.eq.JACOBIAN) then
      if (.not.bus%is_isolated()) then
#ifdef LARGE_MATRIX
        if (.not.get_reference_bus()) then
          values(1) = -p_qinj - p_ybusi * p_v * p_v
          values(2) = p_pinj - p_ybusr * p_v * p_v
          values(3) = p_pinj/p_v + p_ybusr *  p_v
          values(4) = p_qinj/p_v - p_ybusi *  p_v
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
        if (.not.get_reference_bus().and.(.not.p_ispv)) then
          values(1) = -p_qinj - p_ybusi * p_v * p_v
          values(2) = p_pinj - p_ybusr * p_v * p_v
          values(3) = p_pinj/p_v + p_ybusr *  p_v
          values(4) = p_qinj/p_v - p_ybusi *  p_v
          bus_matrix_diag_values = .true.
        else if (.not.get_reference_bus().and.p_ispv) then
          values(1) = -p_qinj - p_ybusi * p_v * p_v
          bus_matrix_diag_values = .true.
        else
          bus_matrix_diag_values = .false.
        endif
        return
#endif
      endif
    else if (p_mode.eq.YBUS) then
      if (.not.bus%is_isolated()) then
        bus_matrix_diag_values = .true.
        values(1) = dcmplx(p_ybusr,p_ybusi)
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
    if (p_mode.eq.RHS.or.p_mode.eq.STATE) then
      if (.not.is_isolated()) then
#ifdef LARGE_MATRIX
        bus_vector_size = .true.
        size = 2
#else
        if (get_reference_bus()) then
          bus_vector_size = .false.
        else if (p_ispv) then
          bus_vector_size = .true.
          size = 1
        else
          bus_vector_size = .true.
          size = 2
        endif
#endif
        return
      else
        bus_vector_size = .false.
        return
      endif
    else if (p_mode.eq.SCAL) then
      bus_vector_size = .true.
      size = 1
    else
      bus_vector_size = .true.
      size = 2
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
    class(application_bus), intent(in) :: bus
    double complex, intent(out) :: values(*)
    double precision retr, reti
    type(application_branch) branch
    integer nbranch, i
    double precision pp, qq, p, q
    bus_vector_values = .true.
    if (p_mode.eq.SCAL) then
      retr = p_v *cos(p_a)
      reti = p_v *sin(p_a)
      values(1) = dcmplx(retr,reti)
    else if (p_mode.eq.STATE) then
      values(1) = dcmplx(p_v,0.0d00)
      values(2) = dcmplx(p_a,0.0d00)
    else if (p_mode.eq.RHS) then
      if (.not.is_isolated()) then
        if (.not.get_reference_bus()) then
          nbranch = bus%bus_get_num_neighbors()
          pp = 0.0d00
          qq = 0.0d00
          do i = 1, nbranch
            branch = get_neighbor_branch(i)
            branch%getPQ(bus,p,q)
            pp = pp + p
            qq = qq + q
          end do
! Also add bus i's own Pi, Qi
          pp = pp + p_v*p_v*p_ybusr
          qq = qq + p_v*p_v*(-p_ybusi)
          p_pinj = pp
          p_qinj = qq
          pp = pp - p_p0
          qq = qq - p_q0
          values(1) = dcmplx(pp,0.0d00)
#ifdef LARGE_MATRIX
          if (.not.p_ispv) then
            values(2) = dcmplx(qq,0.0d00)
          else
            values(2) = dcmplx(0.0d00,0.0d00)
          endif
#else
          if (.not.p_ispv) then
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
  end function bus_vector_values
!
! Check QLIM
! @return false if no violations exist, true otherwise
! @param bus GridPACK bus object
!
  logical function bus_chk_q_lim(bus)
    implicit none
    class(application_bus) :: bus
    class(application_branch), pointer :: branch_ptr
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
        branch_ptr => bus%bus_get_neighbor_branch(i)
        call branch%get_pq(bus,p,q)
        pp = pp + p
        qq = qq + q
      end do
      write(6,'(a,i4,3(a,f16.8))') ' Gen ',bus%bus_get_original_index(), &
        ': Q = ',-qq,', p_QL = ',p_ql,', Q+p_Q0 = ',-qq+bus%p_ql,', QMAX = ',qmax
      if (-qq+bus%p_ql.gt.qmax) then
        write(6,'(a,i4,2(a,f16.8))') ' Gen ',bus%bus_get_original_index(), &
          ' exceeds the QMAX limit ',-qq+bus%p_ql,' vs ',qmax
        bus%p_ql = bus%p_ql+qmax
        bus%p_ispv = 0
        do i = 1, bus%p_ngen
          bus%p_gstatus(i) = 0
        end do
        bus_chk_q_lim = .true.
        return
      else if (-qq+bus%p_ql.lt.qmin) then
        write(6,'(a,i4,2(a,f16.8))') ' Gen ',bus%bus_get_original_index(), &
          ' exceeds the QMIN limit ',-qq+bus%p_ql,' vs ',qmin
        bus%p_ql = bus%p_ql+qmin
        bus%p_ispv = 0
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
    if (.not.bus%is_pv) then
      bus%p_v = bus%p_v - real(values(2))
    endif
#endif
    bus%xc_data%v_ang = bus%p_a
    bus%xc_data%v_mag = bus%p_v
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
    p_isolated = .false.
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
      call bus%set_reference_bus(.true.)
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
    if (data%get_int_value('GENERATOR_NUMBER',bus%p_ngen)) then
      ngen = bus%p_ngen
      allocate(bus%p_pg(ngen))
      allocate(bus%p_qg(ngen))
      allocate(bus%p_gstatus(ngen))
      allocate(bus%p_qmin(ngen))
      allocate(bus%p_qmax(ngen))
      allocate(bus%p_gid(ngen))
      icnt = 0
      do i = 1, ngen
        lgen = .true.
        lgen = lgen.and.data%get_double_value('GENERATOR_PG',pg,i)
        lgen = lgen.and.data%get_double_value('GENERATOR_QG',qg,i)
        lgen = lgen.and.data%get_double_value('GENERATOR_VS',vs,i)
        lgen = lgen.and.data%get_int_value('GENERATOR_STAT',gstatus,i)
        lgen = lgen.and.data%get_double_value('GENERATOR_QMIN',qmin,i)
        lgen = lgen.and.data%get_double_value('GENERATOR_QMAX',qmax,i)
        if (lgen) then
          icnt = icnt + 1
          bus%p_pg(icnt) = pg 
          bus%p_qg(icnt) = qg 
          bus%p_gstatus(icnt) = gstatus 
          bus%p_qmin(icnt) = qmin 
          bus%p_qmax(icnt) = qmax 
          if (gstatus.eq.1) then
            bus%p_v = vs
            if (itype.eq.2) then p_ispv = .true.
          endif
          id = '-1'
          ok = data%get_string_value('GENERATOR_ID',id,i)
          bus%p_gid = id
        endif
      end do
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
    class(application_bus) :: bus
    class(application_branch), pointer :: branch
    double complex ybus
    integer i,nnghbr
    ybus = dcmplx(0.0d00,0.0d00)
    do i = 1, nnghbr
      branch => bus%bus_get_neighbor_branch(i)
      ybus = ybus + branch%get_admittance()
      ybus = ybus + branch%get_transformer(bus)
      ybus = ybus + branch%get_shunt(bus)
    end do
    if (bus%p_shunt) then
      ybus = ybus + dcmplx(bus%p_shunt_gs,bus%p_shunt_bs)
    endif
    bus%p_ybusr = real(ybus)
    bus%p_ybusi =imag(ybus)
    return
  end subroutine bus_set_y_matrix
!
! Get values of the Y-matrix.
! @param bus GridPACK bus object
! @return value of y-matrix for this bus
!
   double complex function bus_get_ymatrix(bus)
     implicit none
     class(application_bus) :: bus
     bus_get_ymatrix = dcmplx(bus%p_ybusr,bus%p_ybusi)
     return
   end function bus_get_ymatrix
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
!  DO NOT EDIT ANYTHING BELOW THIS LINE. THESE FUNCTIONS MUST BE INCLUDED IN
!  THIS FILE BUT SHOULD NOT BE MODIFIED BY THE APPLICATION DEVELOPER
!
! This function converts a C-style pointer coming from the network to a fortran
! bus pointer
! @param ptr C-pointer to fortran wrapper
! @return fortran pointer bus object
!
  function bus_cast(ptr) result(bus)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_bus), pointer :: bus
    type(application_bus_wrapper), pointer :: wbus
    call C_F_POINTER(ptr,wbus)
    bus => wbus%bus
    return
  end function bus_cast
!
! Return the size of the buffer needed for the data exchanges. Note that this
! must be the same size for all buses even if all buses do not require the same
! parameters. Thus, the buffer must be big enough to exchange all variables a
! bus might need, even if individual buses do use all the variables
! @param bus GridPACK bus object
! @return size of buffer in bytes
!
  integer function bus_get_xc_buf_size(bus)
    implicit none
    class(application_bus), intent(in) :: bus
    bus_get_xc_buf_size = c_sizeof(bus%xc_buf)
  end function bus_get_xc_buf_size
!
! Return the location of the data exchange buffer.
! @param bus GridPACK bus object
! @param pointer to exchange buffer
!
  subroutine bus_get_xc_buf(bus, buf)
    implicit none
    class(application_bus), target, intent(in) :: bus
    type(C_PTR), intent(out) :: buf
    buf = c_loc(bus%xc_buf)
  end subroutine bus_get_xc_buf
!
! Get get pointer to branch that is attached to calling bus
! @param bus GridPACK bus object
! @param idx index of neighboring branch (value is between 0 and number of
! neighbors -1)
! @return pointer to branch
!
  function bus_get_neighbor_branch(bus, idx) result(branch_ptr)
    use, intrinsic :: iso_c_binding
    implicit none
    class(application_bus), value, intent(in) :: bus
    integer, value, intent(in) :: idx
    type(application_branch), pointer :: branch_ptr
    type(application_branch_wrapper), pointer :: wbranch
    type(C_PTR) ptr
    integer(C_INT) c_idx
    c_idx = idx
    ptr = p_bus_get_neighbor_branch(bus%c_this,c_idx)
    call C_F_POINTER(ptr,wbranch)
    branch_ptr => wbranch%branch
    return
  end function bus_get_neighbor_branch
!
! Get get pointer to bus that is attached to calling bus via a branch
! @param bus GridPACK bus object
! @param idx index of neighboring bus (value is between 0 and number of
! neighbors -1)
! @return pointer to bus
!
  function bus_get_neighbor_bus(bus, idx) result(bus_ptr)
    use, intrinsic :: iso_c_binding
    implicit none
    class(application_bus), value, intent(in) :: bus
    integer, value, intent(in) :: idx
    class(application_bus), pointer :: bus_ptr
    type(application_bus_wrapper), pointer :: wbus
    type(C_PTR) ptr
    integer(C_INT) c_idx
    c_idx = idx
    ptr = p_bus_get_neighbor_bus(bus%c_this,c_idx)
    call C_F_POINTER(ptr,wbus)
    bus_ptr => wbus%bus
    return
  end function bus_get_neighbor_bus
!
! Allocate Fortran bus wrapper objects and return a C pointer to that
! object
! @return pointer to bus wrapper object
!
  type(C_PTR) function bus_allocate(this) bind(c)
    implicit none
    type(C_PTR), value, intent(in) :: this
    type(application_bus_wrapper), pointer :: wbus
    type(application_bus), pointer :: bus
    allocate(bus)
    allocate(wbus)
    bus%c_this = this
    wbus%bus => bus
    bus_allocate = C_LOC(wbus)
    bus%xc_ptr = C_LOC(bus%xc_buf)
    return
  end function bus_allocate
!
! Deallocate Fortran bus objects
!
  subroutine bus_deallocate(ptr) bind(c)
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_bus_wrapper), pointer :: wbus
    call C_F_POINTER(ptr,wbus)
    if (associated(wbus)) then
      deallocate(wbus%bus)
      deallocate(wbus)
    endif
    return
  end subroutine bus_deallocate
!
! Return the size of the buffer needed for the data exchanges. Note that this
! must be the same size for all buses even if all buses do not require the same
! parameters. Thus, the buffer must be big enough to exchange all variables a
! bus might need, even if individual buses do use all the variables
! @param ptr C pointer to GridPACK bus object
! @return size of buffer in bytes
!
  integer(C_INT) function p_bus_get_xc_buf_size(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_bus_wrapper), pointer :: bus
    integer f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_get_xc_buf_size()
    p_bus_get_xc_buf_size = f_ret
    return
  end function p_bus_get_xc_buf_size
!
! Assign the location of the data exchange buffer.
! @param ptr C pointer to GridPACK bus object
! @param buf C pointer to exchange buffer
!
  subroutine p_bus_get_xc_buf(ptr, buf) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_bus_wrapper), pointer :: bus
    type(C_PTR), intent(out) :: buf
    call C_F_POINTER(ptr,bus)
    call bus%bus%bus_get_xc_buf(buf)
    return
  end subroutine p_bus_get_xc_buf
#if 1
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
    integer(C_INT), intent(out) :: isize
    integer f_isize
    logical f_ret
    call C_F_POINTER(ptr,bus)
    f_ret = bus%bus%bus_vector_size(f_isize)
    isize = f_isize
    p_bus_vector_size = f_ret
    return
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
    type(application_bus_wrapper), pointer :: bus
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
      signal, signal_len, write_len) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_bus_wrapper), pointer :: bus
    character(C_CHAR), intent(out) :: string(*)
    integer(C_INT), value, intent(in) :: bufsize, signal_len
    character(C_CHAR), intent(in) :: signal(*)
    integer(C_INT), intent(out) :: write_len
    character(len=bufsize) fstring
    character(len=signal_len) fsignal
    logical f_ret
    integer f_bufsize, slen, i, fslen
    call C_F_POINTER(ptr,bus)
    f_bufsize = bufsize
    fslen = signal_len
    do i = 1, fslen
      fsignal(i:i) = signal(i)
    end do
    f_ret = bus%bus%bus_serial_write(fstring,f_bufsize,fsignal)
    slen = len(trim(fstring))
    slen = min(bufsize-1,slen)
    i = 1
    do while (i.le.slen.and.f_ret)
      string(i) = fstring(i:i)
      i = i + 1
    end do
    string(slen+1) = C_NULL_CHAR
    write_len = slen
    p_bus_serial_write = f_ret
    return
  end function p_bus_serial_write
#endif
!
! This function converts a C-style pointer coming from the network to a fortran
! branch pointer
! @param ptr C-pointer to fortran wrapper
! @return fortran pointer branch object
!
  function branch_cast(ptr) result(branch)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch), pointer :: branch
    type(application_branch_wrapper), pointer :: wbranch
    call C_F_POINTER(ptr,wbranch)
    branch => wbranch%branch
    return
  end function branch_cast
!
! Return the size of the buffer needed for the data exchanges. Note that this
! must be the same size for all branches even if all branches do not require the same
! parameters. Thus, the buffer must be big enough to exchange all variables a
! branch might need, even if individual branches do use all the variables
! @param branch GridPACK branch object
! @return size of buffer in bytes
!
  integer function branch_get_xc_buf_size(branch)
    implicit none
    class(application_branch), intent(in) :: branch
    branch_get_xc_buf_size = c_sizeof(branch%xc_buf)
  end function branch_get_xc_buf_size
!
! Return the location of the data exchange buffer.
! @param branch GridPACK branch object
! @param pointer to exchange buffer
!
  subroutine branch_get_xc_buf(branch, buf)
    implicit none
    class(application_branch), target, intent(in) :: branch
    type(C_PTR), intent(out) :: buf
    buf = c_loc(branch%xc_buf)
  end subroutine branch_get_xc_buf
!
! Get pointer to bus that is attached to "from" end of branch
! @param branch GridPACK branch object
! @return pointer to bus wrapper
!
  function branch_get_bus1(branch) result(bus_ptr)
    use, intrinsic :: iso_c_binding
    implicit none
    class(application_branch), value, intent(in) :: branch
    class(application_bus), pointer :: bus_ptr
    type(application_bus_wrapper), pointer :: wbus
    type(C_PTR) ptr
    ptr = p_branch_get_bus1(branch%c_this)
    call C_F_POINTER(ptr,wbus)
    bus_ptr => wbus%bus
    return
  end function branch_get_bus1
!
! Get pointer to bus that is attached to "to" end of branch
! @param branch GridPACK branch object
! @return pointer to bus wrapper
!
  function branch_get_bus2(branch) result(bus_ptr)
    use, intrinsic :: iso_c_binding
    implicit none
    class(application_branch), value, intent(in) :: branch
    class(application_bus), pointer :: bus_ptr
    type(application_bus_wrapper), pointer :: wbus
    type(C_PTR) ptr
    ptr = p_branch_get_bus2(branch%c_this)
    call C_F_POINTER(ptr,wbus)
    bus_ptr => wbus%bus
    return
  end function branch_get_bus2
!
! Allocate Fortran branch wrapper objects and return a C pointer to that
! object
! @return pointer to branch wrapper object
!
  type(C_PTR) function branch_allocate(this) bind(c)
    implicit none
    type(C_PTR), value, intent(in) :: this
    type(application_branch_wrapper), pointer :: wbranch
    type(application_branch), pointer :: branch
    allocate(branch)
    allocate(wbranch)
    branch%c_this = this
    wbranch%branch => branch
    branch_allocate = C_LOC(wbranch)
    branch%xc_ptr = C_LOC(branch%xc_buf)
    return
  end function branch_allocate
!
! Deallocate Fortran branch objects
!
  subroutine branch_deallocate(ptr) bind(c)
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: wbranch
    call C_F_POINTER(ptr,wbranch)
    if (associated(wbranch)) then
      deallocate(wbranch%branch)
      deallocate(wbranch)
    endif
    return
  end subroutine branch_deallocate
!
! Return the size of the buffer needed for the data exchanges. Note that this
! must be the same size for all branches even if all branches do not require the
! same parameters. Thus, the buffer must be big enough to exchange all variables
! a branch might need, even if individual branches do use all the variables
! @param ptr C pointer to GridPACK branch object
! @return size of buffer in bytes
!
  integer(C_INT) function p_branch_get_xc_buf_size(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer f_ret
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_get_xc_buf_size()
    p_branch_get_xc_buf_size = f_ret
    return
  end function p_branch_get_xc_buf_size
!
! Assign the location of the data exchange buffer.
! @param ptr C pointer to GridPACK branch object
! @param buf C pointer to exchange buffer
!
  subroutine p_branch_get_xc_buf(ptr, buf) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    type(C_PTR), intent(out) :: buf
    call C_F_POINTER(ptr,branch)
    call branch%branch%branch_get_xc_buf(buf)
    return
  end subroutine p_branch_get_xc_buf
#if 1
!
! Return size of matrix block on the diagonal contributed by branch
! @param ptr C pointer to GridPACK branch object
! @param isize,jsize number of rows and columns in block
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function p_branch_matrix_diag_size(ptr, isize, jsize) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer(C_INT), intent(out) :: isize, jsize
    logical f_ret
    integer f_isize, f_jsize
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_matrix_diag_size(f_isize, f_jsize)
    isize = f_isize
    jsize = f_isize
    p_branch_matrix_diag_size = f_ret
    return
  end function p_branch_matrix_diag_size
!
! Return the values for a diagonal matrix block. The values are returned in
! row-major order.
! @param ptr C pointer to GridPACK branch object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function p_branch_matrix_diag_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    logical f_ret
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_matrix_diag_values(values)
    p_branch_matrix_diag_values = f_ret
    return
  end function p_branch_matrix_diag_values
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the forward direction
! @param ptr C pointer to GridPACK branch object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function p_branch_matrix_forward_size(ptr, isize, jsize) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer(C_INT), intent(out) :: isize, jsize
    logical f_ret
    integer f_isize, f_jsize
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_matrix_forward_size(f_isize, f_jsize)
    isize = f_isize
    jsize = f_jsize
    p_branch_matrix_forward_size = f_ret
    return
  end function p_branch_matrix_forward_size
!
! Return size of off-diagonal matrix block contributed by component. The values
! are in the reverse direction
! @param ptr C pointer to GridPACK branch object
! @param isize,jsize number of rows and columns of matrix block
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function p_branch_matrix_reverse_size(ptr, isize, jsize) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer(C_INT), intent(out) :: isize, jsize
    logical f_ret
    integer f_isize, f_jsize
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_matrix_reverse_size(f_isize, f_jsize)
    isize = f_isize
    jsize = f_jsize
    p_branch_matrix_reverse_size = f_ret
    return
  end function p_branch_matrix_reverse_size
!
! Return the values of off-diagonal matrix block. The values are for the forward
! direction and are returned in row-major order.
! @param ptr C pointer to GridPACK branch object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function p_branch_matrix_forward_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    logical f_ret
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_matrix_forward_values(values)
    p_branch_matrix_forward_values = f_ret
    return
  end function p_branch_matrix_forward_values
!
! Return the values of off-diagonal matrix block. The values are for the reverse
! direction and are returned in row-major order.
! @param ptr C pointer to GridPACK branch object
! @param values array that contains matrix values
! @return false if network component does not contribute matrix element
!
  logical(C_BOOL) function p_branch_matrix_reverse_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    logical f_ret
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_matrix_reverse_values(values)
    p_branch_matrix_reverse_values = f_ret
    return
  end function p_branch_matrix_reverse_values
!
! Return the size of vector block contributed by component
! @param ptr C pointer to GridPACK branch object
! @param isize number of vector elements
! @return false if network component does not contribute vector element
!
  logical(C_BOOL) function p_branch_vector_size(ptr, isize) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer(C_INT), intent(out) :: isize
    integer f_isize
    logical f_ret
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_vector_size(f_isize)
    isize = f_isize
    p_branch_vector_size = f_ret
    return
  end function p_branch_vector_size
!
! Return the values of the vector block
! @param ptr C pointer to GridPACK branch object
! @param values array that contains vector values
! @return false if network component does not contribute vector element
! 
  logical(C_BOOL) function p_branch_vector_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    logical f_ret
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_vector_values(values)
    p_branch_vector_values = f_ret
  end function p_branch_vector_values
!
! Set values in network component based on values in a vector or matrix
! @param ptr C pointer to GridPACK branch object
! @param values array that contains vector values
! 
  subroutine p_branch_set_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    complex(C_DOUBLE_COMPLEX), intent(in) :: values(*)
    call C_F_POINTER(ptr,branch)
    call branch%branch%branch_set_values(values)
  end subroutine p_branch_set_values
!
! Return the number of rows in matrix from component
! @param ptr C pointer to GridPACK branch object
! @return number of rows from component
!
  integer(C_INT) function p_branch_matrix_num_rows(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer f_ret
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_matrix_num_rows()
    p_branch_matrix_num_rows = f_ret
  end function p_branch_matrix_num_rows
!
! Return the number of columns in matrix from component
! @param ptr C pointer to GridPACK branch object
! @return number of columns from component
!
  integer(C_INT) function p_branch_matrix_num_cols(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer f_ret
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_matrix_num_cols()
    p_branch_matrix_num_cols = f_ret
  end function p_branch_matrix_num_cols
!
! Set row indices corresponding to the rows contributed by this component
! @param ptr C pointer to GridPACK branch object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @param idx matrix index of row irow
!
  subroutine p_branch_matrix_set_row_index(ptr, irow, r_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer(C_INT), value, intent(in) :: irow, r_idx
    integer f_row, f_idx
    call C_F_POINTER(ptr,branch)
    f_row = irow
    f_idx = r_idx
    call branch%branch%branch_matrix_set_row_index(f_row, f_idx)
    return
  end subroutine p_branch_matrix_set_row_index
!
! Set column indices corresponding to the columns contributed by this component
! @param ptr C pointer to GridPACK branch object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @param idx matrix index of column icol
!
  subroutine p_branch_matrix_set_col_index(ptr, icol, c_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer(C_INT), value, intent(in) :: icol, c_idx
    integer f_col, f_idx
    call C_F_POINTER(ptr,branch)
    f_col = icol
    f_idx = c_idx
    call branch%branch%branch_matrix_set_col_index(f_col, f_idx)
    return
  end subroutine p_branch_matrix_set_col_index
!
! Get the row index for the rows contributed by this component
! @param ptr C pointer to GridPACK branch object
! @param irow index of row contributed by this component (e.g. if component
! contributes 3 rows then irow is between 0 and 2)
! @return matrix index of row irow
!
  integer(C_INT) function p_branch_matrix_get_row_index(ptr, irow) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer(C_INT), value, intent(in) :: irow
    integer f_row, f_ret
    call C_F_POINTER(ptr,branch)
    f_row = irow
    f_ret = branch%branch%branch_matrix_get_row_index(f_row)
    p_branch_matrix_get_row_index = f_ret
    return
  end function p_branch_matrix_get_row_index
!
! Get the col index for the columns contributed by this component
! @param ptr C pointer to GridPACK branch object
! @param icol index of column contributed by this component (e.g. if component
! contributes 3 columns then icol is between 0 and 2)
! @return matrix index of column icol
!
  integer(C_INT) function p_branch_matrix_get_col_index(ptr, icol) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer(C_INT), value, intent(in) :: icol
    integer f_col, f_ret
    call C_F_POINTER(ptr,branch)
    f_col = icol
    f_ret = branch%branch%branch_matrix_get_col_index(f_col)
    p_branch_matrix_get_col_index = f_ret
  end function p_branch_matrix_get_col_index
!
! Return the number of matrix values contributed by this component
! @param ptr C pointer to GridPACK branch object
! @return number of matrix values
!
  integer(C_INT) function p_branch_matrix_num_values(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer f_ret
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_matrix_num_values()
    p_branch_matrix_num_values = f_ret
  end function p_branch_matrix_num_values
!
! Get a list of matrix values contributed by this components along with their
! matrix indices
! @param ptr C pointer to GridPACK branch object
! @param values list of matrix element values
! @param rows row indices for the matrix elements
! @param cols column indices for the matrix elements
!
  subroutine p_branch_matrix_get_values(ptr, values, rows, cols) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    integer(C_INT), intent(out) :: rows(*), cols(*)
    integer nval, i
    integer, allocatable :: f_rows(:), f_cols(:)
    call C_F_POINTER(ptr,branch)
    nval = branch%branch%branch_matrix_num_values()
    allocate(f_rows(nval))
    allocate(f_cols(nval))
    call branch%branch%branch_matrix_get_values(values,f_rows,f_cols)
    do i = 1, nval
      rows(i) = f_rows(i)
      cols(i) = f_cols(i)
    end do
    deallocate(f_rows)
    deallocate(f_cols)
    return
  end subroutine p_branch_matrix_get_values
!
! Return the number of elements in vector coming from component
! @param ptr C pointer to GridPACK branch object
! @return number of elements contributed from component
!
  integer(C_INT) function p_branch_vector_num_elements(ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer f_ret
    call C_F_POINTER(ptr,branch)
    f_ret = branch%branch%branch_vector_num_elements()
    p_branch_vector_num_elements = f_ret
    return
  end function p_branch_vector_num_elements
!
! Set indices corresponding to the elements contributed by this component
! @param ptr C pointer to GridPACK branch object
! @param ielem index of element contributed by this component (e.g. if component
! contributes 3 elements then ielem is between 0 and 2)
! @param idx vector index of element ielem
!
  subroutine p_branch_vector_set_element_index(ptr, ielem, e_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer(C_INT), value, intent(in) :: ielem, e_idx
    integer f_elem, f_idx
    call C_F_POINTER(ptr,branch)
    f_elem = ielem
    f_idx = e_idx
    call branch%branch%branch_vector_set_element_index(f_elem, f_idx)
    return
  end subroutine p_branch_vector_set_element_index
!
! Get a list of vector indices from the component
! @param ptr C pointer to GridPACK branch object
! @param idx list of indices the component maps onto
!
  subroutine p_branch_vector_get_element_indices(ptr, e_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer(C_INT), intent(out) :: e_idx(*)
    integer, allocatable :: f_idx(:)
    integer nval, i
    call C_F_POINTER(ptr,branch)
    nval = branch%branch%branch_vector_num_elements()
    allocate(f_idx(nval))
    call branch%branch%branch_vector_get_element_indices(f_idx)
    do i = 1, nval
      e_idx(i) = f_idx(i)
    end do
    deallocate(f_idx)
    return
  end subroutine p_branch_vector_get_element_indices
!
! Get a list of vector values contributed by this component and their indices
! @param ptr C pointer to GridPACK branch object
! @param values list of vector element values
! @param idx indices of the vector elements
!
  subroutine p_branch_vector_get_element_values(ptr, values, e_idx) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    integer(C_INT), intent(out) :: e_idx(*)
    integer, allocatable :: f_idx(:)
    integer nval, i
    call C_F_POINTER(ptr,branch)
    nval = branch%branch%branch_vector_num_elements()
    allocate(f_idx(nval))
    call branch%branch%branch_vector_get_element_values(values, f_idx)
    do i = 1, nval
      e_idx(i) = f_idx(i)
    end do
    deallocate(f_idx)
  end subroutine p_branch_vector_get_element_values
!
! Transfer vector values to component
! @param ptr C pointer to GridPACK branch object
! @param values list of vector element values
!
  subroutine p_branch_vector_set_element_values(ptr, values) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    complex(C_DOUBLE_COMPLEX), intent(out) :: values(*)
    call C_F_POINTER(ptr,branch)
    call branch%branch%branch_vector_set_element_values(values)
    return
  end subroutine p_branch_vector_set_element_values
!
! Load data from DataCollection object into corresponding component.
! @param ptr C pointer to GridPACK branch object
! @param data DataCollection object associated with component
!
  subroutine p_branch_load(ptr, data_ptr) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(C_PTR), value, intent(in) :: data_ptr
    type(application_branch_wrapper), pointer :: branch
    type(data_collection), pointer :: data
    call C_F_POINTER(ptr,branch)
    call C_F_POINTER(data_ptr,data)
    call branch%branch%branch_load(data)
    return
  end subroutine p_branch_load
!
! Set an internal variable that can be used to control the behavior of the
! component. This function doesn't need to be implemented, but if needed it can
! be used to change the behavior of the network in different phases of the
! calculation. For example, if a different matrix needs to be generated at
! different times, the mode can be change to get different values from the
! matrix-vector interface functions.
! @param ptr C pointer to GridPACK branch object
! @param mode integer indicating which mode should be used
!
  subroutine p_branch_set_mode(ptr, mode) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    integer(C_INT), value, intent(in) :: mode
    integer f_mode
    call C_F_POINTER(ptr,branch)
    f_mode = mode
    call branch%branch%branch_set_mode(f_mode)
    return
  end subroutine p_branch_set_mode
!
! Copy a string for output into buffer. The behavior of this method can be
! altered by inputting different values for the signal string
! @param ptr C pointer to GridPACK branch object
! @param string buffer containing string to be written to output
! @param bufsize size of string buffer in bytes
! @param signal string to control behavior of routine (e.g. what
! properties to write
! @return true if component is writing a contribution, false otherwise
!
  logical(C_BOOL) function p_branch_serial_write(ptr, string, bufsize, &
      signal, signal_len, write_len) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), value, intent(in) :: ptr
    type(application_branch_wrapper), pointer :: branch
    character(C_CHAR), intent(out) :: string(*)
    integer(C_INT), value, intent(in) :: bufsize, signal_len
    character(C_CHAR), intent(in) :: signal(*)
    integer(C_INT), intent(out) :: write_len 
    character(len=bufsize) fstring
    character(len=signal_len) fsignal
    logical f_ret
    integer f_bufsize, slen, i, fslen
    call C_F_POINTER(ptr,branch)
    f_bufsize = bufsize
    fslen = signal_len
    do i = 1, fslen
      fsignal(i:i) = signal(i)
    end do
    f_ret = branch%branch%branch_serial_write(fstring,f_bufsize,fsignal)
    slen = len(trim(fstring))
    slen = min(bufsize-1,slen)
    i = 1
    do while (i.le.slen.and.f_ret)
      string(i) = fstring(i:i)
      i = i + 1
    end do
    string(slen+1) = C_NULL_CHAR
    write_len = slen
    p_branch_serial_write = f_ret
    return
  end function p_branch_serial_write
#endif
end module application_components
