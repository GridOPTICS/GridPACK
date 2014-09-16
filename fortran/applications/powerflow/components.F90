#define LARGE_MATRIX
module application_components
  use gridpack_component
  use gridpack_data_collection
  use iso_c_binding
  implicit none
  integer JACOBIAN, YBUS, RHS, STATE, SCAL
  parameter(JACOBIAN=1, YBUS=2, RHS=3, STATE=4, SCAL=5)
!
!  Define bus component
!
  private
  type, bind(c), public :: bus_xc_data
    real(C_DOUBLE) :: v_ang
    real(C_DOUBLE) :: v_mag
  end type
!
  type, public :: application_bus
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
    double precision p_pl, p_ql, p_sbase, p_pinj, p_qinj
    logical p_ispv, p_saveispv
    type (bus_xc_data) xc_buf
    contains
!
!  Y-matrix methods
!
    procedure :: bus_matrix_diag_size
    procedure :: bus_matrix_diag_values
    procedure :: bus_set_y_matrix !set_ybus
    procedure :: bus_get_y_matrix !get_ybus
    procedure :: bus_is_isolated
    procedure :: bus_set_isolated
!
!  Powerflow methods
!
    procedure :: bus_matrix_diag_size
    procedure :: bus_matrix_diag_values
    procedure :: bus_vector_size
    procedure :: bus_vector_values
    procedure :: bus_set_values
    procedure :: bus_load
    procedure :: bus_set_mode
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
    procedure :: bus_serial_write
    procedure :: bus_chk_q_lim
  end type
!
! Return the size of matrix block on the diagonal contributed by bus
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
! Return the values of the block for the bus. The values are returned in
! row-major order
! @param bus GridPACK bus object
! @param values array of matrix block values
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
          bus_matrix_diag_values = .falsee.
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
  end function bus_matrix_diag_values
!
! Return the size of the block that this bus contributes to vector
! @param bus GridPACK bus object
! @param size size of vector block
! @return false if network component does not contribute vector element
!
  bool function bus_vector_size(bus, size)
    implicit none
    class(application_bus), intent(in) :: bus
    integer, intent(out) :: size
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
! Return the values of the vector block contributed by this bus
! @param bus GridPACK bus object
! @param values array of vector values
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
          nbranch = bus_get_num_neighbors
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
end module application_components
