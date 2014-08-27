!
!  Fortran data collection component
!
module gridpack_data_collection
  use, intrinsic :: iso_c_binding
  implicit none
!
!  Define data collection type
!
  private
  type, public :: data_collection
    type (C_PTR) :: p_data
    integer :: p_idx
    contains
    procedure::get_int_value
  end type
!
  type, public :: data_wrapper
    type(data_collection), pointer :: data
  end type
!
!  Interface declaration to C calls
!
  interface
!
! Retrieve current value of existing data element in
! data_collection object
! @param data pointer to GridPACK data collection object
! @param name name of data element
! @param value current value of data element
! @return false if no element of the correct name and type exists in
! data_collection object
!
    logical(C_BOOL) function data_collection_get_int_value(data,name,value) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: data
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT), intent(out) :: value
    end function data_collection_get_int_value
  end interface
  contains
!
! Retrieve current value of existing data element in
! data_collection object
! @param p_data pointer to GridPACK data collection object
! @param name name of data element
! @param value current value of data element
! @return false if no element of the correct name and type exists in
! data_collection object
!
  logical function get_int_value(p_data,name,value)
    use, intrinsic :: iso_c_binding
    implicit none
    class(data_collection), intent(in) :: p_data
    character, intent(in) :: name(*)
    integer, intent(out) :: value
    logical c_ret
    integer(C_INT) c_value
    c_ret = data_collection_get_int_value(p_data%p_data,name,c_value)
    value = c_value
    get_int_value = c_ret
    return
  end function get_int_value
end module gridpack_data_collection
