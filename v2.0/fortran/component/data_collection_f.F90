! ----------------------------------------------------------------
! file: data_collection_f.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created September 30, 2014 by Bruce Palmer
! Last Change: 2014-09-30 14:07:10 d3g293
! ----------------------------------------------------------------
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
    procedure::get_logical_value
    procedure::get_string_value
    procedure::get_real_value
    procedure::get_double_value
    procedure::get_int_indexed_value
    procedure::get_logical_indexed_value
    procedure::get_string_indexed_value
    procedure::get_real_indexed_value
    procedure::get_double_indexed_value
    procedure::bind
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
    logical(C_BOOL) function data_collection_get_logical_value(data,name,value) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: data
      character(C_CHAR), intent(in) :: name(*)
      logical(C_BOOL), intent(out) :: value
    end function data_collection_get_logical_value
    logical(C_BOOL) function data_collection_get_string_value(data,name,value) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: data
      character(C_CHAR), intent(in) :: name(*)
      character(C_CHAR), intent(out) :: value(*)
    end function data_collection_get_string_value
    logical(C_BOOL) function data_collection_get_real_value(data,name,value) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: data
      character(C_CHAR), intent(in) :: name(*)
      real(C_FLOAT), intent(out) :: value
    end function data_collection_get_real_value
    logical(C_BOOL) function data_collection_get_double_value(data,name,value) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: data
      character(C_CHAR), intent(in) :: name(*)
      real(C_DOUBLE), intent(out) :: value
    end function data_collection_get_double_value
!
! Retrieve current value of existing data element in
! data_collection object. Assume that the item appears in data_collection with
! tag "name:idx"
! @param data pointer to GridPACK data collection object
! @param name name of data element
! @param value current value of data element
! @param idx index of value
! @return false if no element of the correct name and type exists in
! data_collection object
!
    logical(C_BOOL) function data_collection_get_int_indexed_value(data, &
        name,value,idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: data
      character(C_CHAR), intent(in) :: name(*)
      integer(C_INT), intent(out) :: value
      integer(C_INT), value, intent(in) :: idx
    end function data_collection_get_int_indexed_value
    logical(C_BOOL) function data_collection_get_logical_indexed_value(data, &
       name,value,idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: data
      character(C_CHAR), intent(in) :: name(*)
      logical(C_BOOL), intent(out) :: value
      integer(C_INT), value, intent(in) :: idx
    end function data_collection_get_logical_indexed_value
    logical(C_BOOL) function data_collection_get_string_indexed_value(data, &
        name,value,idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: data
      character(C_CHAR), intent(in) :: name(*)
      character(C_CHAR), intent(inout) :: value(*)
      integer(C_INT), value, intent(in) :: idx
    end function data_collection_get_string_indexed_value
    logical(C_BOOL) function data_collection_get_real_indexed_value(data, &
        name,value,idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: data
      character(C_CHAR), intent(in) :: name(*)
      real(C_FLOAT), intent(out) :: value
      integer(C_INT), value, intent(in) :: idx
    end function data_collection_get_real_indexed_value
    logical(C_BOOL) function data_collection_get_double_indexed_value(data, &
        name,value,idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: data
      character(C_CHAR), intent(in) :: name(*)
      real(C_DOUBLE), intent(out) :: value
      integer(C_INT), value, intent(in) :: idx
    end function data_collection_get_double_indexed_value
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
    character(len=*), intent(in) :: name
    integer, intent(out) :: value
    logical c_ret
    integer(C_INT) c_value
    character(C_CHAR) :: cstring(512)
    call fstring_convert(name,cstring)
    c_ret = data_collection_get_int_value(p_data%p_data,cstring,c_value)
    value = c_value
    get_int_value = c_ret
    return
  end function get_int_value
  logical function get_logical_value(p_data,name,value)
    use, intrinsic :: iso_c_binding
    implicit none
    class(data_collection), intent(in) :: p_data
    character(len=*), intent(in) :: name
    logical, intent(out) :: value
    logical c_ret
    logical(C_BOOL) c_value
    character(C_CHAR) :: cstring(512)
    call fstring_convert(name,cstring)
    c_ret = data_collection_get_logical_value(p_data%p_data,cstring,c_value)
    value = c_value
    get_logical_value = c_ret
    return
  end function get_logical_value
  logical function get_string_value(p_data,name,value)
    use, intrinsic :: iso_c_binding
    implicit none
    class(data_collection), intent(in) :: p_data
    character(len=*), intent(in) :: name
    character(C_CHAR), intent(inout) :: value(*)
    logical c_ret
    character(C_CHAR) :: cstring(512)
    call fstring_convert(name,cstring)
    c_ret = data_collection_get_string_value(p_data%p_data,cstring,value)
    get_string_value = c_ret
    return
  end function get_string_value
  logical function get_real_value(p_data,name,value)
    use, intrinsic :: iso_c_binding
    implicit none
    class(data_collection), intent(in) :: p_data
    character(len=*), intent(in) :: name
    real, intent(out) :: value
    logical c_ret
    real(C_FLOAT) c_value
    character(C_CHAR) :: cstring(512)
    call fstring_convert(name,cstring)
    c_ret = data_collection_get_real_value(p_data%p_data,cstring,value)
    value = c_value
    get_real_value = c_ret
    return
  end function get_real_value
  logical function get_double_value(p_data,name,value)
    use, intrinsic :: iso_c_binding
    implicit none
    class(data_collection), intent(in) :: p_data
    character(len=*), intent(in) :: name
    double precision, intent(out) :: value
    logical c_ret
    real(C_DOUBLE) c_value
    character(C_CHAR) :: cstring(512)
    call fstring_convert(name,cstring)
    c_ret = data_collection_get_double_value(p_data%p_data,cstring,c_value)
    value = c_value
    get_double_value = c_ret
    return
  end function get_double_value
!
! Retrieve current value of existing data element in
! data_collection object. Assume that the item appears in data_collection with
! tag "name:idx"
! @param p_data pointer to GridPACK data collection object
! @param name name of data element
! @param value current value of data element
! @param idx index of value
! @return false if no element of the correct name and type exists in
! data_collection object
!
  logical function get_int_indexed_value(p_data,name,value,idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(data_collection), intent(in) :: p_data
    character(len=*), intent(in) :: name
    integer, intent(out) :: value
    integer, intent(in) :: idx
    logical c_ret
    integer(C_INT) c_idx
    integer(C_INT) c_value
    character(C_CHAR) :: cstring(512)
    call fstring_convert(name,cstring)
    c_idx = idx
    c_ret = data_collection_get_int_indexed_value(p_data%p_data,cstring,c_value,c_idx)
    value = c_value
    get_int_indexed_value = c_ret
    return
  end function get_int_indexed_value
  logical function get_logical_indexed_value(p_data,name,value,idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(data_collection), intent(in) :: p_data
    character(len=*), intent(in) :: name
    logical, intent(out) :: value
    integer, intent(in) :: idx
    logical c_ret
    integer(C_INT) c_idx
    logical(C_BOOL) c_value
    character(C_CHAR) :: cstring(512)
    call fstring_convert(name,cstring)
    c_idx = idx
    c_ret = data_collection_get_logical_indexed_value(p_data%p_data,cstring,c_value,c_idx)
    value = c_value
    get_logical_indexed_value = c_ret
    return
  end function get_logical_indexed_value
  logical function get_string_indexed_value(p_data,name,value,idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(data_collection), intent(in) :: p_data
    character(len=*), intent(in) :: name
    character, intent(inout) :: value(*)
    integer, intent(in) :: idx
    logical c_ret
    integer(C_INT) c_idx
    character(C_CHAR) :: cstring(512)
    call fstring_convert(name,cstring)
    c_idx = idx
    c_ret = data_collection_get_string_indexed_value(p_data%p_data,cstring,value,c_idx)
    get_string_indexed_value = c_ret
    return
  end function get_string_indexed_value
  logical function get_real_indexed_value(p_data,name,value,idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(data_collection), intent(in) :: p_data
    character(len=*), intent(in) :: name
    real, intent(out) :: value
    integer, intent(in) :: idx
    logical c_ret
    integer(C_INT) c_idx
    real(C_FLOAT) c_value
    character(C_CHAR) :: cstring(512)
    call fstring_convert(name,cstring)
    c_idx = idx
    c_ret = data_collection_get_real_indexed_value(p_data%p_data,cstring,c_value,c_idx)
    value = c_value
    get_real_indexed_value = c_ret
    return
  end function get_real_indexed_value
  logical function get_double_indexed_value(p_data,name,value,idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(data_collection), intent(in) :: p_data
    character(len=*), intent(in) :: name
    double precision, intent(out) :: value
    integer, intent(in) :: idx
    logical c_ret
    integer(C_INT) c_idx
    real(C_DOUBLE) c_value
    character(C_CHAR) :: cstring(512)
    call fstring_convert(name,cstring)
    c_idx = idx
    c_ret = data_collection_get_double_indexed_value(p_data%p_data,cstring,c_value,c_idx)
    value = c_value
    get_double_indexed_value = c_ret
    return
  end function get_double_indexed_value
!
! This subroutine binds a pointer to the C DataCollection object to the
! corresponding Fortran object
! @param ptr pointer to C DataCollection object
! @return Fortran pointer to data collection object
!
  subroutine bind(p_data, ptr)
    use, intrinsic :: iso_c_binding
    implicit none
    class(data_collection), intent(inout) :: p_data
    type(C_PTR), value, intent(in) :: ptr
    p_data%p_data = ptr
    return
  end subroutine bind
!
! Utility routine to add C_NULL_CHAR to end of fortran strings
! @param fstring string coming from Fortran calling program
! @param cstring string with appended null character
!
   subroutine fstring_convert(fstring,cstring)
     use, intrinsic :: iso_c_binding
     implicit none
     character(len=*), intent(in) :: fstring
     character(C_CHAR), intent(inout):: cstring(*)
     integer flen, i
     flen = len(trim(fstring))
     do i = 1, flen
       cstring(i) = fstring(i:i)
     end do
     cstring(flen+1) = C_NULL_CHAR
   end subroutine fstring_convert
end module gridpack_data_collection
