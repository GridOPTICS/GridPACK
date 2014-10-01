! ----------------------------------------------------------------
! file: parser_f.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created September 17, 2014 by Bruce Palmer
! ----------------------------------------------------------------
!
!  Fortran parser functions
!
module gridpack_parser
  use, intrinsic :: iso_c_binding
  use gridpack_network
  implicit none
!
! Define parser types
!
  private
!
  type, public :: pti23_parser
    type(C_PTR) :: p_parser
    contains
    procedure::create
    procedure::destroy
    procedure::parse
  end type
!
  interface
!
! Create a PTI23 parser
! @param parser pointer to Fortran pti23 parser object
! @param network pointer to Fortran network object
!
    subroutine pti23_parser_create(parser, network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: parser
      type(C_PTR), value, intent(in) :: network
    end subroutine pti23_parser_create
!
! Destroy a PTI23 parser
! @param parser pointer to Fortran pti23 parser object
!
    subroutine pti23_parser_destroy(parser) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: parser
    end subroutine pti23_parser_destroy
!
! Parser network configuration file and create network
! @param parser pointer to Fortran pti23 parser object
! @param name file name of network file
!
    subroutine pti23_parser_parse(parser, string) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: parser
      character(C_CHAR), intent(in) :: string(*)
    end subroutine pti23_parser_parse
  end interface
  contains
!
! Create a PTI23 parser
! @param p_parser pointer to Fortran pti23 parser object
! @param p_network pointer to Fortran network object
!
  subroutine create(p_parser, p_network)
    use, intrinsic :: iso_c_binding
    use gridpack_network
    implicit none
    class(pti23_parser), intent(inout) :: p_parser
    class(network), intent(in) :: p_network
    call pti23_parser_create(p_parser%p_parser,p_network%p_network)
    return
  end subroutine create
!
! Destroy a PTI23 parser
! @param p_parser pointer to Fortran pti23 parser object
!
  subroutine destroy(p_parser)
    use, intrinsic :: iso_c_binding
    implicit none
    class(pti23_parser), intent(inout) :: p_parser
    call pti23_parser_destroy(p_parser%p_parser)
    return
  end subroutine destroy
!
! Parse network configuration file and create network
! @param p_parser pointer to Fortran pti23 parser object
! @param string name of network file
!
  subroutine parse(p_parser, string)
    use, intrinsic :: iso_c_binding
    implicit none
    class(pti23_parser), intent(in) :: p_parser
    character(len=*), intent(in) :: string
    character(C_CHAR) str_dummy(512)
    integer i,slen
    slen = len(trim(string))
    if (slen.gt.511) then
      write(6,'(a)') 'File name too long in pti23_parser%parse'
    endif
    do i=1, slen
      str_dummy(i) = string(i:i)
    end do
    str_dummy(slen+1) = C_NULL_CHAR
    call pti23_parser_parse(p_parser%p_parser, str_dummy)
    return
  end subroutine parse
end module gridpack_parser
