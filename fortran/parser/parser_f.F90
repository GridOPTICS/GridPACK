!
!  Fortran parser functions
!
module gridpack_mapper
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
      type(C_PTR), intent(in) :: network
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
      type(C_PTR), intent(in) :: parser
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
    character, intent(in) :: string(*)
    call pti23_parser_parse(p_parser%p_parser, string)
    return
  end subroutine parse
end module gridpack_mapper
