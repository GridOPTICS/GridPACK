! ----------------------------------------------------------------
! file: mapper_factory_f.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created September 15, 2014 by Bruce Palmer
! ----------------------------------------------------------------
!
!  Fortran factory functions
!
module application_factory
  use, intrinsic :: iso_c_binding
  use gridpack_network
  use gridpack_factory
  implicit none
!
!  Define factory type
!
  private

  type, extends(factory), public :: app_factory
    class(network), pointer :: p_network_int
!
!  Add factory parameters here
!
    contains
!
!  Add application-specific methods here
!
!
!  Do not modify or remove the create and destroy methods
!
    procedure::create
    procedure::destroy
  end type
!
  interface
!
! Create a new factory
! @param factory new GridPACK factory object
! @param network pointer to GridPACK network object
!
    subroutine factory_create(p_factory, p_network) bind(c)
      use, intrinsic :: iso_c_binding
      use gridpack_network
      implicit none
      type(C_PTR), intent(inout) :: p_factory
      type(C_PTR), value, intent(in) :: p_network
    end subroutine factory_create
!
! Clean up old factory
! @param factory old GridPACK factory object
!
    subroutine factory_destroy(p_factory) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: p_factory
    end subroutine factory_destroy
  end interface
!
  contains
!
!  Add application-specific methods here
!
!
!  DO NOT REMOVE THIS INCLUDE FILE. IT CONTAINS FUNCTIONS THAT MUST
!  BE INCLUDED IN THIS FILE BUT SHOULD NOT BE MODIFIED BY THE
!  APPLICATION DEVELOPER
!
  include 'factory_inc.F90'
end module application_factory
