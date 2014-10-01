! ----------------------------------------------------------------
! file: factory_inc.f90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created October 1, 2014 by Bruce Palmer
! ----------------------------------------------------------------
!
!  Fortran application factory include file. These are functions
!  that need to included in the application factory file but
!  should not be modified by the developer
!
!  DO NOT EDIT ANYTHING IN THIS FILE. THESE FUNCTIONS MUST BE INCLUDED IN
!  THIS FILE BUT SHOULD NOT BE MODIFIED BY THE APPLICATION DEVELOPER
!
! Create a new factory
! @param p_factory new GridPACK factory object
! @param p_network pointer to GridPACK network object
!
  subroutine create(p_factory, p_network)
    use, intrinsic :: iso_c_binding
    use gridpack_network
    implicit none
    class(app_factory), intent(inout) :: p_factory
    class(network), target, intent(in) :: p_network
    p_factory%p_network_int => p_network
    call factory_create(p_factory%p_factory,p_network%p_network)
    return
  end subroutine create
!
! Clean up old factory
! @param p_factory old GridPACK factory object
!
  subroutine destroy(p_factory)
    use, intrinsic :: iso_c_binding
    implicit none
    class(app_factory), intent(inout) :: p_factory
    call factory_destroy(p_factory%p_factory)
    return
  end subroutine destroy
