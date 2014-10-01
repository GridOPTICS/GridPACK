! ----------------------------------------------------------------
! file: pf_factory.f90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
!
!  Fortran factory functions
!
module application_factory
  use, intrinsic :: iso_c_binding
  use gridpack_network
  use gridpack_factory
  use application_components
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
    procedure::set_y_bus
    procedure::set_s_bus
    procedure::set_pq
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
! Initialize parameters for the admittance (Y-bus) matrix
! @param factory Fortran factory object
!
  subroutine set_y_bus(factory)
    class(app_factory), intent(in) :: factory
    class(application_bus), pointer :: bus
    class(application_branch), pointer :: branch
    class(network), pointer :: grid
    integer nbus, nbranch, i
    grid => factory%p_network_int
    nbus = grid%num_buses()
    nbranch = grid%num_branches()
    do i = 1, nbus
      bus => bus_cast(grid%get_bus(i))
      call bus%bus_set_y_matrix()
    end do
    do i = 1, nbranch
      branch => branch_cast(grid%get_branch(i))
      call branch%branch_set_y_matrix()
    end do
    return
  end subroutine set_y_bus
!
! Construct S-bus vector
! @param factory Fortran factory object
!
  subroutine set_s_bus(factory)
    class(app_factory), intent(in) :: factory
    class(application_bus), pointer :: bus
    class(network), pointer :: grid
    integer nbus, i
    grid => factory%p_network_int
    nbus = grid%num_buses()
    do i = 1, nbus
      bus => bus_cast(grid%get_bus(i))
      call bus%bus_set_sbus()
    end do
    return
  end subroutine set_s_bus
!
! Construct the PQ vector
! @param factory Fortran factory object
!
  subroutine set_pq(factory)
    class(app_factory), intent(in) :: factory
    class(application_bus), pointer :: bus
    class(network), pointer :: grid
    integer nbus, i
    double complex v(2)
    logical ok
    grid => factory%p_network_int
    nbus = grid%num_buses()
    do i = 1, nbus
      bus => bus_cast(grid%get_bus(i))
      ok = bus%bus_vector_values(v)
    end do
    return
  end subroutine set_pq
!
!  DO NOT REMOVE THIS INCLUDE FILE. IT CONTAINS FUNCTIONS THAT MUST
!  BE INCLUDED IN THIS FILE BUT SHOULD NOT BE MODIFIED BY THE
!  APPLICATION DEVELOPER
!
  include 'factory_inc.F90'
end module application_factory
