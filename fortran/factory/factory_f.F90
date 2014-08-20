!
!  Fortran factory functions
!
module gridpack_factory
  use, intrinsic :: iso_c_binding
  implicit none
!
!  Define factory type
!
  private

  type, public :: factory
    type(C_PTR) :: p_factory
    contains
    procedure::create
    procedure::destroy
    procedure::set_components
    procedure::load
    procedure::set_exchange
    procedure::set_mode
    procedure::check_true
  end type
  interface
!
! Create a new factory
! @param factory new GridPACK factory object
!
    subroutine factory_create(factory, network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: factory
      type(C_PTR), intent(in) :: network
    end subroutine factory_create
!
! Clean up old factory
! @param factory old GridPACK factory object
!
    subroutine factory_destroy(factory) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: factory
    end subroutine factory_destroy
!
! Set pointers in each bus and branch component so that it points to
! connected buses and branches.
! @param factory GridPACK factory object
!
    subroutine factory_set_components(factory) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: factory
    end subroutine factory_set_components
!
! Generic method that invokes the "load" method on all branches and buses
! to move data from the data collection objects on the network into the
! corresponding buses and branches
! @param factory GridPACK factory object
!
    subroutine factory_load(factory) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: factory
    end subroutine factory_load
!
! Set up the exchange buffers so that they work correctly. This should only
! be called after the network topology has been specified
! @param factory GridPACK factory object
!
    subroutine factory_set_exchange(factory) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: factory
    end subroutine factory_set_exchange
!
! Set the mode for all network component objects in the network.
! @param factory GridPACK factory object
! @param mode integer representing desired mode
!
    subroutine factory_set_mode(factory, mode) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: factory
      integer(C_INT), value, intent(in) :: mode
    end subroutine factory_set_mode
!
! A convenience function that checks to see if something is true on all
! processors
! @param factory GridPACK factory object
! @param flag boolean flag on each processor
! @return true if flag is true on all processors, false otherwise
!
    logical(C_BOOL) function factory_check_true(factory, flag) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: factory
      logical(C_BOOL), value, intent(in) :: flag
    end function factory_check_true
  end interface
  contains
!
! Create a new factory
! @param p_factory new GridPACK factory object
!
  subroutine create(p_factory, p_network)
    use, intrinsic :: iso_c_binding
    use gridpack_network
    implicit none
    class(factory), intent(inout) :: p_factory
    class(network), intent(in) :: p_network
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
    class(factory), intent(inout) :: p_factory
    call factory_destroy(p_factory%p_factory)
    return
  end subroutine destroy
!
! Set pointers in each bus and branch component so that it points to
! connected buses and branches.
! @param p_factory GridPACK factory object
!
  subroutine set_components(p_factory)
    use, intrinsic :: iso_c_binding
    implicit none
    class(factory), value, intent(in) :: p_factory
    call factory_set_components(p_factory%p_factory)
    return
  end subroutine set_components
!
! Generic method that invokes the "load" method on all branches and buses
! to move data from the data collection objects on the network into the
! corresponding buses and branches
! @param p_factory GridPACK factory object
!
  subroutine load(p_factory)
    use, intrinsic :: iso_c_binding
    implicit none
    class(factory), value, intent(in) :: p_factory
    call factory_load(p_factory%p_factory)
    return
  end subroutine load
!
! Set up the exchange buffers so that they work correctly. This should only
! be called after the network topology has been specified
! @param p_factory GridPACK factory object
!
  subroutine set_exchange(p_factory)
    use, intrinsic :: iso_c_binding
    implicit none
    class(factory), value, intent(in) :: p_factory
    call factory_set_exchange(p_factory%p_factory)
    return
  end subroutine set_exchange
!
! Set the mode for all network component objects in the network.
! @param p_factory GridPACK factory object
! @param mode integer representing desired mode
!
  subroutine set_mode(p_factory, mode)
    use, intrinsic :: iso_c_binding
    implicit none
    class(factory), value, intent(in) :: p_factory
    integer, value, intent(in) :: mode
    integer(C_INT) c_mode
    c_mode = mode
    call factory_set_mode(p_factory%p_factory, c_mode)
    return
  end subroutine set_mode
!
! A convenience function that checks to see if something is true on all
! processors
! @param p_factory GridPACK factory object
! @param flag boolean flag on each processor
! @return true if flag is true on all processors, false otherwise
!
  logical(C_BOOL) function check_true(p_factory, flag)
    use, intrinsic :: iso_c_binding
    implicit none
    class(factory), value, intent(in) :: p_factory
    logical, value, intent(in) :: flag
    logical(C_BOOL) c_flag
    c_flag = flag
    check_true = factory_check_true(p_factory%p_factory, c_flag)
    return
  end function check_true
end module gridpack_factory
