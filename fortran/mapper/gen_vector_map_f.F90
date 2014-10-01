! ----------------------------------------------------------------
! file: gen_vector_map_f.F90
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
!  Fortran mapper functions
!
module gridpack_gen_vector_map
  use, intrinsic :: iso_c_binding
  use gridpack_network
  use gridpack_math
  implicit none
!
!  Define mapper types
!
  private
!
  type, public :: gen_vector_map
    type(C_PTR) :: p_mapper
    contains
    procedure::create
    procedure::destroy
    procedure::map_to_vector
    procedure::remap_to_vector
    procedure::map_to_bus
  end type
!
  interface
!
! Create a generic vector map
! @param mapper pointer to Fortran generic vector map object
! @param network pointer to Fortran network object
!
    subroutine p_gen_vector_map_create(mapper, network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: mapper
      type(C_PTR), value, intent(in) :: network
    end subroutine p_gen_vector_map_create
!
! Destroy a generic vector map
! @param mapper pointer to Fortran generic vector map object
!
    subroutine p_gen_vector_map_destroy(mapper) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: mapper
    end subroutine p_gen_vector_map_destroy
!
! Create a vector from the current bus state
! @param mapper pointer to mapper
! @return pointer to Fortran vector object
!
    type(C_PTR) function p_gen_vector_map_map_to_vector(mapper) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: mapper
    end function p_gen_vector_map_map_to_vector
!
! Reset a vector from the current bus state (vector should be created with the
! same mapper)
! @param mapper pointer to mapper
! @param vector pointer to Fortran vector object
!
    subroutine p_gen_vector_map_remap_to_vector(mapper, vector) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: mapper
      type(C_PTR), value, intent(in) :: vector
    end subroutine p_gen_vector_map_remap_to_vector
!
! Push data from a vector to the network buses
! @param mapper pointer to mapper
! @param vector pointer to Fortran vector object
!
    subroutine p_gen_vector_map_map_to_bus(mapper, vector) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: mapper
      type(C_PTR), value, intent(in) :: vector
    end subroutine p_gen_vector_map_map_to_bus
  end interface
  contains
!
! Create a generic vector map
! @param p_mapper pointer to Fortran generic vector map object
! @param p_network pointer to Fortran network object
!
  subroutine create(p_mapper, p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_vector_map), intent(inout) :: p_mapper
    class(network), intent(in) :: p_network
    call p_gen_vector_map_create(p_mapper%p_mapper, p_network%p_network)
    return
  end subroutine create
!
! Destroy a generic vector map
! @param p_mapper pointer to Fortran generic vector map object
!
  subroutine destroy(p_mapper)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_vector_map), intent(inout) :: p_mapper
    call p_gen_vector_map_destroy(p_mapper%p_mapper)
    return
  end subroutine destroy
!
! Create a vector from the current bus state
! @param p_mapper pointer to mapper
! @return pointer to Fortran vector object
!
  function map_to_vector(p_mapper) result(fvector)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_vector_map), intent(in) :: p_mapper
    type(vector) fvector
    fvector%vec = p_gen_vector_map_map_to_vector(p_mapper%p_mapper)
    return
  end function map_to_vector
!
! Reset a vector from the current bus state (vector should be created with the
! same mapper)
! @param p_mapper pointer to mapper
! @param p_vector pointer to Fortran vector object
!
  subroutine remap_to_vector(p_mapper, p_vector)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_vector_map), intent(in) :: p_mapper
    class(vector), intent(in) :: p_vector
    call p_gen_vector_map_remap_to_vector(p_mapper%p_mapper, p_vector%vec)
    return
  end subroutine remap_to_vector
!
! Push data from a vector to the network buses
! @param p_mapper pointer to mapper
! @param p_vector pointer to Fortran vector object
!
  subroutine map_to_bus(p_mapper, p_vector)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_vector_map), intent(in) :: p_mapper
    class(vector), intent(in) :: p_vector
    call p_gen_vector_map_map_to_bus(p_mapper%p_mapper, p_vector%vec)
    return
  end subroutine map_to_bus
end module gridpack_gen_vector_map
