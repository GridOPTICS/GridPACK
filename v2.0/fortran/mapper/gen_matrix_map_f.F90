! ----------------------------------------------------------------
! file: gen_matrix_map_f.F90
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
module gridpack_gen_matrix_map
  use, intrinsic :: iso_c_binding
  use gridpack_network
  use gridpack_math
  implicit none
!
!  Define mapper types
!
  private
!
  type, public :: gen_matrix_map
    type(C_PTR) :: p_mapper
    contains
    procedure::create
    procedure::destroy
    procedure::map_to_matrix
    procedure::remap_to_matrix
  end type
!
  interface
!
! Create a generic matrix map
! @param mapper pointer to Fortran generic matrix map object
! @param network pointer to Fortran network object
!
    subroutine p_gen_matrix_map_create(mapper, network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: mapper
      type(C_PTR), value, intent(in) :: network
    end subroutine p_gen_matrix_map_create
!
! Destroy a generic matrix map
! @param mapper pointer to Fortran generic matrix map object
!
    subroutine p_gen_matrix_map_destroy(mapper) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: mapper
    end subroutine p_gen_matrix_map_destroy
!
! Create a matrix from the network
! @param mapper pointer to mapper
! @return pointer to Fortran matrix object
!
    type(C_PTR) function p_gen_matrix_map_map_to_matrix(mapper) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: mapper
    end function p_gen_matrix_map_map_to_matrix
!
! Update a matrix from the network
! @param mapper pointer to mapper
! @param matrix pointer to Fortran matrix object
!
    subroutine p_gen_matrix_map_remap_to_matrix(mapper, matrix) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: mapper
      type(C_PTR), value, intent(in) :: matrix
    end subroutine p_gen_matrix_map_remap_to_matrix
  end interface
  contains
!
! Create a generic matrix map
! @param p_mapper pointer to Fortran generic matrix map object
! @param p_network pointer to Fortran network object
!
  subroutine create(p_mapper, p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_matrix_map), intent(inout) :: p_mapper
    class(network), intent(in) :: p_network
    call p_gen_matrix_map_create(p_mapper%p_mapper, p_network%p_network)
  end subroutine create
!
! Destroy a generic matrix map
! @param p_mapper pointer to Fortran generic matrix map object
!
  subroutine destroy(p_mapper)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_matrix_map), intent(inout) :: p_mapper
    call p_gen_matrix_map_destroy(p_mapper%p_mapper)
    return
  end subroutine destroy
!
! Create a matrix from the network
! @param p_mapper pointer to mapper
! @return pointer to Fortran matrix object
!
  function map_to_matrix(p_mapper) result(fmatrix)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_matrix_map), intent(in) :: p_mapper
    type(matrix) fmatrix
    fmatrix%mat = p_gen_matrix_map_map_to_matrix(p_mapper%p_mapper)
    return
  end function map_to_matrix
!
! Update a matrix from the network
! @param p_mapper pointer to mapper
! @param p_matrix pointer to Fortran matrix object
!
  subroutine remap_to_matrix(p_mapper, p_matrix)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_matrix_map), intent(in) :: p_mapper
    class(matrix), intent(in) :: p_matrix
    call p_gen_matrix_map_remap_to_matrix(p_mapper%p_mapper, &
       p_matrix%mat)
    return
  end subroutine remap_to_matrix
end module gridpack_gen_matrix_map
