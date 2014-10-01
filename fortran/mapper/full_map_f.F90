! ----------------------------------------------------------------
! file: full_map_f.F90
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
module gridpack_full_matrix_map
  use, intrinsic :: iso_c_binding
  use gridpack_network
  use gridpack_math
  implicit none
!
!  Define mapper types
!
  private
!
  type, public :: full_matrix_map
    type(C_PTR) :: p_mapper
    contains
    procedure::create
    procedure::destroy
    procedure::map_to_matrix
    procedure::remap_to_matrix
    procedure::overwrite_matrix
    procedure::increment_matrix
    procedure::check
  end type
!
  interface
!
! Create a full matrix map
! @param mapper pointer to Fortran full matrix map object
! @param network pointer to Fortran network object
!
    subroutine p_full_matrix_map_create(mapper, network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: mapper
      type(C_PTR), value, intent(in) :: network
    end subroutine p_full_matrix_map_create
!
! Destroy a full matrix map
! @param mapper pointer to Fortran full matrix map object
!
    subroutine p_full_matrix_map_destroy(mapper) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: mapper
    end subroutine p_full_matrix_map_destroy
!
! Create a matrix from the network
! @param mapper pointer to mapper
! @return pointer to Fortran matrix object
!
    type(C_PTR) function p_full_matrix_map_map_to_matrix(mapper) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: mapper
    end function p_full_matrix_map_map_to_matrix
!
! Update a matrix from the network
! @param mapper pointer to mapper
! @param matrix pointer to Fortran matrix object
!
    subroutine p_full_matrix_map_remap_to_matrix(mapper, matrix) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: mapper
      type(C_PTR), value, intent(in) :: matrix
    end subroutine p_full_matrix_map_remap_to_matrix
!
! Overwrite selected elements of existing matrix.
! @param mapper pointer to mapper
! @param matrix pointer to Fortran matrix object
!
    subroutine p_full_matrix_map_overwrite_matrix(mapper, matrix) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: mapper
      type(C_PTR), value, intent(in) :: matrix
    end subroutine p_full_matrix_map_overwrite_matrix
!
! Increment selected elements of existing matrix.
! @param mapper pointer to mapper
! @param matrix pointer to Fortran matrix object
!
    subroutine p_full_matrix_map_increment_matrix(mapper, matrix) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: mapper
      type(C_PTR), value, intent(in) :: matrix
    end subroutine p_full_matrix_map_increment_matrix
!
! Check to see if matrix looks well formed. This method runs through all
! branches and verifies that the dimensions of the branch contributions match
! the dimensions of the bus contributions at each end. If there is a
! discrepancy, an error message is generated.
! @param mapper pointer to mapper
! @return true if no discrepancy found
!
    logical(C_BOOL) function p_full_matrix_map_check(mapper) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: mapper
    end function p_full_matrix_map_check
  end interface
  contains
!
! Create a full matrix map
! @param p_mapper pointer to Fortran full matrix map object
! @param p_network pointer to Fortran network object
!
  subroutine create(p_mapper, p_network)
    use, intrinsic :: iso_c_binding
    use gridpack_network
    implicit none
    class(full_matrix_map), intent(inout) :: p_mapper
    class(network), intent(in) :: p_network
    call p_full_matrix_map_create(p_mapper%p_mapper, p_network%p_network)
    return
  end subroutine create
!
! Destroy a full matrix map
! @param p_mapper pointer to Fortran full matrix map object
!
  subroutine destroy(p_mapper)
    use, intrinsic :: iso_c_binding
    implicit none
    class(full_matrix_map), intent(inout) :: p_mapper
    call p_full_matrix_map_destroy(p_mapper%p_mapper)
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
    class(full_matrix_map), intent(in) :: p_mapper
    type(matrix) fmatrix
    fmatrix%mat = p_full_matrix_map_map_to_matrix(p_mapper%p_mapper)
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
    class(full_matrix_map), intent(in) :: p_mapper
    class(matrix), intent(in) :: p_matrix
    call p_full_matrix_map_remap_to_matrix(p_mapper%p_mapper, &
      p_matrix%mat)
    return
  end subroutine remap_to_matrix
!
! Overwrite selected elements of existing matrix.
! @param p_mapper pointer to mapper
! @param p_matrix pointer to Fortran matrix object
!
  subroutine overwrite_matrix(p_mapper, p_matrix)
    use, intrinsic :: iso_c_binding
    implicit none
    class(full_matrix_map), intent(in) :: p_mapper
    class(matrix), intent(in) :: p_matrix
    call p_full_matrix_map_overwrite_matrix(p_mapper%p_mapper, &
      p_matrix%mat)
    return
  end subroutine overwrite_matrix
!
! Increment selected elements of existing matrix.
! @param p_mapper pointer to mapper
! @param p_matrix pointer to Fortran matrix object
!
  subroutine increment_matrix(p_mapper, p_matrix)
    use, intrinsic :: iso_c_binding
    implicit none
    class(full_matrix_map), intent(in) :: p_mapper
    class(matrix), intent(in) :: p_matrix
    call p_full_matrix_map_increment_matrix(p_mapper%p_mapper, &
      p_matrix%mat)
    return
  end subroutine increment_matrix
!
! Check to see if matrix looks well formed. This method runs through all
! branches and verifies that the dimensions of the branch contributions match
! the dimensions of the bus contributions at each end. If there is a
! discrepancy, an error message is generated.
! @param p_mapper pointer to mapper
! @return true if no discrepancy found
!
  logical function check(p_mapper)
    use, intrinsic :: iso_c_binding
    implicit none
    class(full_matrix_map), intent(in) :: p_mapper
    logical(C_BOOL) ret
    ret = p_full_matrix_map_check(p_mapper%p_mapper)
    check = ret
    return
  end function check
end module gridpack_full_matrix_map
