!
!  Fortran mapper functions
!
module gridpack_mapper
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
    procedure::full_matrix_map_create
    procedure::full_matrix_map_destroy
    procedure::full_matrix_map_map_to_matrix
    procedure::full_matrix_map_remap_to_matrix
    procedure::full_matrix_map_overwrite_matrix
    procedure::full_matrix_map_increment_matrix
    procedure::full_matrix_map_check
  end type
!
  type, public :: bus_vector_map
    type(C_PTR) :: p_mapper
    contains
    procedure::bus_vector_map_create
    procedure::bus_vector_map_destroy
    procedure::bus_vector_map_map_to_vector
    procedure::bus_vector_map_remap_to_vector
    procedure::bus_vector_map_map_to_bus
  end type
!
  type, public :: gen_matrix_map
    type(C_PTR) :: p_mapper
    contains
    procedure::gen_matrix_map_create
    procedure::gen_matrix_map_destroy
    procedure::gen_matrix_map_map_to_matrix
    procedure::gen_matrix_map_remap_to_matrix
  end type
!
  type, public :: gen_vector_map
    type(C_PTR) :: p_mapper
    contains
    procedure::gen_vector_map_create
    procedure::gen_vector_map_destroy
    procedure::gen_vector_map_map_to_vector
    procedure::gen_vector_map_remap_to_vector
    procedure::gen_vector_map_map_to_bus
  end type
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
      type(C_PTR), intent(in) :: network
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
      type(C_PTR), intent(in) :: mapper
    end function p_full_matrix_map_map_to_matrix
!
! Update a matrix from the network
! @param mapper pointer to mapper
! @param matrix pointer to Fortran matrix object
!
    subroutine p_full_matrix_map_remap_to_matrix(mapper, matrix) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(in) :: mapper
      type(C_PTR), intent(in) :: matrix
    end subroutine p_full_matrix_map_remap_to_matrix
!
! Overwrite selected elements of existing matrix.
! @param mapper pointer to mapper
! @param matrix pointer to Fortran matrix object
!
    subroutine p_full_matrix_map_overwrite_matrix(mapper, matrix) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(in) :: mapper
      type(C_PTR), intent(in) :: matrix
    end subroutine p_full_matrix_map_overwrite_matrix
!
! Increment selected elements of existing matrix.
! @param mapper pointer to mapper
! @param matrix pointer to Fortran matrix object
!
    subroutine p_full_matrix_map_increment_matrix(mapper, matrix) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(in) :: mapper
      type(C_PTR), intent(in) :: matrix
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
      type(C_PTR), intent(in) :: mapper
    end function p_full_matrix_map_check
!
! Create a bus vector map
! @param mapper pointer to Fortran bus vector map object
! @param network pointer to Fortran network object
!
    subroutine p_bus_vector_map_create(mapper, network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: mapper
      type(C_PTR), intent(in) :: network
    end subroutine p_bus_vector_map_create
!
! Destroy a bus vector map
! @param mapper pointer to Fortran bus vector map object
!
    subroutine p_bus_vector_map_destroy(mapper) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: mapper
    end subroutine p_bus_vector_map_destroy
!
! Create a vector from the current bus state
! @param mapper pointer to mapper
! @return pointer to Fortran vector object
!
    type(C_PTR) function p_bus_vector_map_map_to_vector(mapper) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(in) :: mapper
    end function p_bus_vector_map_map_to_vector
!
! Reset a vector from the current bus state (vector should be created with the
! same mapper)
! @param mapper pointer to mapper
! @param vector pointer to Fortran vector object
!
    subroutine p_bus_vector_map_remap_to_vector(mapper, vector) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(in) :: mapper
      type(C_PTR), intent(in) :: vector
    end subroutine p_bus_vector_map_remap_to_vector
!
! Push data from a vector to the network buses
! @param mapper pointer to mapper
! @param vector pointer to Fortran vector object
!
    subroutine p_bus_vector_map_map_to_bus(mapper, vector) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(in) :: mapper
      type(C_PTR), intent(in) :: vector
    end subroutine p_bus_vector_map_map_to_bus
!
! Create a generic matrix map
! @param mapper pointer to Fortran generic matrix map object
! @param network pointer to Fortran network object
!
    subroutine p_gen_matrix_map_create(mapper, network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: mapper
      type(C_PTR), intent(in) :: network
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
      type(C_PTR), intent(in) :: mapper
    end function p_gen_matrix_map_map_to_matrix
!
! Update a matrix from the network
! @param mapper pointer to mapper
! @param matrix pointer to Fortran matrix object
!
    subroutine p_gen_matrix_map_remap_to_matrix(mapper, matrix) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(in) :: mapper
      type(C_PTR), intent(in) :: matrix
    end subroutine p_gen_matrix_map_remap_to_matrix
!
! Create a generic vector map
! @param mapper pointer to Fortran generic vector map object
! @param network pointer to Fortran network object
!
    subroutine p_gen_vector_map_create(mapper, network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: mapper
      type(C_PTR), intent(in) :: network
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
      type(C_PTR), intent(in) :: mapper
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
      type(C_PTR), intent(in) :: mapper
      type(C_PTR), intent(in) :: vector
    end subroutine p_gen_vector_map_remap_to_vector
!
! Push data from a vector to the network buses
! @param mapper pointer to mapper
! @param vector pointer to Fortran vector object
!
    subroutine p_gen_vector_map_map_to_bus(mapper, vector) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(in) :: mapper
      type(C_PTR), intent(in) :: vector
    end subroutine p_gen_vector_map_map_to_bus
  end interface
  contains
!
! Create a full matrix map
! @param p_mapper pointer to Fortran full matrix map object
! @param p_network pointer to Fortran network object
!
  subroutine full_matrix_map_create(p_mapper, p_network)
    use, intrinsic :: iso_c_binding
    use gridpack_network
    implicit none
    class(full_matrix_map), intent(inout) :: p_mapper
    class(network), intent(in) :: p_network
    call p_full_matrix_map_create(p_mapper%p_mapper, p_network%p_network)
    return
  end subroutine full_matrix_map_create
!
! Destroy a full matrix map
! @param p_mapper pointer to Fortran full matrix map object
!
  subroutine full_matrix_map_destroy(p_mapper)
    use, intrinsic :: iso_c_binding
    implicit none
    class(full_matrix_map), intent(inout) :: p_mapper
    call p_full_matrix_map_destroy(p_mapper%p_mapper)
    return
  end subroutine full_matrix_map_destroy
!
! Create a matrix from the network
! @param p_mapper pointer to mapper
! @return pointer to Fortran matrix object
!
  function full_matrix_map_map_to_matrix(p_mapper) result(fmatrix)
    use, intrinsic :: iso_c_binding
    implicit none
    class(full_matrix_map), intent(in) :: p_mapper
    type(matrix) fmatrix
    fmatrix%mat = p_full_matrix_map_map_to_matrix(p_mapper%p_mapper)
    return
  end function full_matrix_map_map_to_matrix
!
! Update a matrix from the network
! @param p_mapper pointer to mapper
! @param p_matrix pointer to Fortran matrix object
!
  subroutine full_matrix_map_remap_to_matrix(p_mapper, p_matrix)
    use, intrinsic :: iso_c_binding
    implicit none
    class(full_matrix_map), intent(in) :: p_mapper
    class(matrix), intent(in) :: p_matrix
    call p_full_matrix_map_remap_to_matrix(p_mapper%p_mapper, &
      p_matrix%mat)
    return
  end subroutine full_matrix_map_remap_to_matrix
!
! Overwrite selected elements of existing matrix.
! @param p_mapper pointer to mapper
! @param p_matrix pointer to Fortran matrix object
!
  subroutine full_matrix_map_overwrite_matrix(p_mapper, p_matrix)
    use, intrinsic :: iso_c_binding
    implicit none
    class(full_matrix_map), intent(in) :: p_mapper
    class(matrix), intent(in) :: p_matrix
    call p_full_matrix_map_overwrite_matrix(p_mapper%p_mapper, &
      p_matrix%mat)
    return
  end subroutine full_matrix_map_overwrite_matrix
!
! Increment selected elements of existing matrix.
! @param p_mapper pointer to mapper
! @param p_matrix pointer to Fortran matrix object
!
  subroutine full_matrix_map_increment_matrix(p_mapper, p_matrix)
    use, intrinsic :: iso_c_binding
    implicit none
    class(full_matrix_map), intent(in) :: p_mapper
    class(matrix), intent(in) :: p_matrix
    call p_full_matrix_map_increment_matrix(p_mapper%p_mapper, &
      p_matrix%mat)
    return
  end subroutine full_matrix_map_increment_matrix
!
! Check to see if matrix looks well formed. This method runs through all
! branches and verifies that the dimensions of the branch contributions match
! the dimensions of the bus contributions at each end. If there is a
! discrepancy, an error message is generated.
! @param p_mapper pointer to mapper
! @return true if no discrepancy found
!
  logical function full_matrix_map_check(p_mapper)
    use, intrinsic :: iso_c_binding
    implicit none
    class(full_matrix_map), intent(in) :: p_mapper
    logical(C_BOOL) ret
    ret = p_full_matrix_map_check(p_mapper%p_mapper)
    full_matrix_map_check = ret
    return
  end function full_matrix_map_check
!
! Create a bus vector map
! @param p_mapper pointer to Fortran bus vector map object
! @param p_network pointer to Fortran network object
!
  subroutine bus_vector_map_create(p_mapper, p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_vector_map), intent(inout) :: p_mapper
    class(network), intent(in) :: p_network
    call p_bus_vector_map_create(p_mapper%p_mapper, p_network%p_network)
    return
  end subroutine bus_vector_map_create
!
! Destroy a bus vector map
! @param p_mapper pointer to Fortran bus vector map object
!
  subroutine bus_vector_map_destroy(p_mapper)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_vector_map), intent(inout) :: p_mapper
    call p_bus_vector_map_destroy(p_mapper%p_mapper)
    return
  end subroutine bus_vector_map_destroy
!
! Create a vector from the current bus state
! @param p_mapper pointer to mapper
! @return pointer to Fortran vector object
!
  function bus_vector_map_map_to_vector(p_mapper) result(fvector)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_vector_map), intent(in) :: p_mapper
    type(vector) :: fvector
    fvector%vec = p_bus_vector_map_map_to_vector(p_mapper%p_mapper)
    return
  end function bus_vector_map_map_to_vector
!
! Reset a vector from the current bus state (vector should be created with the
! same mapper)
! @param p_mapper pointer to mapper
! @param p_vector pointer to Fortran vector object
!
  subroutine bus_vector_map_remap_to_vector(p_mapper, p_vector)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_vector_map), intent(in) :: p_mapper
    class(vector), intent(in) :: p_vector
    call p_bus_vector_map_remap_to_vector(p_mapper%p_mapper,p_vector%vec)
    return
  end subroutine bus_vector_map_remap_to_vector
!
! Push data from a vector to the network buses
! @param p_mapper pointer to mapper
! @param p_vector pointer to Fortran vector object
!
  subroutine bus_vector_map_map_to_bus(p_mapper, p_vector)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_vector_map), intent(in) :: p_mapper
    class(vector), intent(in) :: p_vector
    call p_bus_vector_map_map_to_bus(p_mapper%p_mapper,p_vector%vec)
    return
  end subroutine bus_vector_map_map_to_bus
!
! Create a generic matrix map
! @param p_mapper pointer to Fortran generic matrix map object
! @param p_network pointer to Fortran network object
!
  subroutine gen_matrix_map_create(p_mapper, p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_matrix_map), intent(inout) :: p_mapper
    class(network), intent(in) :: p_network
    call p_gen_matrix_map_create(p_mapper%p_mapper, p_network%p_network)
  end subroutine gen_matrix_map_create
!
! Destroy a generic matrix map
! @param p_mapper pointer to Fortran generic matrix map object
!
  subroutine gen_matrix_map_destroy(p_mapper)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_matrix_map), intent(inout) :: p_mapper
    call p_gen_matrix_map_destroy(p_mapper%p_mapper)
    return
  end subroutine gen_matrix_map_destroy
!
! Create a matrix from the network
! @param p_mapper pointer to mapper
! @return pointer to Fortran matrix object
!
  function gen_matrix_map_map_to_matrix(p_mapper) result(fmatrix)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_matrix_map), intent(in) :: p_mapper
    type(matrix) fmatrix
    fmatrix%mat = p_gen_matrix_map_map_to_matrix(p_mapper%p_mapper)
    return
  end function gen_matrix_map_map_to_matrix
!
! Update a matrix from the network
! @param p_mapper pointer to mapper
! @param p_matrix pointer to Fortran matrix object
!
  subroutine gen_matrix_map_remap_to_matrix(p_mapper, p_matrix)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_matrix_map), intent(in) :: p_mapper
    class(matrix), intent(in) :: p_matrix
    call p_gen_matrix_map_remap_to_matrix(p_mapper%p_mapper, &
       p_matrix%mat)
    return
  end subroutine gen_matrix_map_remap_to_matrix
!
! Create a generic vector map
! @param p_mapper pointer to Fortran generic vector map object
! @param p_network pointer to Fortran network object
!
  subroutine gen_vector_map_create(p_mapper, p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_vector_map), intent(inout) :: p_mapper
    class(network), intent(in) :: p_network
    call p_gen_vector_map_create(p_mapper%p_mapper, p_network%p_network)
    return
  end subroutine gen_vector_map_create
!
! Destroy a generic vector map
! @param p_mapper pointer to Fortran generic vector map object
!
  subroutine gen_vector_map_destroy(p_mapper)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_vector_map), intent(inout) :: p_mapper
    call p_gen_vector_map_destroy(p_mapper%p_mapper)
    return
  end subroutine gen_vector_map_destroy
!
! Create a vector from the current bus state
! @param p_mapper pointer to mapper
! @return pointer to Fortran vector object
!
  function gen_vector_map_map_to_vector(p_mapper) result(fvector)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_vector_map), intent(in) :: p_mapper
    type(vector) fvector
    fvector%vec = p_gen_vector_map_map_to_vector(p_mapper%p_mapper)
    return
  end function gen_vector_map_map_to_vector
!
! Reset a vector from the current bus state (vector should be created with the
! same mapper)
! @param p_mapper pointer to mapper
! @param p_vector pointer to Fortran vector object
!
  subroutine gen_vector_map_remap_to_vector(p_mapper, p_vector)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_vector_map), intent(in) :: p_mapper
    class(vector), intent(in) :: p_vector
    call p_gen_vector_map_remap_to_vector(p_mapper%p_mapper, p_vector%vec)
    return
  end subroutine gen_vector_map_remap_to_vector
!
! Push data from a vector to the network buses
! @param p_mapper pointer to mapper
! @param p_vector pointer to Fortran vector object
!
  subroutine gen_vector_map_map_to_bus(p_mapper, p_vector)
    use, intrinsic :: iso_c_binding
    implicit none
    class(gen_vector_map), intent(in) :: p_mapper
    class(vector), intent(in) :: p_vector
    call p_gen_vector_map_map_to_bus(p_mapper%p_mapper, p_vector%vec)
    return
  end subroutine gen_vector_map_map_to_bus
end module gridpack_mapper
