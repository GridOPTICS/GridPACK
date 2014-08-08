!
! Fortran network functions
!
module gridpack_network
  use, intrinsic :: iso_c_binding
!
! Define network type
!
  implicit none
  
  private

  type, public :: network
    type (c_ptr) :: p_network
    contains
    procedure::create_network
!    procedure::add_bus
!    procedure::add_branch
  end type
  interface
!
! Create new network and return handle 
! @param comm GridPACK communicator
! @return handle for new network
!
    subroutine p_create_network(network,comm)
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: network
      type(c_ptr), intent(in) :: comm
    end subroutine p_create_network
  end interface
  contains
!
! Create new network and return handle 
! @param comm GridPACK communicator
! @return handle for new network
!
  subroutine p_create_network(network,comm)
    use gridpack_communicator
    use, intrinsic :: iso_c_binding
    implicit none
    type(C_PTR), intent(inout) :: network
    type(C_PTR), value, intent(in) :: comm
    p_create_network(network%p_network,comm%comm)
    return
  end function create_network
#if 0
!
! Add a bus locally to the network
! @param n_handle network handle
! @param idx original index of bus
!
  subroutine add_bus(idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer(C_INT) c_handle, c_idx
    c_handle = n_handle
    c_idx = idx
    call p_add_bus(p_network,c_idx)
    return
  end subroutine add_bus
!
! Add a branch locally to the network
! @param n_handle network handle
! @param idx1 original index of "from" bus
! @param idx2 original index of "to" bus
!
  subroutine add_branch(n_handle, idx1, idx2)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx1, idx2
    integer(C_INT) c_handle, c_idx1, c_idx2
    c_handle = n_handle
    c_idx1 = idx1
    c_idx2 = idx2
    call p_add_branch(c_handle, c_idx1, c_idx2)
    return
  end subroutine add_branch
!
! Number of local buses (both active and inactive) on processor
! @param n_handle network handle
! @return number of local buses
! 
  integer function num_buses(n_handle)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle
    integer(C_INT) c_handle, nbus
    c_handle = n_handle
    nbus = p_num_buses(c_handle)
    num_buses = nbus
    return
  end function num_buses
!
! Return the total number of buses in the entire network
! @param n_handle network handle
! @return total number of buses in network
! 
  integer function total_buses(n_handle)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle
    integer(C_INT) c_handle, nbus
    c_handle = n_handle
    nbus = p_total_buses(c_handle)
    total_buses = nbus
    return
  end function total_buses
!
! Number of local branches (both active and inactive) on processor
! @param n_handle network handle
! @return number of local branches
! 
  integer function num_branches(n_handle)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle
    integer(C_INT) c_handle, nbranch
    c_handle = n_handle
    nbranch = p_num_branches(c_handle)
    num_branches = nbranch
    return
  end function num_branches
!
! Return the total number of branches in the entire network
! @param n_handle network handle
! @return total number of branches in network
! 
  integer function total_branches(n_handle)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle
    integer(C_INT) c_handle, nbranch
    c_handle = n_handle
    nbranch = p_total_branches(c_handle)
    total_branches = nbranch
    return
  end function total_branches
!
! Designate a bus as a reference bus
! @param n_handle network handle
! @param idx local bus index
!
  subroutine set_reference_bus(n_handle, idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer(C_INT) c_handle, c_idx
    c_handle = n_handle
    c_idx = idx
    call p_set_reference_bus(c_handle, c_idx)
    return
  end subroutine set_reference_bus
!
! Return the index of reference bus. Return -1 if reference bus is not on this
! processor
! @param n_handle network handle
! @return local index of reference bus
!
  integer function get_reference_bus(n_handle)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle
    integer(C_INT) c_handle, c_idx
    c_handle = n_handle
    c_idx = p_get_reference_bus(c_handle)
    get_reference_bus = c_idx
    return
  end function get_reference_bus
!
! Set the original index of the bus (from configuration file)
! @param n_handle network handle
! @param idx local bus index
! @param o_idx original bus index
! @return false if bus does not exist at that value of idx
!
  logical function set_original_bus_index(n_handle, idx, o_idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer, intent(out) :: o_idx
    integer(C_INT) c_handle, c_idx, c_o_idx
    logical(C_BOOL) c_flag
    c_handle = n_handle
    c_idx = idx
    c_flag = p_set_original_bus_index(c_handle, c_idx, c_o_idx)
    set_original_bus_index = c_flag
    o_idx = c_o_idx
    return
  end function set_original_bus_index
!
! Set the global index of the bus
! @param n_handle network handle
! @param idx local bus index
! @param o_idx global bus index
! @return false if bus does not exist at that value of idx
!
  logical function set_global_bus_index(n_handle, idx, g_idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer, intent(out) :: g_idx
    integer(C_INT) c_handle, c_idx, c_g_idx
    logical(C_BOOL) c_flag
    c_handle = n_handle
    c_idx = idx
    c_flag = p_set_global_bus_index(c_handle, c_idx, c_g_idx)
    set_global_bus_index = c_flag
    g_idx = c_g_idx
    return
  end function set_global_bus_index
!
! Set the global index of the branch
! @param n_handle network handle
! @param idx local branch index
! @param o_idx global branch index
! @return false if branch does not exist at that value of idx
!
  logical function set_global_branch_index(n_handle, idx, g_idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer, intent(out) :: g_idx
    integer(C_INT) c_handle, c_idx, c_g_idx
    logical(C_BOOL) c_flag
    c_handle = n_handle
    c_idx = idx
    c_flag = p_set_global_branch_index(c_handle, c_idx, c_g_idx)
    set_global_branch_index = c_flag
    g_idx = c_g_idx
    return
  end function set_global_branch_index
!
! Set the origina1 index of the bus at the "from" end of branch
! @param n_handle network handle
! @param idx local branch index
! @param b_idx original index of "from" bus
! @return false if branch does not exist at that value of idx
!
  logical function set_original_bus_index1(n_handle, idx, b_idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx, b_idx
    integer(C_INT) c_handle, c_idx, c_b_idx
    logical(C_BOOL) c_flag
    c_handle = n_handle
    c_idx = idx
    c_b_idx = b_idx
    c_flag = p_set_original_bus_index1(c_handle, c_idx, c_b_idx)
    set_original_bus_index1 = c_flag
    return
  end function set_original_bus_index1
!
! Set the origina1 index of the bus at the "to" end of branch
! @param n_handle network handle
! @param idx local branch index
! @param b_idx original index of "to" bus
! @return false if branch does not exist at that value of idx
!
  logical function set_original_bus_index2(n_handle, idx, b_idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx, b_idx
    integer(C_INT) c_handle, c_idx, c_b_idx
    logical(C_BOOL) c_flag
    c_handle = n_handle
    c_idx = idx
    c_b_idx = b_idx
    c_flag = p_set_original_bus_index2(c_handle, c_idx, c_b_idx)
    set_original_bus_index2 = c_flag
    return
  end function set_original_bus_index2
!
! Set the globa1 index of the bus at the "from" end of branch
! @param n_handle network handle
! @param idx local branch index
! @param b_idx global index of "from" bus
! @return false if branch does not exist at that value of idx
!
  logical function set_global_bus_index1(n_handle, idx, b_idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx, b_idx
    integer(C_INT) c_handle, c_idx, c_b_idx
    logical(C_BOOL) c_flag
    c_handle = n_handle
    c_idx = idx
    c_b_idx = b_idx
    c_flag = p_set_global_bus_index1(c_handle, c_idx, c_b_idx)
    set_global_bus_index1 = c_flag
    return
  end function set_global_bus_index1
!
! Set the globa1 index of the bus at the "to" end of branch
! @param n_handle network handle
! @param idx local branch index
! @param b_idx global index of "to" bus
! @return false if branch does not exist at that value of idx
!
  logical function set_global_bus_index2(n_handle, idx, b_idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx, b_idx
    integer(C_INT) c_handle, c_idx, c_b_idx
    logical(C_BOOL) c_flag
    c_handle = n_handle
    c_idx = idx
    c_b_idx = b_idx
    c_flag = p_set_global_bus_index2(c_handle, c_idx, c_b_idx)
    set_global_bus_index2 = c_flag
    return
  end function set_global_bus_index2
!
! Set the loca1 index of the bus at the "from" end of branch
! @param n_handle network handle
! @param idx local branch index
! @param b_idx local index of "from" bus
! @return false if branch does not exist at that value of idx
!
  logical function set_local_bus_index1(n_handle, idx, b_idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx, b_idx
    integer(C_INT) c_handle, c_idx, c_b_idx
    logical(C_BOOL) c_flag
    c_handle = n_handle
    c_idx = idx
    c_b_idx = b_idx
    c_flag = p_set_local_bus_index1(c_handle, c_idx, c_b_idx)
    set_local_bus_index1 = c_flag
    return
  end function set_local_bus_index1
!
! Set the loca1 index of the bus at the "to" end of branch
! @param n_handle network handle
! @param idx local branch index
! @param b_idx local index of "to" bus
! @return false if branch does not exist at that value of idx
!
  logical function set_local_bus_index2(n_handle, idx, b_idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx, b_idx
    integer(C_INT) c_handle, c_idx, c_b_idx
    logical(C_BOOL) c_flag
    c_handle = n_handle
    c_idx = idx
    c_b_idx = b_idx
    c_flag = p_set_local_bus_index2(c_handle, c_idx, c_b_idx)
    set_local_bus_index2 = c_flag
    return
  end function set_local_bus_index2
!
! Set the active flag of the bus
! @param n_handle network handle
! @param idx local bus index
! @param flag active status of bus
! @return false if bus does not exist at that value of idx
! 
  logical function set_active_bus(n_handle, idx, flag)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    logical flag
    integer(C_INT) c_handle, c_idx
    logical(C_BOOL) c_flag, c_ret
    c_handle = n_handle
    c_idx = idx
    c_flag = flag
    c_ret = p_set_active_bus(c_handle, c_idx, c_flag)
    set_active_bus = c_ret
    return
  end function set_active_bus
!
! Set the active flag of the branch
! @param n_handle network handle
! @param idx local branch index
! @param flag active status of branch
! @return false if branch does not exist at that value of idx
! 
  logical function set_active_branch(n_handle, idx, flag)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    logical flag
    integer(C_INT) c_handle, c_idx
    logical(C_BOOL) c_flag, c_ret
    c_handle = n_handle
    c_idx = idx
    c_flag = flag
    c_ret = p_set_active_branch(c_handle, c_idx, c_flag)
    set_active_branch = c_ret
    return
  end function set_active_branch
!
! Clear the list of neighbors for the bus at idx
! @param n_handle network handle
! @param idx local branch index
! @return false if bus does not exist at that value of idx
! 
  logical function clear_branch_neighbors(n_handle, idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer(C_INT) c_handle, c_idx
    logical(C_BOOL) c_ret
    c_handle = n_handle
    c_idx = idx
    c_ret = p_clear_branch_neighbors(c_handle, c_idx)
    clear_branch_neighbors = c_ret
    return
  end function clear_branch_neighbors
!
! Add a local index for a branch attached to the bus at idx
! @param n_handle network handle
! @param idx local bus index
! @param br_idx local branch index
! @return false if bus does not exist at that value of idx
! 
  logical function add_branch_neighbor(n_handle, idx, br_idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx, br_idx
    integer(C_INT) c_handle, c_idx, c_br_idx
    logical(C_BOOL) c_ret
    c_handle = n_handle
    c_idx = idx
    c_br_idx = br_idx
    c_ret = p_add_branch_neighbor(c_handle, c_idx, c_br_idx)
    add_branch_neighbor = c_ret
    return
  end function add_branch_neighbor
!
! Get status of the bus (local or ghosted)
! @param n_handle network handle
! @param idx local bus index
! @return status of bus
!
  logical function get_active_bus(n_handle, idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer(C_INT) c_handle, c_idx
    logical(C_BOOL) c_ret
    c_handle = n_handle
    c_idx = idx
    c_ret = p_get_active_bus(c_handle, c_idx)
    get_active_bus = c_ret
    return
  end function get_active_bus
!
! Get original index of bus
! @param n_handle network handle
! @param idx local bus index
! @return original index of bus
!
  integer function get_original_bus_index(n_handle, idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer(C_INT) c_handle, c_idx
    integer(C_INT) c_ret
    c_handle = n_handle
    c_idx = idx
    c_ret = p_get_original_bus_index(c_handle, c_idx)
    get_original_bus_index = c_ret
    return
  end function get_original_bus_index
!
! Get global index of bus
! @param n_handle network handle
! @param idx local bus index
! @return global index of bus
!
  integer function get_global_bus_index(n_handle, idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer(C_INT) c_handle, c_idx
    integer(C_INT) c_ret
    c_handle = n_handle
    c_idx = idx
    c_ret = p_get_global_bus_index(c_handle, c_idx)
    get_global_bus_index = c_ret
    return
  end function get_global_bus_index
!
! Get status of the branch (local or ghosted)
! @param n_handle network handle
! @param idx local branch index
! @return status of branch
!
  logical function get_active_branch(n_handle, idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer(C_INT) c_handle, c_idx
    logical(C_BOOL) c_ret
    c_handle = n_handle
    c_idx = idx
    c_ret = p_get_active_branch(c_handle, c_idx)
    get_active_branch = c_ret
    return
  end function get_active_branch
!
! Get global index of branch
! @param n_handle network handle
! @param idx local branch index
! @return global index of bus
!
  integer function get_global_branch_index(n_handle, idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer(C_INT) c_handle, c_idx
    integer(C_INT) c_ret
    c_handle = n_handle
    c_idx = idx
    c_ret = p_get_global_branch_index(c_handle, c_idx)
    get_global_branch_index = c_ret
    return
  end function get_global_branch_index
!
! Get original indices of buses at each end of the branch
! @param n_handle network handle
! @param idx local branch index
! @param idx1 original index of "from" bus
! @param idx2 original index of "to" bus
!
  subroutine get_original_branch_endpoints(n_handle, idx, idx1, idx2)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer, intent(out) :: idx1, idx2
    integer(C_INT) c_handle, c_idx
    integer(C_INT) c_idx1, c_idx2
    c_handle = n_handle
    c_idx = idx
    call p_get_original_branch_endpoints(c_handle, c_idx, c_idx1, c_idx2)
    idx1 = c_idx1
    idx2 = c_idx2
    return
  end subroutine get_original_branch_endpoints
!
! Return the number of branches connected to a bus. This can be used to allocate
! arrays to hold the branch indices
! @param n_handle network handle
! @param idx local bus index
! @return number of connected branches
!
  integer function get_num_connected_branches(n_handle, idx)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer(C_INT) c_handle, c_idx
    integer(C_INT) c_ret
    c_handle = n_handle
    c_idx = idx
    c_ret = p_get_num_connected_branches(c_handle, c_idx)
    get_num_connected_branches = c_ret
    return
  end function get_num_connected_branches
!
! Return list of branches connected to bus
! @param n_handle network handle
! @param idx local bus index
! @param branches list of local branch indices
!
  subroutine get_connected_branches(n_handle, idx, branches)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer, intent(out) :: branches(*)
    integer(C_INT) c_handle, c_idx, c_num
    integer(C_INT), allocatable :: c_branches(:)
    integer i
    c_handle = n_handle
    c_idx = idx
    c_num = p_get_num_connected_branches(c_handle, c_idx)
    allocate(c_branches(c_num))
    call p_get_connected_branches(c_handle, c_idx, c_branches)
    do i = 1, c_num
      branches(i) = c_branches(i)
    end do
    deallocate(c_branches)
    return
  end subroutine get_connected_branches
!
! Return list of buses connected to bus via one branch
! @param n_handle network handle
! @param idx local bus index
! @param buses list of local bus indices
!
  subroutine get_connected_buses(n_handle, idx, buses)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer, intent(out) :: buses(*)
    integer(C_INT) c_handle, c_idx, c_num
    integer(C_INT), allocatable :: c_buses(:)
    integer i
    c_handle = n_handle
    c_idx = idx
    c_num = p_get_num_connected_branches(c_handle, c_idx)
    allocate(c_buses(c_num))
    call p_get_connected_buses(c_handle, c_idx, c_buses)
    do i = 1, c_num
      buses(i) = c_buses(i)
    end do
    deallocate(c_buses)
    return
  end subroutine get_connected_buses
!
! Return the indices of the buses at either end of the branch
! @param n_handle network handle
! @param idx local branch index
! @param idx1 local index of "from" bus
! @param idx2 local index of "to" bus
!
  subroutine get_branch_endpoints(n_handle, idx, idx1, idx2)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer, intent(out) :: idx1, idx2
    integer(C_INT) c_handle, c_idx, c_idx1, c_idx2
    c_handle = n_handle
    c_idx = idx
    call p_get_branch_endpoints(c_handle, c_idx, c_idx1, c_idx2)
    idx1 = c_idx1
    idx2 = c_idx2
    return
  end subroutine get_branch_endpoints
!
! Partition network over available processes
! @param n_handle network handle
!
  subroutine partition_network(n_handle)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle
    integer(C_INT) c_handle
    c_handle = n_handle
    call p_partition_network(c_handle)
    return
  end subroutine partition_network
!
! Clean all ghost buses and branches from the system. This can be used before
! repartioning the network
! @param n_handle network handle
!
  subroutine clean_network(n_handle)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle
    integer(C_INT) c_handle
    c_handle = n_handle
    call p_clean_network(c_handle)
    return
  end subroutine clean_network
!
! Store the location of externally allocated buffer within a network
! @param n_handle network handle
! @param idx local bus index
! @param pointer to buffer
!
  subroutine set_xc_bus_buffer(n_handle, idx, buf)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer(C_INT) c_handle, c_idx
    type(C_PTR) buf
    c_handle = n_handle
    c_idx = idx
    call p_set_xc_bus_buffer(c_handle, c_idx, buf)
    return
  end subroutine set_xc_bus_buffer
!
! Store the location of externally allocated buffer within a network
! @param n_handle network handle
! @param idx local branch index
! @param pointer to buffer
!
  subroutine set_xc_branch_buffer(n_handle, idx, buf)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer(C_INT) c_handle, c_idx
    type(C_PTR) buf
    c_handle = n_handle
    c_idx = idx
    call p_set_xc_branch_buffer(c_handle, c_idx, buf)
    return
  end subroutine set_xc_branch_buffer
!
! This subroutine must be called before calling the update bus routine. It
! initializes data structures for the bus update
! @param n_handle network handle
!
  subroutine init_bus_update(n_handle)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle
    integer(C_INT) c_handle
    c_handle = n_handle
    call p_init_bus_update(c_handle) 
    return
  end subroutine init_bus_update
!
! Update the ghost bus values
! @param n_handle network handle
!
  subroutine update_buses(n_handle)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle
    integer(C_INT) c_handle
    c_handle = n_handle
    call p_update_buses(c_handle) 
    return
  end subroutine update_buses
!
! This subroutine must be called before calling the update branch routine. It
! initializes data structures for the branch update
! @param n_handle network handle
!
  subroutine init_branch_update(n_handle)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle
    integer(C_INT) c_handle
    c_handle = n_handle
    call p_init_branch_update(c_handle) 
    return
  end subroutine init_branch_update
!
! Update the ghost branch values
! @param n_handle network handle
!
  subroutine update_branches(n_handle)
    use, intrinsic :: iso_c_binding
    use p_gridpack_network
    implicit none
    integer, value, intent(in) :: n_handle
    integer(C_INT) c_handle
    c_handle = n_handle
    call p_update_branches(c_handle) 
    return
  end subroutine update_branches
#endif
end module gridpack_network
