! Wrappers for Fortran network functions
!
! Create new network and return handle 
!
integer(C_INT) function p_create_network(comm) bind(c, name="p_create_network")
  use gridpack_communicator
  use, intrinsic :: iso_c_binding
  implicit none
  type(C_PTR), value, intent(in) :: comm
end function p_create_network
!
! Add a bus locally to the network
!
subroutine p_add_bus(n_handle, idx) bind(c, name="p_add_bus")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
end subroutine p_add_bus
!
! Add a branch locally to the network
!
subroutine p_add_branch(n_handle, idx1, idx2) bind(c, name="p_add_branch")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx1, idx2
end subroutine p_add_branch
!
! Number of local buses (both active and inactive) on processor
!
integer(C_INT) function p_num_buses(n_handle) bind(c, name="p_num_buses")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle
end function p_num_buses
!
! Return the total number of buses in the entire network
!
integer(C_INT) function p_total_buses(n_handle) bind(c, name="p_total_buses")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle
end function p_total_buses
!
! Number of local branches (both active and inactive) on processor
!
integer(C_INT) function p_num_branches(n_handle) bind(c, name="p_num_branches")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle
end function p_num_branches
!
! Return the total number of buses in the entire network
!
integer(C_INT) function p_total_branches(n_handle) bind(c, name="p_total_branches")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle
end function p_total_branches
!
! Designate a bus as a reference bus
!
subroutine p_set_reference_bus(n_handle, idx) bind(c, name="p_set_reference_bus")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
end subroutine p_set_reference_bus
!
! Return index of reference bus. Return -1 if reference bus is not on this
! processor
!
integer(C_INT) function p_get_reference_bus(n_handle) &
   bind(c, name="p_get_reference_bus")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle
end function p_get_reference_bus
!
! Set the original index of the bus (from configuration file)
!
logical(C_BOOL) function p_set_original_bus_index(n_handle, idx, o_idx) &
   bind(c, name="p_set_original_bus_index")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx, o_idx
end function p_set_original_bus_index
!
! Set the global index of the bus
!
logical(C_BOOL) function p_set_global_bus_index(n_handle, idx, g_idx) &
   bind(c, name="p_set_global_bus_index")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx, g_idx
end function p_set_global_bus_index
!
! Set the global index of the branch
!
logical(C_BOOL) function p_set_global_branch_index(n_handle, idx, g_idx) &
   bind(c, name="p_set_global_branch_index")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx, g_idx
end function p_set_global_branch_index
!
! Set the original index of the bus at the "from" end of the branch
!
logical(C_BOOL) function p_set_original_bus_index1(n_handle, idx, b_idx) &
   bind(c, name="p_set_original_bus_index1")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx, b_idx
end function p_set_original_bus_index1
!
! Set the original index of the bus at the "to" end of the branch
!
logical(C_BOOL) function p_set_original_bus_index2(n_handle, idx, b_idx) &
   bind(c, name="p_set_original_bus_index2")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx, b_idx
end function p_set_original_bus_index2
!
! Set the global index of the bus at the "from" end of the branch
!
logical(C_BOOL) function p_set_global_bus_index1(n_handle, idx, b_idx) &
   bind(c, name="p_set_global_bus_index1")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx, b_idx
end function p_set_global_bus_index1
!
! Set the global index of the bus at the "to" end of the branch
!
logical(C_BOOL) function p_set_global_bus_index2(n_handle, idx, b_idx) &
   bind(c, name="p_set_global_bus_index2")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx, b_idx
end function p_set_global_bus_index2
!
! Set the local index of the bus at the "from" end of the branch
!
logical(C_BOOL) function p_set_local_bus_index1(n_handle, idx, b_idx) &
   bind(c, name="p_set_local_bus_index1")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx, b_idx
end function p_set_local_bus_index1
!
! Set the local index of the bus at the "to" end of the branch
!
logical(C_BOOL) function p_set_local_bus_index2(n_handle, idx, b_idx) &
   bind(c, name="p_set_local_bus_index2")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx, b_idx
end function p_set_local_bus_index2
!
! Set the active flag of the bus
!
logical(C_BOOL) function p_set_active_bus(n_handle, idx, flag) &
   bind(c, name="p_set_active_bus")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
  logical(C_BOOL), value, intent(in) :: flag
end function p_set_active_bus
!
! Set the active flag of the branch
!
logical(C_BOOL) function p_set_active_branch(n_handle, idx, flag) &
   bind(c, name="p_set_active_branch")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
  logical(C_BOOL), value, intent(in) :: flag
end function p_set_active_branch
!
! Clear the list of neighbors for the bus at idx
!
logical(C_BOOL) function p_clear_branch_neighbors(n_handle, idx) &
   bind(c, name="p_clear_branch_neighbors")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
end function p_clear_branch_neighbors
!
! Add a local index for a branch attached to the bus at idx
!
logical(C_BOOL) function p_add_branch_neighbor(n_handle, idx, br_idx) &
   bind(c, name="p_add_branch_neighbor")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx, br_idx
end function p_add_branch_neighbor
!
! Get status of the bus (local or ghosted)
!
logical(C_BOOL) function p_get_active_bus(n_handle, idx) &
   bind(c, name="p_get_active_bus")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
end function p_get_active_bus
!
! Get original index of the bus
!
integer(C_INT) function p_get_original_bus_index(n_handle, idx) &
   bind(c, name="p_get_original_bus_index")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
end function p_get_original_bus_index
!
! Get global index of the bus
!
integer(C_INT) function p_get_global_bus_index(n_handle, idx) &
   bind(c, name="p_get_global_bus_index")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
end function p_get_global_bus_index
!
! Get status of the branch (local or ghosted)
!
logical(C_BOOL) function p_get_active_branch(n_handle, idx) &
   bind(c, name="p_get_active_branch")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
end function p_get_active_branch
!
! Get the global index of the branch
!
integer(C_INT) function p_get_global_branch_index(n_handle, idx) &
   bind(c, name="p_get_global_branch_index")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
end function p_get_global_branch_index
!
! Get original indices of the two buses at each end of the branch
!
subroutine p_get_original_branch_endpoints(n_handle, idx, idx1, idx2) &
   bind(c, name="p_get_original_branch_endpoints")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
  integer(C_INT), intent(out) :: idx1, idx2
end subroutine p_get_original_branch_endpoints
!
! Return the number of branches connected to a bus. This can be used to allocate
! arrays to hold branch indices
!
integer(C_INT) function p_get_num_connected_branches(n_handle, idx) &
   bind(c, name="p_get_num_connected_branches")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
end function p_get_num_connected_branches
!
! Return list of branches connected to bus
!
subroutine p_get_connected_branches(n_handle, idx, branches) &
   bind(c, name="p_get_connected_branches")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
  integer(C_INT), intent(out) :: branches(*)
end subroutine p_get_connected_branches
!
! Return list of buses connected to bus via one branch
!
subroutine p_get_connected_buses(n_handle, idx, buses) &
   bind(c, name="p_get_connected_buses")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
  integer(C_INT), intent(out) :: buses(*)
end subroutine p_get_connected_buses
!
! Return the indices of the buses at either end of the branch
!
subroutine p_get_branch_endpoints(n_handle, idx, idx1, idx2) &
   bind(c, name="p_get_branch_endpoints")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
  integer(C_INT), intent(out) :: idx1, idx2
end subroutine p_get_branch_endpoints
!
! Partition the network over the available processes
!
subroutine p_partition(n_handle) bind(c, name="p_partition")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle
end subroutine p_partition
!
! Clean all ghost buses and branches from the system. This can be used before
! repartitioning the network
!
subroutine p_clean(n_handle) bind(c, name="p_clean")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle
end subroutine p_clean
!
! Store the location of externally allocated buffer within a network
!
subroutine p_set_xc_bus_buffer(n_handle, idx, buf) &
   bind(c, name="p_set_xc_bus_buffer")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
  type(C_PTR), intent(in) :: buf
end subroutine p_set_xc_bus_buffer
!
! Store the location of externally allocated buffer within a network
!
subroutine p_set_xc_branch_buffer(n_handle, idx, buf) &
   bind(c, name="p_set_xc_branch_buffer")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle, idx
  type(C_PTR), intent(in) :: buf
end subroutine p_set_xc_branch_buffer
!
! This subroutine must be called before calling the update bus routine. It
! initializes data structures for the bus update
!
subroutine p_init_bus_update(n_handle) bind(c, name="p_init_bus_update")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle
end subroutine p_init_bus_update
!
! Update the bus ghost values.
!
subroutine p_update_buses(n_handle) bind(c, name="p_update_buses")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle
end subroutine p_update_buses
!
! This subroutine must be called before calling the update branch routine. It
! initializes data structures for the branch update
!
subroutine p_init_branch_update(n_handle) bind(c, name="p_init_branch_update")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle
end subroutine p_init_branch_update
!
! Update the branch ghost values.
!
subroutine p_update_branches(n_handle) bind(c, name="p_update_branches")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value, intent(in) :: n_handle
end subroutine p_update_branches
