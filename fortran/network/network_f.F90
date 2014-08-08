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
    type (C_PTR) :: p_network
    contains
    procedure::create
    procedure::destroy
    procedure::add_bus
    procedure::add_branch
    procedure::num_buses
    procedure::total_buses
    procedure::num_branches
    procedure::total_branches
    procedure::set_reference_bus
    procedure::get_reference_bus
    procedure::set_original_bus_index
    procedure::set_global_bus_index
    procedure::set_global_branch_index
    procedure::set_original_bus_index1
    procedure::set_original_bus_index2
    procedure::set_global_bus_index1
    procedure::set_global_bus_index2
    procedure::set_local_bus_index1
    procedure::set_local_bus_index2
    procedure::set_active_bus
    procedure::set_active_branch
    procedure::clear_branch_neighbors
    procedure::add_branch_neighbor
    procedure::get_active_bus
    procedure::get_original_bus_index
    procedure::get_global_bus_index
    procedure::get_active_branch
    procedure::get_global_branch_index
    procedure::get_original_branch_endpoints
    procedure::get_num_connected_branches
    procedure::get_connected_branches
    procedure::get_connected_buses
    procedure::get_branch_endpoints
    procedure::partition
    procedure::clean
    procedure::init_bus_update
    procedure::update_buses
    procedure::init_branch_update
    procedure::update_branches
  end type
  interface
!
! Create new network and return handle 
! @param network new GridPACK network object
! @param comm GridPACK communicator
!
    subroutine network_create(network,comm) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: network
      type(C_PTR), intent(in) :: comm
    end subroutine network_create
!
! Clean up old network 
! @param network old GridPACK network object
! @return handle for new network
!
    subroutine network_destroy(network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
    end subroutine network_destroy
!
! Add a bus locally to the network
! @param network  GridPACK fortran network
! @param idx original index of bus
!
    subroutine network_add_bus(network,idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), intent(in) :: idx
    end subroutine network_add_bus
!
! Add a branch locally to the network
! @param network  GridPACK fortran network
! @param idx1 original index of "from" bus
! @param idx2 original index of "to" bus
!
    subroutine network_add_branch(network,idx1,idx2) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), intent(in) :: idx1, idx2
    end subroutine network_add_branch
!
! Number of local buses (both active and inactive) on processor
! @param network GridPACK network object
! @return number of local buses
! 
    integer(C_INT) function network_num_buses(network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
    end function network_num_buses
!
! Return the total number of buses in the entire network
! @param network GridPACK network object
! @return total number of buses in network
! 
    integer(C_INT) function network_total_buses(network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
    end function network_total_buses
!
! Number of local branches (both active and inactive) on processor
! @param network GridPACK network object
! @return number of local branches
! 
    integer(C_INT) function network_num_branches(network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
    end function network_num_branches
!
! Return the total number of branches in the entire network
! @param network GridPACK network object
! @return total number of branches in network
! 
    integer(C_INT) function network_total_branches(network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
    end function network_total_branches
!
! Designate a bus as a reference bus
! @param network GridPACK network object
! @param idx local bus index
!
    subroutine network_set_reference_bus(network, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
    end subroutine network_set_reference_bus
!
! Return the index of reference bus. Return -1 if reference bus is not on this
! processor
! @param network GridPACK network object
! @return local index of reference bus
!
    integer(C_INT) function network_get_reference_bus(network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
    end function network_get_reference_bus
!
! Set the original index of the bus (from configuration file)
! @param network GridPACK network object
! @param idx local bus index
! @param o_idx original bus index
! @return false if bus does not exist at that value of idx
!
    logical(C_BOOL) function network_set_original_bus_index(network, &
        idx, o_idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx, o_idx
    end function network_set_original_bus_index
!
! Set the global index of the bus
! @param network GridPACK network object
! @param idx local bus index
! @param o_idx global bus index
! @return false if bus does not exist at that value of idx
!
    logical(C_BOOL) function network_set_global_bus_index(network, &
        idx, g_idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx, g_idx
    end function network_set_global_bus_index
!
! Set the global index of the branch
! @param network GridPACK network object
! @param idx local branch index
! @param o_idx global branch index
! @return false if branch does not exist at that value of idx
!
    logical(C_BOOL) function network_set_global_branch_index(network, &
        idx, g_idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx, g_idx
    end function network_set_global_branch_index
!
! Set the origina1 index of the bus at the "from" end of branch
! @param network GridPACK network object
! @param idx local branch index
! @param b_idx original index of "from" bus
! @return false if branch does not exist at that value of idx
!
    logical(C_BOOL) function network_set_original_bus_index1(network, &
        idx, b_idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx, b_idx
    end function network_set_original_bus_index1
!
! Set the origina1 index of the bus at the "to" end of branch
! @param network GridPACK network object
! @param idx local branch index
! @param b_idx original index of "to" bus
! @return false if branch does not exist at that value of idx
!
    logical(C_BOOL) function network_set_original_bus_index2(network, &
        idx, b_idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx, b_idx
    end function network_set_original_bus_index2
!
! Set the globa1 index of the bus at the "from" end of branch
! @param network GridPACK network object
! @param idx local branch index
! @param b_idx global index of "from" bus
! @return false if branch does not exist at that value of idx
!
    logical(C_BOOL) function network_set_global_bus_index1(network, &
        idx, b_idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx, b_idx
    end function network_set_global_bus_index1
!
! Set the globa1 index of the bus at the "to" end of branch
! @param network GridPACK network object
! @param idx local branch index
! @param b_idx global index of "to" bus
! @return false if branch does not exist at that value of idx
!
    logical(C_BOOL) function network_set_global_bus_index2(network, &
        idx, b_idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx, b_idx
    end function network_set_global_bus_index2
!
! Set the loca1 index of the bus at the "from" end of branch
! @param network GridPACK network object
! @param idx local branch index
! @param b_idx local index of "from" bus
! @return false if branch does not exist at that value of idx
!
    logical(C_BOOL) function network_set_local_bus_index1(network, &
        idx, b_idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx, b_idx
    end function network_set_local_bus_index1
!
! Set the loca1 index of the bus at the "to" end of branch
! @param network GridPACK network object
! @param idx local branch index
! @param b_idx local index of "to" bus
! @return false if branch does not exist at that value of idx
!
    logical(C_BOOL) function network_set_local_bus_index2(network, &
        idx, b_idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx, b_idx
    end function network_set_local_bus_index2
!
! Set the active flag of the bus
! @param network GridPACK network object
! @param idx local bus index
! @param flag active status of bus
! @return false if bus does not exist at that value of idx
! 
    logical(C_BOOL) function network_set_active_bus(network, idx, flag) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
      logical(C_BOOL), value, intent(in) :: flag
    end function network_set_active_bus
!
! Set the active flag of the branch
! @param network GridPACK network object
! @param idx local branch index
! @param flag active status of branch
! @return false if branch does not exist at that value of idx
! 
    logical(C_BOOL) function network_set_active_branch(network, idx, flag) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
      logical(C_BOOL), value, intent(in) :: flag
    end function network_set_active_branch
!
! Clear the list of neighbors for the bus at idx
! @param p_network GridPACK network object
! @param idx local branch index
! @return false if bus does not exist at that value of idx
! 
    logical(C_BOOL) function network_clear_branch_neighbors(network, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
    end function network_clear_branch_neighbors
!
! Add a local index for a branch attached to the bus at idx
! @param network GridPACK network object
! @param idx local bus index
! @param br_idx local branch index
! @return false if bus does not exist at that value of idx
! 
    logical(C_BOOL) function network_add_branch_neighbor(network, &
        idx, br_idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx, br_idx
    end function network_add_branch_neighbor
!
! Get status of the bus (local or ghosted)
! @param network GridPACK network object
! @param idx local bus index
! @return status of bus
!
    logical(C_BOOL) function network_get_active_bus(network, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
    end function network_get_active_bus
!
! Get original index of bus
! @param network GridPACK network object
! @param idx local bus index
! @return original index of bus
!
    integer(C_INT) function network_get_original_bus_index(network, &
        idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
    end function network_get_original_bus_index
!
! Get global index of bus
! @param network GridPACK network object
! @param idx local bus index
! @return global index of bus
!
    integer(C_INT) function network_get_global_bus_index(network, &
        idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
    end function network_get_global_bus_index
!
! Get status of the branch (local or ghosted)
! @param network GridPACK network object
! @param idx local branch index
! @return status of branch
!
    logical(C_BOOL) function network_get_active_branch(network, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
    end function network_get_active_branch
!
! Get global index of branch
! @param network GridPACK network object
! @param idx local branch index
! @return global index of bus
!
    integer(C_INT) function network_get_global_branch_index(network, idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
    end function network_get_global_branch_index
!
! Get original indices of buses at each end of the branch
! @param network GridPACK network object
! @param idx local branch index
! @param idx1 original index of "from" bus
! @param idx2 original index of "to" bus
!
    subroutine network_get_original_branch_endpoints(network, &
        idx, idx1, idx2) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
      integer(C_INT), intent(out) :: idx1, idx2
    end subroutine network_get_original_branch_endpoints
!
! Return the number of branches connected to a bus. This can be used to allocate
! arrays to hold the branch indices
! @param network GridPACK network object
! @param idx local bus index
! @return number of connected branches
!
    integer(C_INT) function network_get_num_connected_branches(network, &
        idx) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
    end function network_get_num_connected_branches
!
! Return list of branches connected to bus
! @param network GridPACK network object
! @param idx local bus index
! @param branches list of local branch indices
!
    subroutine network_get_connected_branches(network, idx, branches) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
      integer(C_INT), intent(out) :: branches(*)
    end subroutine network_get_connected_branches
!
! Return list of buses connected to bus via one branch
! @param network GridPACK network object
! @param idx local bus index
! @param buses list of local bus indices
!
    subroutine network_get_connected_buses(network, idx, buses) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
      integer(C_INT), intent(out) :: buses(*)
    end subroutine network_get_connected_buses
!
! Return the indices of the buses at either end of the branch
! @param network GridPACK network object
! @param idx local branch index
! @param idx1 local index of "from" bus
! @param idx2 local index of "to" bus
!
    subroutine network_get_branch_endpoints(network, idx, idx1, idx2) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
      integer(C_INT), value, intent(in) :: idx
      integer(C_INT), intent(out) :: idx1, idx2
    end subroutine network_get_branch_endpoints
!
! Partition network over available processes
! @param network GridPACK network object
!
    subroutine network_partition(network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
    end subroutine network_partition
!
! Clean all ghost buses and branches from the system. This can be used before
! repartioning the network
! @param network GridPACK network object
!
    subroutine network_clean(network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
    end subroutine network_clean
!
! This subroutine must be called before calling the update bus routine. It
! initializes data structures for the bus update
! @param network GridPACK network object
!
    subroutine network_init_bus_update(network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
    end subroutine network_init_bus_update
!
! Update the ghost bus values
! @param network GridPACK network object
!
    subroutine network_update_buses(network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
    end subroutine network_update_buses
!
! This subroutine must be called before calling the update branch routine. It
! initializes data structures for the branch update
! @param network GridPACK network object
!
    subroutine network_init_branch_update(network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
    end subroutine network_init_branch_update
!
! Update the ghost branch values
! @param network GridPACK network object
!
    subroutine network_update_branches(network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: network
    end subroutine network_update_branches
  end interface
  contains
!
! Create new network and return handle 
! @param p_network new GridPACK network object
! @param comm GridPACK communicator
!
  subroutine create(p_network,p_comm)
    use gridpack_communicator
    implicit none
    class(network), intent(inout) :: p_network
    class(communicator), value, intent(in) :: p_comm
    call network_create(p_network%p_network,p_comm%comm)
    return
  end subroutine create
!
! Clean up old network
! @param p_network new GridPACK network object
! @param comm GridPACK communicator
!
  subroutine destroy(p_network)
    implicit none
    class(network), intent(inout) :: p_network
    call network_destroy(p_network%p_network)
    return
  end subroutine destroy
!
! Add a bus locally to the network
! @param p_network GridPACK network object
! @param idx original index of bus
!
  subroutine add_bus(p_network,idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer(C_INT) c_idx
    c_idx = idx
    call network_add_bus(p_network%p_network,c_idx)
    return
  end subroutine add_bus
!
! Add a branch locally to the network
! @param p_network GridPACK network object
! @param idx1 original index of "from" bus
! @param idx2 original index of "to" bus
!
  subroutine add_branch(p_network, idx1, idx2)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx1, idx2
    integer(C_INT) c_idx1, c_idx2
    c_idx1 = idx1
    c_idx2 = idx2
    call network_add_branch(p_network%p_network, c_idx1, c_idx2)
    return
  end subroutine add_branch
!
! Number of local buses (both active and inactive) on processor
! @param p_network GridPACK network object
! @return number of local buses
! 
  integer function num_buses(p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer(C_INT) nbus
    nbus = network_num_buses(p_network%p_network)
    num_buses = nbus
    return
  end function num_buses
!
! Return the total number of buses in the entire network
! @param p_network GridPACK network object
! @return total number of buses in network
! 
  integer function total_buses(p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer(C_INT) nbus
    nbus = network_total_buses(p_network%p_network)
    total_buses = nbus
    return
  end function total_buses
!
! Number of local branches (both active and inactive) on processor
! @param p_network GridPACK network object
! @return number of local branches
! 
  integer function num_branches(p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer(C_INT) nbranch
    nbranch = network_num_branches(p_network%p_network)
    num_branches = nbranch
    return
  end function num_branches
!
! Return the total number of branches in the entire network
! @param p_network GridPACK network object
! @return total number of branches in network
! 
  integer function total_branches(p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer(C_INT) nbranch
    nbranch = network_total_branches(p_network%p_network)
    total_branches = nbranch
    return
  end function total_branches
!
! Designate a bus as a reference bus
! @param p_network GridPACK network object
! @param idx local bus index
!
  subroutine set_reference_bus(p_network, idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer(C_INT) c_idx
    c_idx = idx
    call network_set_reference_bus(p_network%p_network, c_idx)
    return
  end subroutine set_reference_bus
!
! Return the index of reference bus. Return -1 if reference bus is not on this
! processor
! @param p_network GridPACK network object
! @return local index of reference bus
!
  integer function get_reference_bus(p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer(C_INT) c_idx
    c_idx = network_get_reference_bus(p_network%p_network)
    get_reference_bus = c_idx
    return
  end function get_reference_bus
!
! Set the original index of the bus (from configuration file)
! @param p_network GridPACK network object
! @param idx local bus index
! @param o_idx original bus index
! @return false if bus does not exist at that value of idx
!
  logical function set_original_bus_index(p_network, idx, o_idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx, o_idx
    integer(C_INT) c_idx, c_o_idx
    logical(C_BOOL) c_flag
    c_idx = idx
    c_o_idx = o_idx
    c_flag = network_set_original_bus_index(p_network%p_network, c_idx, c_o_idx)
    set_original_bus_index = c_flag
    return
  end function set_original_bus_index
!
! Set the global index of the bus
! @param p_network GridPACK network object
! @param idx local bus index
! @param o_idx global bus index
! @return false if bus does not exist at that value of idx
!
  logical function set_global_bus_index(p_network, idx, g_idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx, g_idx
    integer(C_INT) c_idx, c_g_idx
    logical(C_BOOL) c_flag
    c_idx = idx
    c_g_idx = g_idx
    c_flag = network_set_global_bus_index(p_network%p_network, c_idx, c_g_idx)
    set_global_bus_index = c_flag
    return
  end function set_global_bus_index
!
! Set the global index of the branch
! @param p_network GridPACK network object
! @param idx local branch index
! @param o_idx global branch index
! @return false if branch does not exist at that value of idx
!
  logical function set_global_branch_index(p_network, idx, g_idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx, g_idx
    integer(C_INT) c_idx, c_g_idx
    logical(C_BOOL) c_flag
    c_idx = idx
    c_g_idx = g_idx
    c_flag = network_set_global_branch_index(p_network%p_network, c_idx, c_g_idx)
    set_global_branch_index = c_flag
    return
  end function set_global_branch_index
!
! Set the origina1 index of the bus at the "from" end of branch
! @param p_network GridPACK network object
! @param idx local branch index
! @param b_idx original index of "from" bus
! @return false if branch does not exist at that value of idx
!
  logical function set_original_bus_index1(p_network, idx, b_idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx, b_idx
    integer(C_INT) c_idx, c_b_idx
    logical(C_BOOL) c_flag
    c_idx = idx
    c_b_idx = b_idx
    c_flag = network_set_original_bus_index1(p_network%p_network, c_idx, c_b_idx)
    set_original_bus_index1 = c_flag
    return
  end function set_original_bus_index1
!
! Set the origina1 index of the bus at the "to" end of branch
! @param p_network GridPACK network object
! @param idx local branch index
! @param b_idx original index of "to" bus
! @return false if branch does not exist at that value of idx
!
  logical function set_original_bus_index2(p_network, idx, b_idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx, b_idx
    integer(C_INT) c_idx, c_b_idx
    logical(C_BOOL) c_flag
    c_idx = idx
    c_b_idx = b_idx
    c_flag = network_set_original_bus_index2(p_network%p_network, c_idx, c_b_idx)
    set_original_bus_index2 = c_flag
    return
  end function set_original_bus_index2
!
! Set the globa1 index of the bus at the "from" end of branch
! @param p_network GridPACK network object
! @param idx local branch index
! @param b_idx global index of "from" bus
! @return false if branch does not exist at that value of idx
!
  logical function set_global_bus_index1(p_network, idx, b_idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx, b_idx
    integer(C_INT) c_idx, c_b_idx
    logical(C_BOOL) c_flag
    c_idx = idx
    c_b_idx = b_idx
    c_flag = network_set_global_bus_index1(p_network%p_network, c_idx, c_b_idx)
    set_global_bus_index1 = c_flag
    return
  end function set_global_bus_index1
!
! Set the globa1 index of the bus at the "to" end of branch
! @param p_network GridPACK network object
! @param idx local branch index
! @param b_idx global index of "to" bus
! @return false if branch does not exist at that value of idx
!
  logical function set_global_bus_index2(p_network, idx, b_idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx, b_idx
    integer(C_INT) c_idx, c_b_idx
    logical(C_BOOL) c_flag
    c_idx = idx
    c_b_idx = b_idx
    c_flag = network_set_global_bus_index2(p_network%p_network, c_idx, c_b_idx)
    set_global_bus_index2 = c_flag
    return
  end function set_global_bus_index2
!
! Set the loca1 index of the bus at the "from" end of branch
! @param p_network GridPACK network object
! @param idx local branch index
! @param b_idx local index of "from" bus
! @return false if branch does not exist at that value of idx
!
  logical function set_local_bus_index1(p_network, idx, b_idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx, b_idx
    integer(C_INT) c_idx, c_b_idx
    logical(C_BOOL) c_flag
    c_idx = idx
    c_b_idx = b_idx
    c_flag = network_set_local_bus_index1(p_network%p_network, c_idx, c_b_idx)
    set_local_bus_index1 = c_flag
    return
  end function set_local_bus_index1
!
! Set the loca1 index of the bus at the "to" end of branch
! @param p_network GridPACK network object
! @param idx local branch index
! @param b_idx local index of "to" bus
! @return false if branch does not exist at that value of idx
!
  logical function set_local_bus_index2(p_network, idx, b_idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx, b_idx
    integer(C_INT) c_idx, c_b_idx
    logical(C_BOOL) c_flag
    c_idx = idx
    c_b_idx = b_idx
    c_flag = network_set_local_bus_index2(p_network%p_network, c_idx, c_b_idx)
    set_local_bus_index2 = c_flag
    return
  end function set_local_bus_index2
!
! Set the active flag of the bus
! @param p_network GridPACK network object
! @param idx local bus index
! @param flag active status of bus
! @return false if bus does not exist at that value of idx
! 
  logical function set_active_bus(p_network, idx, flag)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    logical, value, intent(in) :: flag
    integer(C_INT) c_idx
    logical(C_BOOL) c_flag, c_ret
    c_idx = idx
    c_flag = flag
    write(6,'(a)') 'set_active_bus'
    c_ret = network_set_active_bus(p_network%p_network, c_idx, c_flag)
    write(6,'(a)') 'completed network_set_active_bus'
    set_active_bus = c_ret
    return
  end function set_active_bus
!
! Set the active flag of the branch
! @param p_network GridPACK network object
! @param idx local branch index
! @param flag active status of branch
! @return false if branch does not exist at that value of idx
! 
  logical function set_active_branch(p_network, idx, flag)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    logical flag
    integer(C_INT) c_idx
    logical(C_BOOL) c_flag, c_ret
    c_idx = idx
    c_flag = flag
    c_ret = network_set_active_branch(p_network%p_network, c_idx, c_flag)
    set_active_branch = c_ret
    return
  end function set_active_branch
!
! Clear the list of neighbors for the bus at idx
! @param p_network GridPACK network object
! @param idx local branch index
! @return false if bus does not exist at that value of idx
! 
  logical function clear_branch_neighbors(p_network, idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer(C_INT) c_idx
    logical(C_BOOL) c_ret
    c_idx = idx
    c_ret = network_clear_branch_neighbors(p_network%p_network, c_idx)
    clear_branch_neighbors = c_ret
    return
  end function clear_branch_neighbors
!
! Add a local index for a branch attached to the bus at idx
! @param p_network GridPACK network object
! @param idx local bus index
! @param br_idx local branch index
! @return false if bus does not exist at that value of idx
! 
  logical function add_branch_neighbor(p_network, idx, br_idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx, br_idx
    integer(C_INT) c_idx, c_br_idx
    logical(C_BOOL) c_ret
    c_idx = idx
    c_br_idx = br_idx
    c_ret = network_add_branch_neighbor(p_network%p_network, c_idx, c_br_idx)
    add_branch_neighbor = c_ret
    return
  end function add_branch_neighbor
!
! Get status of the bus (local or ghosted)
! @param p_network GridPACK network object
! @param idx local bus index
! @return status of bus
!
  logical function get_active_bus(p_network, idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer(C_INT) c_idx
    logical(C_BOOL) c_ret
    c_idx = idx
    c_ret = network_get_active_bus(p_network%p_network, c_idx)
    get_active_bus = c_ret
    return
  end function get_active_bus
!
! Get original index of bus
! @param p_network GridPACK network object
! @param idx local bus index
! @return original index of bus
!
  integer function get_original_bus_index(p_network, idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer(C_INT) c_idx
    integer(C_INT) c_ret
    c_idx = idx
    c_ret = network_get_original_bus_index(p_network%p_network, c_idx)
    get_original_bus_index = c_ret
    return
  end function get_original_bus_index
!
! Get global index of bus
! @param p_network GridPACK network object
! @param idx local bus index
! @return global index of bus
!
  integer function get_global_bus_index(p_network, idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer(C_INT) c_idx
    integer(C_INT) c_ret
    c_idx = idx
    c_ret = network_get_global_bus_index(p_network%p_network, c_idx)
    get_global_bus_index = c_ret
    return
  end function get_global_bus_index
!
! Get status of the branch (local or ghosted)
! @param p_network GridPACK network object
! @param idx local branch index
! @return status of branch
!
  logical function get_active_branch(p_network, idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer(C_INT) c_idx
    logical(C_BOOL) c_ret
    c_idx = idx
    c_ret = network_get_active_branch(p_network%p_network, c_idx)
    get_active_branch = c_ret
    return
  end function get_active_branch
!
! Get global index of branch
! @param p_network GridPACK network object
! @param idx local branch index
! @return global index of bus
!
  integer function get_global_branch_index(p_network, idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer(C_INT) c_idx
    integer(C_INT) c_ret
    c_idx = idx
    c_ret = network_get_global_branch_index(p_network%p_network, c_idx)
    get_global_branch_index = c_ret
    return
  end function get_global_branch_index
!
! Get original indices of buses at each end of the branch
! @param p_network GridPACK network object
! @param idx local branch index
! @param idx1 original index of "from" bus
! @param idx2 original index of "to" bus
!
  subroutine get_original_branch_endpoints(p_network, idx, idx1, idx2)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer, intent(out) :: idx1, idx2
    integer(C_INT) c_idx
    integer(C_INT) c_idx1, c_idx2
    c_idx = idx
    call network_get_original_branch_endpoints(p_network%p_network, &
      c_idx, c_idx1, c_idx2)
    idx1 = c_idx1
    idx2 = c_idx2
    return
  end subroutine get_original_branch_endpoints
!
! Return the number of branches connected to a bus. This can be used to allocate
! arrays to hold the branch indices
! @param p_network GridPACK network object
! @param idx local bus index
! @return number of connected branches
!
  integer function get_num_connected_branches(p_network, idx)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer(C_INT) c_idx
    integer(C_INT) c_ret
    c_idx = idx
    c_ret = network_get_num_connected_branches(p_network%p_network, c_idx)
    get_num_connected_branches = c_ret
    return
  end function get_num_connected_branches
!
! Return list of branches connected to bus
! @param p_network GridPACK network object
! @param idx local bus index
! @param branches list of local branch indices
!
  subroutine get_connected_branches(p_network, idx, branches)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer, intent(out) :: branches(*)
    integer(C_INT) c_idx, c_num
    integer(C_INT), allocatable :: c_branches(:)
    integer i
    c_idx = idx
    c_num = network_get_num_connected_branches(p_network%p_network, c_idx)
    allocate(c_branches(c_num))
    call network_get_connected_branches(p_network%p_network, c_idx, c_branches)
    do i = 1, c_num
      branches(i) = c_branches(i)
    end do
    deallocate(c_branches)
    return
  end subroutine get_connected_branches
!
! Return list of buses connected to bus via one branch
! @param p_network GridPACK network object
! @param idx local bus index
! @param buses list of local bus indices
!
  subroutine get_connected_buses(p_network, idx, buses)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer, intent(out) :: buses(*)
    integer(C_INT) c_idx, c_num
    integer(C_INT), allocatable :: c_buses(:)
    integer i
    c_idx = idx
    c_num = network_get_num_connected_branches(p_network%p_network, c_idx)
    allocate(c_buses(c_num))
    call network_get_connected_buses(p_network%p_network, c_idx, c_buses)
    do i = 1, c_num
      buses(i) = c_buses(i)
    end do
    deallocate(c_buses)
    return
  end subroutine get_connected_buses
!
! Return the indices of the buses at either end of the branch
! @param p_network GridPACK network object
! @param idx local branch index
! @param idx1 local index of "from" bus
! @param idx2 local index of "to" bus
!
  subroutine get_branch_endpoints(p_network, idx, idx1, idx2)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    integer, value, intent(in) :: idx
    integer, intent(out) :: idx1, idx2
    integer(C_INT) c_idx, c_idx1, c_idx2
    c_idx = idx
    call network_get_branch_endpoints(p_network%p_network, c_idx, c_idx1, c_idx2)
    idx1 = c_idx1
    idx2 = c_idx2
    return
  end subroutine get_branch_endpoints
!
! Partition network over available processes
! @param p_network GridPACK network object
!
  subroutine partition(p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    call network_partition(p_network%p_network)
    return
  end subroutine partition
!
! Clean all ghost buses and branches from the system. This can be used before
! repartioning the network
! @param p_network GridPACK network object
!
  subroutine clean(p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    call network_clean(p_network%p_network)
    return
  end subroutine clean
#if 0
!
! Store the location of externally allocated buffer within a network
! @param n_handle network handle
! @param idx local bus index
! @param pointer to buffer
!
  subroutine set_xc_bus_buffer(n_handle, idx, buf)
    use, intrinsic :: iso_c_binding
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
    implicit none
    integer, value, intent(in) :: n_handle, idx
    integer(C_INT) c_handle, c_idx
    type(C_PTR) buf
    c_handle = n_handle
    c_idx = idx
    call p_set_xc_branch_buffer(c_handle, c_idx, buf)
    return
  end subroutine set_xc_branch_buffer
#endif
!
! This subroutine must be called before calling the update bus routine. It
! initializes data structures for the bus update
! @param p_network GridPACK network object
!
  subroutine init_bus_update(p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    call network_init_bus_update(p_network%p_network) 
    return
  end subroutine init_bus_update
!
! Update the ghost bus values
! @param p_network GridPACK network object
!
  subroutine update_buses(p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    call network_update_buses(p_network%p_network) 
    return
  end subroutine update_buses
!
! This subroutine must be called before calling the update branch routine. It
! initializes data structures for the branch update
! @param p_network GridPACK network object
!
  subroutine init_branch_update(p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    call network_init_branch_update(p_network%p_network) 
    return
  end subroutine init_branch_update
!
! Update the ghost branch values
! @param p_network GridPACK network object
!
  subroutine update_branches(p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(network), intent(in) :: p_network
    call network_update_branches(p_network%p_network) 
    return
  end subroutine update_branches
end module gridpack_network
