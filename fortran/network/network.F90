! Fortran network functions
!
! Create new network and return handle 
!
integer function create_network(comm)
  use gridpack_communicator
  use iso_c_binding
  implicit none
  type(C_PTR), value, intent(in) :: comm
  integer(C_INT) p_create_network
  integer(C_INT) handle
  handle = p_create_network(comm)
  create_network = handle
  return
end function create_network
