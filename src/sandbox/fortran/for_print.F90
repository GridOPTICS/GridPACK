subroutine for_print(i,v) bind(c)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(C_INT), value :: i ! i is being passed into routine, add value attribute
  real(C_DOUBLE), intent(out) :: v ! v is being passed out of routine
  v = dble(i)
  write(6,'(a,i4,a,f6.1)') 'From Fortran: integer: ',i,' double: ',v
  return
end subroutine for_print
