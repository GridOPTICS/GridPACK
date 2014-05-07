!
!  Declare data type in a module so that it can be used across multiple code
!  units
!
module data_type
  use, intrinsic :: iso_c_binding
  implicit none
!
!  Derived type that is used to construct the bus exchange buffer. Note that the
!  contents of this buffers are defined by the application, although they should
!  follow the general formatting rules for constructing compatible types using
!  iso_c_binding. The name field is the case sensitive name that can be used to
!  identify this data type inside C++ code. It should not be changed.
!
  type, bind(c) :: busxcbuf
   integer(C_INT) :: i1, i2
   real(C_DOUBLE) :: d1, d2
   complex(C_DOUBLE_COMPLEX) :: c1, c2
   character(kind=C_CHAR,len=32) :: str
  end type
end module data_type
!
program f_test
  use, intrinsic :: iso_c_binding
  use data_type
  implicit none
!
!  interface declaration for subroutine overwrite
!
  interface
    subroutine overwrite(bsize,a,b)
      use, intrinsic :: iso_c_binding
      use data_type
      implicit none
      type(busxcbuf) :: a, b
      integer bsize
    end subroutine
  end interface
!
!  declare a vector containing ten elements of type busxcbuf
!
  type(busxcbuf) data(10)
  integer i, bsize
!
!  Initialize data array
!
  do i = 1, 10
    data(i)%i1 = i
    data(i)%i2 = i + 1
    data(i)%d1 = dble(i)
    data(i)%d2 = dble(i)+1.0d00
    data(i)%c1 = dcmplx(dble(i),dble(i+1))
    data(i)%c2 = dcmplx(dble(i+2),dble(i+3))
    write(data(i)%str,'(a,i2)') 'Storing value on element ',i
  end do
!
!   Copy data from array element 2 to array element 6 using a function written
!   C.
!
  bsize = 2*4 + 2*8 + 2*16 + 32
  call overwrite(bsize,data(2),data(6));
!
!   Print out values of data elements 2 and 6
!
  write(6,'(a,2i4)') 'Element 2 integers: ',data(2)%i1,data(2)%i2
  write(6,'(a,2f6.1)') 'Element 2 doubles: ',data(2)%d1,data(2)%d2
  write(6,'(a,f6.1,a,f6.1,a,f6.2,a,f6.1,a)') 'Element 2 double complex: (', &
    real(data(2)%c1),',',imag(data(2)%c1),'), (', &
    real(data(2)%c2),',',imag(data(2)%c2),')'
  write(6,'(a,a)') 'Element 2 string: ',data(2)%str
  write(6,*)
  write(6,'(a,2i4)') 'Element 6 integers: ',data(6)%i1,data(6)%i2
  write(6,'(a,2f6.1)') 'Element 6 doubles: ',data(6)%d1,data(6)%d2
  write(6,'(a,f6.1,a,f6.1,a,f6.2,a,f6.1,a)') 'Element 6 double complex: (', &
    real(data(6)%c1),',',imag(data(6)%c1),'), (', &
    real(data(6)%c2),',',imag(data(6)%c2),')'
  write(6,'(a,a)') 'Element 6 string: ',data(6)%str
  stop
end program f_test
subroutine overwrite(bsize,a,b)
  use, intrinsic :: iso_c_binding
  use data_type
  implicit none
  integer bsize
  type(busxcbuf), target :: a, b
  type(busxcbuf), pointer :: a_fptr, b_fptr
  integer(C_INT) :: esize
  type(C_PTR) :: a_ptr, b_ptr
  type(busxcbuf) t(2)
!
!  Interface declaration that guarantees that Fortran function matches up with
!  implementation in C/C++. The subroutine c_overwrite is the actual implementation
!  of this subroutine in C.
!
  interface
    subroutine c_overwrite(esize,a_ptr,b_ptr) bind(c, name="c_overwrite")
      use, intrinsic :: iso_c_binding
      use data_type
      implicit none
      integer(C_INT) :: esize
      type(C_PTR) :: a_ptr, b_ptr
    end subroutine c_overwrite
  end interface
!
!   Copy data from a to b using a function written C. The first argument is
!   the size of an array element in bytes
!
  a_fptr => a
  b_fptr => b
  esize = bsize
  a_ptr = c_loc(a) 
  b_ptr = c_loc(b) 
  call c_overwrite(esize, a_ptr, b_ptr)
  call c_f_pointer(a_ptr, a_fptr)
  call c_f_pointer(b_ptr, b_fptr)
  return
end subroutine overwrite
