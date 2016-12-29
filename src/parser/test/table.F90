program main
  implicit none
  integer maxdim, nval
  parameter (maxdim = 100, nval = 100)
  integer i, j, one, stride
  double precision rval(nval)
  stride = 3
  one = 1
  do i=1, maxdim
    if (mod(i,stride).eq.0) then
      do j = 1, nval
        rval(j) = dble((i-1)*nval+j-1)
      end do 
      write(6,'(i4,i2,100f12.1)')i,one,(rval(j),j=1,nval)
    endif
  end do
  stop
end program
