program main
  implicit none
  integer maxdim, nval
  parameter (maxdim = 100, nval = 100)
  integer i, j, one, stride
  double precision rval(nval)
  stride = 3
  one = 1
  write(6,'(a)') '# 1_2_3_4_5_6_7_8_9'
  write(6,'(a)') '# 1 2 3 4 5 6 7 8 9'
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
