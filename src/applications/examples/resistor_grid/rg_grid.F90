program main
  implicit none
  integer maxdim
  parameter (maxdim = 10)
  integer i, j, n, nbus, icnt,one,two,three
  double precision zero,load_pl,load_ql,v_mag,v_ang,base_kv
  double precision r,x,b
  character(len=32) bstring
!
  load_pl = 0.0d00
  load_ql = 0.0d00
  v_mag = 0.0d00
  v_ang = 0.0d00
  base_kv = 1.0d00
!
  r = 0.05d00
  x = 0.0d00
  b = 0.0d00
!
  zero = 0.0d00
  one = 1
  two = 2
  three = 3
!
!  Figure out number of digits in bus name string
!
   nbus = maxdim*maxdim
   icnt = 0
   do while (nbus.gt.0)
     icnt = icnt + 1
     nbus = nbus/10
   end do
   if (icnt.lt.10) then
     write(bstring,'(a,i1,a,i1,a)') '(i',icnt,'.',icnt,',a)'
   else
     write(bstring,'(a,i2,a,i2,a)') '(i',icnt,'.',icnt,',a)'
   endif
!
!  Print out case data
!
  write(6,'(a)') '0  100.000'
  write(6,*)
  write(6,*)
!
!  print out bus block
!
  nbus = maxdim**2
  do j = 1, maxdim
    do i = 1, maxdim
      n = (j-1)*maxdim + i
      write(6,'(i9)',advance='no') n
      if ((i.eq.1.and.j.eq.1).or.(i.eq.maxdim.and.j.eq.maxdim)) then
        write(6,'(a,i4)',advance='no') ',',two
        write(6,'(a,f8.3)',advance='no') ',',zero
        write(6,'(a,f8.3)',advance='no') ',',zero
      else
        write(6,'(a,i4)',advance='no') ',',one
        write(6,'(a,f8.3)',advance='no') ',',load_pl
        write(6,'(a,f8.3)',advance='no') ',',load_ql
      endif
      write(6,'(a,f8.3)',advance='no') ',',zero
      write(6,'(a,f8.3)',advance='no') ',',zero
      write(6,'(a,i4)',advance='no') ',',one
      write(6,'(a,f8.3)',advance='no') ',',v_mag
      write(6,'(a,f8.3)',advance='no') ',',zero
      write(6,'(a)',advance='no') ',''BUS-'  ! Double apostrophe to print "'"
      write(6,bstring,advance='no') n,'   '''
      if (i.eq.1.and.j.eq.1) then
        write(6,'(a,f10.4)',advance='no') ',',base_kv
      else if (i.eq.maxdim.and.j.eq.maxdim) then
        write(6,'(a,f10.4)',advance='no') ',',-base_kv
      else
        write(6,'(a,f10.4)',advance='no') ',',zero
      end if
      write(6,'(a,i4)') ',',two
    end do
  end do
  write(6,'(a)') '0 / END OF BUS DATA, BEGIN GENERATOR DATA'
  write(6,'(a)') '0 / END OF GENERATOR DATA, BEGIN BRANCH DATA'
!
!  print out branch block. Start with branches in i-direction
!
  do j = 1, maxdim
    do i = 2, maxdim
      write(6,'(i9)',advance='no') (j-1)*maxdim + i-1
      write(6,'(a,i9)',advance='no') ',',(j-1)*maxdim + i
      write(6,'(a)',advance='no') ',''BL'''
      write(6,'(a,f10.5)',advance='no') ',',r
      write(6,'(a,f10.5)',advance='no') ',',x
      write(6,'(a,f10.5)',advance='no') ',',b
      write(6,'(a,f6.2)',advance='no') ',',zero
      write(6,'(a,f6.2)',advance='no') ',',zero
      write(6,'(a,f6.2)',advance='no') ',',zero
      write(6,'(a,f8.5)',advance='no') ',',zero
      write(6,'(a,f8.3)',advance='no') ',',zero
      write(6,'(a,f8.5)',advance='no') ',',zero
      write(6,'(a,f8.5)',advance='no') ',',zero
      write(6,'(a,f8.5)',advance='no') ',',zero
      write(6,'(a,f8.5)',advance='no') ',',zero
      write(6,'(a,i4)') ',',one
    end do
  end do
!
!  Branches in j-direction
!
  do j = 2, maxdim
    do i = 1, maxdim
      write(6,'(i9)',advance='no') (j-2)*maxdim + i
      write(6,'(a,i9)',advance='no') ',',(j-1)*maxdim + i
      write(6,'(a)',advance='no') ',''BL'''
      write(6,'(a,f10.5)',advance='no') ',',r
      write(6,'(a,f10.5)',advance='no') ',',x
      write(6,'(a,f10.5)',advance='no') ',',b
      write(6,'(a,f6.2)',advance='no') ',',zero
      write(6,'(a,f6.2)',advance='no') ',',zero
      write(6,'(a,f6.2)',advance='no') ',',zero
      write(6,'(a,f8.5)',advance='no') ',',zero
      write(6,'(a,f8.3)',advance='no') ',',zero
      write(6,'(a,f8.5)',advance='no') ',',zero
      write(6,'(a,f8.5)',advance='no') ',',zero
      write(6,'(a,f8.5)',advance='no') ',',zero
      write(6,'(a,f8.5)',advance='no') ',',zero
      write(6,'(a,i4)') ',',one
    end do
  end do
  write(6,'(a)') '0 / END OF BRANCH DATA, BEGIN TRANSFORMER ADJUSTMENT DATA'
  write(6,'(a)') '0 / END OF TRANSFORMER ADJUSTMENT DATA, BEGIN AREA DATA'
  write(6,'(a)') '0 / END OF AREA DATA, BEGIN TWO-TERMINAL DC DATA'
  write(6,'(a)') '0 / END OF TWO-TERMINAL DATA, BEGIN SWITCHED SHUNT DATA'
  write(6,'(a)') '0 / END OF SWITCHED SHUNT DATA'
  stop
end program
