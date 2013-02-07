program convert
use dict_mod
implicit none

integer i,j,k,i1,i2,ii,jj,kk,iii,kkk,kkkk,PROJECTING,ibad,changes
real*8 a1,a2, ratio, dum,dum2,dum3,aa,some
character*1 :: char1,char1b,char1c
character*3 :: char3
character*4 :: char4
character*6 :: char6
character*7 :: char7
character*9 :: char9,char9b
character*12:: char12,char12b
! Input reading
integer blanks
character(len=80)                :: string80
character(len=80), dimension(40) :: tokens
integer, dimension(40)           :: sublen
type(DICT_STRUCT), pointer :: dict_x
type(DICT_STRUCT), pointer :: dict_y
character(len=80) :: char_last
character(len=10) :: obj_name

logical foundone
logical file_exists

integer*8 count

type growable
   integer length
   integer used
   real*8 b ! The right hand side B value
   real*8, dimension(:), allocatable :: vec
   integer equals ! Is is equals or less than or greater than
   integer,dimension(:), allocatable :: ind
end type growable

type bounds
   real*8 upper
   real*8 lower
   logical bound ! Does it have bounds at all
end type bounds

type(growable), dimension(:), allocatable :: array ! A matrix, element zero will be objective function
type(bounds),  dimension(:), allocatable :: limits ! The limits on X

character(1024) :: filename
integer ncolA ! Number of X
integer nrowA ! Number of B
integer ranges
integer :: nrowA_real,slack,ncolA_real

real*8, allocatable,  dimension(:) :: slack_work
real*8, allocatable, dimension(:) :: temp_wrt
integer, allocatable, dimension(:) :: adjuster
!$OMP THREADPRIVATE (/saver/) ! This just evil code to make sure each thread gets its own copy
integer, parameter :: ncola_max=1280000
real*8,  dimension(ncola_max) :: work_vec
integer, dimension(ncola_max) :: work_ind
common /saver/ work_vec, work_ind
integer, allocatable, dimension(:) :: signs

!
!
!
! Allocate memory & setup
!
!
!

if (iargc() .ne. 1) then
   PRINT *,'Usage: convert filename'
   stop
end if
CALL getarg(1,filename)

! Open files
INQUIRE(FILE=TRIM(filename)//'.out', EXIST=file_exists)
if (file_exists) call abrt(__LINE__, TRIM(filename)//'.out'//' exists') ! Avoid segfault if file exists
open(unit=66,file=TRIM(filename)//'.out',status='new',form='unformatted',access='stream') ! open early, so we bomb early if it exists
open(unit=77, file=TRIM(filename),STATUS='OLD',ACCESS='sequential',form='formatted')

print*,' Determining size of input file ',TRIM(filename)
call flush(6)

read(77,'(a4)') char4
if (char4 .ne. 'NAME') call abrt(__LINE__)
read(77,'(a4)') char4
if (char4 .ne. 'ROWS') call abrt(__LINE__)

print*,' Counting ROWS'

nrowA = 0
do
  read(77,'(a80)') string80
  call get_strings(string80, tokens, sublen)
  if (tokens(1) .eq. 'COLUMNS') exit
  nrowA = nrowA + 1
enddo

print*,' Counting COLUMNS'
ncola = 0
char_last=''
do
  read(77,'(a80)') string80
  call get_strings(string80, tokens, sublen)
  if (sublen(1) .eq. 0) call abrt(__LINE__)
  if (tokens(1) .eq. 'RHS') exit
  if (tokens(1) .eq. char_last) then
    ! Still on same one
  else
    char_last = tokens(1) ! New one
    ncola = ncola+1
  endif
enddo
ncola = ncola+1 ! Extra one for slack

ranges = 0
do
  read(77,'(a80)') string80
  call get_strings(string80, tokens, sublen)
  if (tokens(1) .eq. 'ENDATA') goto 972
  if (tokens(1) .eq. 'RANGES') then
    char_last = ''
    do ! No way to know how long
      read(77,'(a80)') string80
      call get_strings(string80, tokens, sublen)
      if (tokens(1)(1:4) .eq. 'BOUN') goto 972
      if (tokens(1) .eq. 'ENDATA') goto 972
      if (sublen(1) .eq. 0) call abrt(__LINE__)
      if (sublen(3) .eq. 0 .or. (sublen(3) .ne. 0 .and. sublen(4) .ne. 0 .and. sublen(5) .eq. 0)) then ! Unnamed section
        if (char_last .eq. '') char_last = 'NOTHING HERE AT ALL'
        if (char_last .ne. 'NOTHING HERE AT ALL') goto 972 ! We only do first section
        blanks = 1
      else
        if (char_last .eq. '') char_last = trim(tokens(1))
        if (trim(tokens(1)) .ne. char_last) goto 972 ! We only do first section
        blanks = 0
      endif
      ranges = ranges+1
      if (sublen(4-blanks) .ne. 0 ) ranges = ranges+1
    enddo
  endif
enddo
972 continue
nrowa = nrowa + ranges

print*,ncola-1,nrowA-1 ! ignore slack and objective function

if(ncola .gt. ncola_max) then
  print*, ' ncola_max is too small. recompile code after increasing to',ncola
  stop
endif

rewind(77)

allocate(slack_work(0:nrowA),temp_wrt(ncolA),signs(ncolA))

signs(:) =1

allocate(array(0:nrowA))
do i = 0, nrowA - ranges
  allocate(array(i)%vec(ncolA/10)) ! Low memory start
  allocate(array(i)%ind(ncolA/10))
  array(i)%length = ncolA/10
  array(i)%used   = 0
  array(i)%b      = 0.0d0 ! Default
  array(i)%equals = -2    ! Invalid choice
enddo
allocate(limits(ncolA))
do i = 1, ncolA
  limits(i)%upper = 1.0D200
  limits(i)%lower = 0.0d0
  limits(i)%bound = .true.
enddo
limits(ncolA)%upper = 0.0D0
limits(ncolA)%lower = 0.0d0
limits(ncolA)%bound = .true. ! Use as a slack variable


!
!
!
! Read in the file
!
!
!

print*,' Reading in input file'
call flush(6)

! Read in NAME
read(77,'(a4)') char4
if (char4 .ne. 'NAME') call abrt(__LINE__)

! Read in ROWS
read(77,'(a4)') char4
if (char4 .ne. 'ROWS') call abrt(__LINE__)

print*,' Reading ROWS'
call flush(6)

call dict_create(dict_x)

obj_name = ''
do i = 1, nrowa-ranges
  read(77,'(a80)') string80
  call get_strings(string80, tokens, sublen)
  if(sublen(2) .eq. 0) call abrt(__LINE__)
  if(sublen(3) .ne. 0) call abrt(__LINE__)
  if(sublen(1) .ne. 1) then
    call abrt(__LINE__)
  else if (tokens(1) .eq. 'E') then
    array(i)%equals = 0 
    call dict_add_key(dict_x, trim(tokens(2)) , i) 
  else if (tokens(1) .eq. 'L') then
    array(i)%equals = -1 
    call dict_add_key(dict_x, trim(tokens(2)) , i)
  else if (tokens(1) .eq. 'G') then
    array(i)%equals = 1 
    call dict_add_key(dict_x, trim(tokens(2)) , i)
  else if (tokens(1) .eq. 'N') then
    if (obj_name .eq. '') then
      obj_name = trim(tokens(2)) ! First one
      array(0)%equals =  1 ! THIS FLIPS SIGN, BECAUSE CPLEX MINIMIZES
    else
      call dict_add_key(dict_x, trim(tokens(2)) , i)
    endif
    array(i)%equals = -4 ! Variable to ignore.  All "N" variables are always in zero position
  else 
    call abrt(__LINE__)
  endif
enddo
if (obj_name .eq. '') call abrt(__LINE__)

print*,' Reading COLUMNS'
call flush(6)

! Read in COLUMNS
call dict_create(dict_y)
read(77,'(a7)') char7
if (char7 .ne. 'COLUMNS') call abrt(__LINE__)
j = 0
do ! No way to know how long
  read(77,'(a80)') string80
  call get_strings(string80, tokens, sublen)
  if(tokens(1) .eq. 'RHS') goto 123
  if(sublen(2) .eq. 8) then
    if (tokens(2)(2:7) .eq. 'MARKER') then
      call abrt(__LINE__,'MARKER detected for mixed integer input') ! Integer not supported
    endif
  endif
  if(j .eq. 0) then
     j = 1
     call dict_add_key(dict_y, trim(tokens(1)), j)
     char_last = tokens(1)
  else if(char_last .eq. tokens(1) ) then
     ! Do nothing still on same line
  else
     j = j + 1
     call dict_add_key(dict_y, trim(tokens(1)), j)
     char_last = tokens(1)
  endif
  if(tokens(2) .eq. obj_name ) then
    i1 = 0
  else
    i1 = dict_get_key( dict_x, trim(tokens(2)) )
    if(i1 .eq. 0) call abrt(__LINE__)
  endif
  read(tokens(3),*) a1
  call add_element(j,i1,a1)

  if(sublen(4) .ne. 0) then
    if(tokens(4) .eq. obj_name ) then
      i1 = 0
    else
      i1 = dict_get_key( dict_x, trim(tokens(4)) )
      if(i1 .eq. 0) call abrt(__LINE__)
    endif
    read(tokens(5),*) a1
    call add_element(j,i1,a1)
  endif
enddo
123 continue

! Read in RHS
print*,' Reading RHS'
call flush(6)

char_last = ''
do ! No way to know how long
  read(77,'(a80)') string80
  call get_strings(string80, tokens, sublen)
  if (tokens(1)(1:4) .eq. 'BOUN') goto 456
  if (tokens(1) .eq. 'ENDATA') goto 789
  if (tokens(1) .eq. 'RANGES') goto 3456
  if (sublen(1) .eq. 0) call abrt(__LINE__)
  if (sublen(3) .eq. 0 .or. (sublen(4) .ne. 0 .and. sublen(5) .eq. 0)) then ! Unnamed section
    if (char_last .eq. '') char_last = 'NOTHING HERE AT ALL'
    if (char_last .ne. 'NOTHING HERE AT ALL') goto 4567 ! We only do first section
    blanks = 1
  else
    if (char_last .eq. '') char_last = trim(tokens(1))
    if (trim(tokens(1)) .ne. char_last) goto 4567 ! We only do first section
    blanks = 0
  endif
  i1 = dict_get_key( dict_x, trim(tokens(2-blanks)) )
  if (i1 .eq. 0) call abrt(__LINE__)
  read(tokens(3-blanks),*) a1
  array(i1)%b = a1
  if (sublen(4-blanks) .ne. 0 ) then
    i2 = dict_get_key( dict_x, trim(tokens(4-blanks)) )
    if (i2 .eq. 0) call abrt(__LINE__)
    read(tokens(5-blanks),*) a2
    array(i2)%b = a2
  endif
  4567 continue
enddo

3456 continue
print*,' Reading RANGES'
call flush(6)

! Read in RANGES
char_last = ''
i = nrowa - ranges ! put them after
do ! No way to know how long.  We know how many, not how many lines
  read(77,'(a80)') string80
  call get_strings(string80, tokens, sublen)
  if (tokens(1)(1:4) .eq. 'BOUN') goto 456
  if (tokens(1) .eq. 'ENDATA') goto 789
  if (sublen(1) .eq. 0) call abrt(__LINE__)
  if (sublen(3) .eq. 0 .or. (sublen(3) .ne. 0 .and. sublen(4) .ne. 0 .and. sublen(5) .eq. 0)) then ! Unnamed section
    if (char_last .eq. '') char_last = 'NOTHING HERE AT ALL'
    if (char_last .ne. 'NOTHING HERE AT ALL') goto 5678 ! We only do first section
    blanks = 1
  else
    if (char_last .eq. '') char_last = trim(tokens(1))
    if (trim(tokens(1)) .ne. char_last) goto 5678 ! We only do first section
    blanks = 0
  endif
  i1 = dict_get_key( dict_x, trim(tokens(2-blanks)) )
  if (i1 .eq. 0) call abrt(__LINE__)
  read(tokens(3-blanks),*) a1
  i = i + 1
  if      (array(i1)%equals .eq. -4) then 
    call abrt(__LINE__)
  else if (array(i1)%equals .eq. -2) then
    call abrt(__LINE__)
  else if (array(i1)%equals .eq. -1) then
     array(i)%equals = 1
     array(i)%b = array(i1)%b-abs(a1)
  else if (array(i1)%equals .eq.  1) then
     array(i)%equals = -1
     array(i)%b = array(i1)%b+abs(a1)
  else if (array(i1)%equals .eq.  0 .and. a1 .ge. 0) then
     array(i)%equals = -1
     array(i1)%equals = 1
     array(i)%b = array(i1)%b+a1
  else if (array(i1)%equals .eq.  0 .and. a1 .le. 0) then
     array(i)%equals = 1
     array(i1)%equals = -1
     array(i)%b = array(i1)%b+a1
  else
    call abrt(__LINE__)
  endif
  array(i)%used = array(i1)%used
  array(i)%length = array(i1)%length
  allocate(array(i)%vec(array(i)%length))
  allocate(array(i)%ind(array(i)%length))
  array(i)%vec(1:array(i1)%used) = array(i1)%vec(1:array(i1)%used)
  array(i)%ind(1:array(i1)%used) = array(i1)%ind(1:array(i1)%used)
  if (sublen(4-blanks) .ne. 0 ) then
    i1 = dict_get_key( dict_x, trim(tokens(4-blanks)) )
    if (i1 .eq. 0) call abrt(__LINE__)
    read(tokens(5-blanks),*) a1
    i = i + 1
    if      (array(i1)%equals .eq. -4) then
      call abrt(__LINE__)
    else if (array(i1)%equals .eq. -2) then
      call abrt(__LINE__)
    else if (array(i1)%equals .eq. -1) then
       array(i)%equals = 1
       array(i)%b = array(i1)%b-abs(a1)
    else if (array(i1)%equals .eq.  1) then
       array(i)%equals = -1
       array(i)%b = array(i1)%b+abs(a1)
    else if (array(i1)%equals .eq.  0 .and. a1 .ge. 0) then
       array(i)%equals = -1
       array(i1)%equals = 1
       array(i)%b = array(i1)%b+a1
    else if (array(i1)%equals .eq.  0 .and. a1 .le. 0) then
       array(i)%equals = 1
       array(i1)%equals = -1
       array(i)%b = array(i1)%b+a1
    else
      call abrt(__LINE__)
    endif
    array(i)%used = array(i1)%used
    array(i)%length = array(i1)%length
    allocate(array(i)%vec(array(i)%length))
    allocate(array(i)%ind(array(i)%length))
    array(i)%vec(1:array(i1)%used) = array(i1)%vec(1:array(i1)%used)
    array(i)%ind(1:array(i1)%used) = array(i1)%ind(1:array(i1)%used)
  endif
  5678 continue
enddo
456 continue

print*,' Reading BOUNDS'
call flush(6)

! Read in BOUNDS

char_last = ''
do ! No way to know how long
  read(77,'(a80)') string80
  call get_strings(string80, tokens, sublen)
  if (sublen(1) .eq. 0) call abrt(__LINE__)
  if (tokens(1) .eq. 'SOS') call abrt(__LINE__) ! We do not so specialy ordered sets
  if (tokens(1) .eq. 'ENDATA') goto 789
  if ((sublen(4) .eq. 0 .and. tokens(1) .ne. 'FR') .or. sublen(3) .eq. 0 ) then ! No name at all
    if (char_last .eq. '') char_last = 'NOTHING HERE AT ALL'
    if (char_last .ne. 'NOTHING HERE AT ALL') goto 7890 ! We only do first section
    blanks = 1
  else
    if (char_last .eq. '') char_last = trim(tokens(2))
    if (trim(tokens(2)) .ne. char_last) goto 7890 ! We only do first section
    blanks = 0
  endif
  i1 = dict_get_key( dict_y, trim(tokens(3-blanks)))
  if(i1 .eq. 0) call abrt(__LINE__)
  if (tokens(1) .eq. 'LO') then
    read(tokens(4-blanks),*) limits(i1)%lower
  else if (tokens(1) .eq. 'UP') then
    read(tokens(4-blanks),*) limits(i1)%upper
  else if (tokens(1) .eq. 'FX') then
    read(tokens(4-blanks),*) limits(i1)%upper
    limits(i1)%lower = limits(i1)%upper
  else if (tokens(1) .eq. 'FR') then
    limits(i1)%upper = 1.0D200
    limits(i1)%lower =-1.0D200
    limits(i1)%bound = .false.
  else
    call abrt(__LINE__)
  endif
  7890 continue
enddo

! Read in ENDDATA
789 continue ! Already done

! Close file
close(unit=77) !TODO add support for various case keywords (not tokens, which are case sensitve!)

!
!
!
!
! How sparse is this thing?
!
!
!
!

count = 0
do i = 1, nrowA
  count = count + array(i)%used
enddo
print*
print*, 'nrowA', nrowA-1 ! Ignore objective function
print*, 'ncolA', ncolA-1 ! Ignore slack
print*, 'non-zero elements', count
print*

!
!
!
! Flip signs of greater thans
!
!
!

print*,' Turning greater than into less than'
call flush(6)

do i = 0, nrowA
  if (array(i)%equals .eq. 1) then
    do j = 1, array(i)%used
      array(i)%vec(j) = -array(i)%vec(j)
    enddo
    array(i)%b = -array(i)%b
    array(i)%equals = -1
  else if (array(i)%equals .eq. -4) then
    if (allocated(array(i)%vec)) deallocate(array(i)%vec,array(i)%ind)
    array(i)%length = 0
    array(i)%used   = 0
    array(i)%b      = 0.0d0
    array(i)%equals = -3
  else if (array(i)%equals .eq. -2) then ! Unset
    call abrt(__LINE__)
  endif
enddo 

!
!
!
!
! Project out fixed variables
!
!
!
!

print*, ' Projecting out fixed variables'
call flush(6)

do j = 1, ncolA
  if(limits(j)%lower .eq. limits(j)%upper) then
    do i = nrowA,1,-1
      do jj = 1, array(i)%used
        if (array(i)%ind(jj) .eq. j)then
          array(i)%b = array(i)%b - array(i)%vec(jj)*limits(j)%lower ! subtract it off
          array(i)%ind(jj:array(i)%used-1) = array(i)%ind(jj+1:array(i)%used)
          array(i)%vec(jj:array(i)%used-1) = array(i)%vec(jj+1:array(i)%used)
          array(i)%used = array(i)%used - 1
          goto 9378 ! Done, no duplicates
        else if (array(i)%ind(jj) .gt. j) then
          goto 9378 ! Done, they are sorted
        endif
      enddo
 9378 continue
    enddo
  endif
enddo

!
!
!
! Project out equalities
!
!
!

print*, ' Projecting out equalities'
call flush(6)

nrowA_real = nrowA
do i = nrowA,1,-1
  if(mod(i,2000) .eq. 0) then
    print*,'start',i
    call flush(6)
  endif
  projecting = -1
  if (array(i)%equals .eq. 0) then ! Found an equality
    do j = 1, array(i)%used
      if(.not. limits(array(i)%ind(j))%bound) then ! Found a free variable that can be projected out 
        if (projecting .eq. -1) then
          projecting = j
        else if (abs(array(i)%vec(j)) .gt. abs(array(i)%vec(projecting)))then
          projecting = j ! project biggest one
        endif
      endif
    enddo
    if(projecting .ne. -1) then
        j = projecting
        projecting = array(i)%ind(projecting)
!$OMP PARALLEL DO PRIVATE(k,foundone,kkk,kkkk,kk,ratio,ii,jj,iii) DEFAULT (NONE) SHARED(i,j,projecting,array,nrowA_real,ncolA) SCHEDULE(DYNAMIC,8)
        do k = 0, nrowA_real ! Loop over all others and project
          if(k .ne. i .and. array(k)%equals .ne. -3 .and. array(k)%equals .ne. -5) then ! Found one that might be projectable
            foundone = .false.
            kkk  = 1
            kkkk = array(k)%used
            do ! Binary search, because we are sorted
             kk = (kkkk + kkk)/2
             if (array(k)%ind(kk) .eq. projecting)  then
                foundone = .true.
                goto 414
             else if (array(k)%ind(kk) .lt. projecting)  then
                kkk = kk + 1
             else
                kkkk = kk - 1
             endif
             if (kkkk .lt. kkk) goto 414
            enddo
 414        continue
            if (foundone) then
               ratio = -array(k)%vec(kk)/array(i)%vec(j)
               ii = 1
               kk = 1
               jj = 1
               iii = -1
               kkk = -1
               do while (iii .lt. 100000000 .or. kkk .lt. 100000000)
                if (ii .le. array(i)%used) then
                  iii = array(i)%ind(ii)
                else
                  iii = 110000000 ! We should be so lucky to handle a system this big  :-)
                endif
                if (kk .le. array(k)%used) then
                  kkk = array(k)%ind(kk)
                else
                  kkk = 110000000
                endif
! Add both k&i to temporary work array (keeps elements in order) ie. sorted
                if (iii .eq. kkk) then
                  if(iii .ge. 100000000) then
                    ! Done
                  else if (iii .eq. projecting) then
                    ! We know that it is zero, don't let noise say otherwise
                    kk = kk + 1
                    ii = ii + 1
                  else
                    work_ind(jj) = iii ! Same as kkk
                    work_vec(jj) = array(i)%vec(ii)*ratio + array(k)%vec(kk)
                    kk = kk + 1
                    ii = ii + 1
                    if(abs(work_vec(jj)) .gt. 1.0D-20) jj = jj + 1 ! Only counts if not zero
                  endif
                else if (iii .gt. kkk) then
                  work_ind(jj) = kkk
                  work_vec(jj) = array(k)%vec(kk)
                  kk = kk + 1
                  if(abs(work_vec(jj)) .gt. 1.0D-20) jj = jj + 1 ! Only counts if not zero
                else if (iii .lt. kkk) then
                  work_ind(jj) = iii
                  work_vec(jj) = array(i)%vec(ii)*ratio
                  ii = ii + 1
                  if(abs(work_vec(jj)) .gt. 1.0D-20) jj = jj + 1 ! Only counts if not zero
                endif
              enddo
! Set k to this work array
              jj = jj - 1
              if (jj .gt. ncolA) call abrt(__LINE__)
              if (jj .gt. array(k)%length ) then !Growing
                deallocate(array(k)%vec,array(k)%ind)
                array(k)%length=max(2*array(k)%length,jj+1000)
                array(k)%length=min(array(k)%length,ncolA)
                allocate(array(k)%vec(array(k)%length),array(k)%ind(array(k)%length))
              endif
              array(k)%vec(1:jj) = work_vec(1:jj)
              array(k)%ind(1:jj) = work_ind(1:jj)
              array(k)%used = jj
! Do the b variables too
              array(k)%b = array(k)%b + ratio*array(i)%b
            endif
          endif
        enddo
        deallocate(array(i)%vec,array(i)%ind)
        array(i)%length = 0
        array(i)%used   = 0
        array(i)%b      = 0.0d0
        array(i)%equals = -3
        limits(projecting)%bound = .true.
        limits(projecting)%upper = 0.0d0
        limits(projecting)%lower = 0.0d0
    else
      if (array(i)%used .ne. 0) then
        print*,'Found equality that cannot be removed. #',i ! Only can remove variables that are free
        call flush(6)
        array(i)%equals = -5
      else if (array(i)%used .ne. 0 .or. array(i)%b .ne. 0.0d0)  then
        print*,'Partially destoyed equality found.  Strange. #',i,array(i)%used,array(i)%b
        call flush(6)
        call abrt(__LINE__) ! Really should never happen
      else
        array(i)%length = 0
        array(i)%used   = 0
        array(i)%b      = 0.0d0
        array(i)%equals = -3
        if(allocated(array(i)%vec)) deallocate(array(i)%vec,array(i)%ind)
      endif
    endif
    if(nrowA_real .eq. i) nrowA_real=i-1 ! if all we have seen is equalities, then we can shorten the searching
  endif
enddo

!
!
!
! Remove negatives from objective function to be maximized (we flipped sign earlier, since input is minimization problem)
!
!
!

print*,' Removing negatives from objective function'
call flush(6)

temp_wrt(:) = 0.0d0 ! Will hold full copy of original objective function

do j=1,array(0)%used ! Objective vector (known as f and c)
  k = array(0)%ind(j)
  dum = array(0)%vec(j)
  temp_wrt(k) = dum
  array(0)%vec(j) = abs(dum)
enddo

do i=1,nrowA ! A Matrix
  if(allocated(array(i)%vec)) then
    do j=1,array(i)%used
      k = array(i)%ind(j)
      if(temp_wrt(k) .lt. 0) then
        array(i)%vec(j) = -array(i)%vec(j)
      endif
    enddo
  endif
enddo

do j=1,ncolA ! Explicit constraints
  if (limits(j)%bound) then
    if (temp_wrt(j) .lt. 0) then
      dum = limits(j)%lower
      limits(j)%lower = -limits(j)%upper
      limits(j)%upper = -dum
    endif
  endif
enddo

do j=1,ncolA
  if (temp_wrt(j) .lt. 0) then
     signs(j) = -1
  else
     signs(j) = 1
  endif
enddo

!
!
!
! Make sure that all 'b' variables are positive definite
!
!
!

print*,' Making b positive definite by adding a fixed slack varible'
call flush(6)
! Find a place to put the slack variable
do j = ncolA,1,-1 ! Use one nearest end, to make insertion easier
  if (limits(j)%bound .and. limits(j)%upper .eq. 0.0d0 .and. limits(j)%lower .eq. 0.0d0)  then
    goto 939
  endif
enddo
call abrt(__LINE__) ! Did not find a variable to use as slack
939 continue
slack = j
limits(slack)%upper = 1.0d0
limits(slack)%lower = 1.0d0
slack_work(0:nrowA) = 0.0d0
! Use it
do i = 1, nrowA
  if (array(i)%equals .ne. -3 .and. array(i)%b .le. 0.0d0) then
     dum = 0.1d0 - array(i)%b
     array(i)%b = array(i)%b + dum
     slack_work(i) = slack_work(i) + dum
  endif
enddo
! Add it
do i = 0, nrowA
  if (slack_work(i) .ne. 0) then
    call insert_element(slack,i,slack_work(i))
  endif
enddo

!
!
!
! Pre-process inequalities
!
!
!


! For example, if x+y+z<10 and x,y,z<0, then we know that x,y,z<10.
! The turns limits that are implicit into ones that are explicit

print*,'Pre-processing inequalities'
call flush(6)

933 continue

some=1.0d10
do while(some .gt. 0.0001d0)
   print*,'CYCLE'
   call flush(6)
   some=0.0d0
   do i=1,nrowA
     dum = 0
     do j=1,array(i)%used
       jj = array(i)%ind(j)
       aa = array(i)%vec(j)
       if (aa .eq. 0) then
         ! Nothing - strange
       else if (aa .gt. 0) then
         if (limits(jj)%lower .lt. -1.0d60) goto 911
         dum = dum + limits(jj)%lower*aa
       else ! aa .lt. 0
         if ( limits(jj)%upper .gt. 1.0d60) goto 911
         dum = dum + limits(jj)%upper*aa
       endif
     enddo
     dum = dum - array(i)%b
     do j=1,array(i)%used
       jj = array(i)%ind(j)
       aa = array(i)%vec(j)
       if (aa .eq. 0) then
         ! Nothing
       else if (aa .gt. 0) then
         dum3 = dum - limits(jj)%lower*aa! dum2 without j part
         dum2 = -dum3/aa
         if (dum2 .lt. limits(jj)%upper) then
           some=some+limits(jj)%upper-dum2
           limits(jj)%upper = dum2
           if(limits(jj)%upper .lt. limits(jj)%lower) then
             print*,'freezing out',jj
             limits(jj)%upper = limits(jj)%lower
           endif
         endif
       else
         dum3 = dum - limits(jj)%upper*aa! dum2 without j part
         dum2 = -dum3/aa
         if (dum2 .gt. limits(jj)%lower) then
           some=some+dum2-limits(jj)%lower
           limits(jj)%lower = dum2
           if(limits(jj)%upper .lt. limits(jj)%lower) then
             print*,'freezing out',jj
             limits(jj)%lower = limits(jj)%upper
           endif
         endif
       endif
     enddo
 911 continue ! Found an infinite one
   enddo
enddo


some=1.0d10
changes=0
do while(some .gt. 0.001d0)
   print*,'CYCLING'
   call flush(6)
   some=0.0d0
   do i=1,nrowA
     dum = 0
     ibad = 0
     do j=1,array(i)%used
       jj = array(i)%ind(j)
       aa = array(i)%vec(j)
       if (aa .eq. 0) then
         ! Nothing - strange
       else if (aa .gt. 0) then
         if (limits(jj)%lower .lt. -1.0d60) then
            if(ibad .ne.0) goto 922 ! More than one
            ibad = j
         else
            dum = dum + limits(jj)%lower*aa
         endif
       else ! aa .lt. 0
         if ( limits(jj)%upper .gt. 1.0d60) then
            if(ibad .ne.0) goto 922 ! More than one
            ibad = j
         else
            dum = dum + limits(jj)%upper*aa
         endif
       endif
     enddo
     if (ibad .ne.0) then
       dum = dum - array(i)%b
       j = ibad
       jj = array(i)%ind(j)
       aa = array(i)%vec(j)
       dum2 = -dum/aa
       if (aa .gt. 0) then
           if (dum2 .lt. limits(jj)%upper) then
           some=some+limits(jj)%upper-dum2
           limits(jj)%upper = dum2
           changes = changes+1
         endif
       else
         if (dum2 .gt. limits(jj)%lower) then
           some=some+dum2-limits(jj)%lower
           limits(jj)%lower = dum2
           changes=changes+1
         endif
       endif
     endif
 922 continue ! Found an infinite one
   enddo
enddo
if(changes .gt. 0) goto 933 ! Time to start all over


!
!
!
! Save out results of Objective function and A Matrix, etc.
!
!
!

print*,' Saving results to output file'
call flush(6)

nrowA_real = 0
do i=1,nrowA
  if(allocated(array(i)%vec)) then
     nrowA_real = nrowA_real + 1
  endif
enddo

ncolA_real = 0
do i = 1,ncolA
  if(limits(i)%upper .ne. 0 .or. limits(i)%lower .ne. 0) then
    ncolA_real = ncolA_real + 1
  endif
enddo

write(66) nrowA_real, ncolA_real

allocate(adjuster(ncola))
adjuster(:) = -999999
k=0
do j=1,ncolA
  if(limits(j)%upper .ne. 0 .or. limits(j)%lower .ne. 0 ) then
    k=k+1
    adjuster(j) = k ! Map ncolA to ncola_real
  endif
enddo  
if (ncola_real .ne. k) call abrt(__LINE__)


! Remove variables that are fixed to zero by limits
do i=1,nrowA
  if(allocated(array(i)%vec)) then
    do j = array(i)%used,1,-1
      if (adjuster(array(i)%ind(j)) .lt. 0) Then
        do jj= j+1, array(i)%used
           array(i)%ind(jj-1) = array(i)%ind(jj)
           array(i)%vec(jj-1) = array(i)%vec(jj)
        enddo
        array(i)%used = array(i)%used - 1
      endif
    enddo
  endif
enddo

! Save out B vector
do i=1,nrowA
  if(allocated(array(i)%vec)) then
    write(66)  array(i)%b
  endif
enddo

! Save out limits
do i = 1, ncolA
  if(limits(i)%upper .ne. 0 .or. limits(i)%lower .ne. 0 ) then
    if (limits(i)%bound) then
      write(66) limits(i)%upper
    else
      write(66) 1.0d120
    endif
  endif
enddo

do i = 1, ncolA
  if(limits(i)%upper .ne. 0 .or. limits(i)%lower .ne. 0 ) then
    if (limits(i)%bound) then
      write(66) limits(i)%lower
    else
      write(66) -1.0d120
    endif
  endif
enddo

! Save out internal sign flips
do i = 1, ncolA
  if(limits(i)%upper .ne. 0 .or. limits(i)%lower .ne. 0 ) then
    write(66) signs(i)
  endif
enddo

! Save out internal variables shifts: in case we add any
do i = 1, ncolA
  if(limits(i)%upper .ne. 0 .or. limits(i)%lower .ne. 0 ) then
    write(66) 0.0d0 ! Shifter - not used, but might be in the future
  endif
enddo

! Save out holes in X
do j=1,ncolA
  if(limits(j)%upper .ne. 0 .or. limits(j)%lower .ne. 0 ) then
    write(66) j
  endif
enddo

! Save out equalities
do i = 1, nrowA
  if(allocated(array(i)%vec)) then
    foundone = array(i)%equals .eq. -5
    write(66) foundone
  endif
enddo

! Write out f
temp_wrt(:) = 0.0d0
do j=1,array(0)%used
  temp_wrt(array(0)%ind(j)) = array(0)%vec(j)
enddo
k=0
do j=1,ncolA
  if(limits(j)%upper .ne. 0 .or. limits(j)%lower .ne. 0 ) then
    k=k+1
    temp_wrt(k) = temp_wrt(j)
  endif
enddo
if(k.ne.ncolA_real) call abrt(__LINE__)
write(66)temp_wrt(1:ncolA_real)

! Save out the data proper
do i=1,nrowA
  if(mod(i,1000) .eq. 0) then
     print*,' writing',i
     call flush(6)
  endif
  if(allocated(array(i)%vec)) then
    if (array(i)%equals .eq. 0) then
      call abrt(__LINE__) ! Should not be any left
    else if (array(i)%equals .eq.-1) then
      ! OK less than 
    else if (array(i)%equals .eq. -5) then
      ! OK equality
    else if (array(i)%equals .eq. 1) then
      call abrt(__LINE__) ! Should have been flipped after reading in file
    else if (array(i)%equals .eq. -2) then
      call abrt(__LINE__) ! How did we not read anything in?
    else
      call abrt(__LINE__) ! What?
    endif
    write(66) array(i)%used
    write(66) adjuster(array(i)%ind(1:array(i)%used))
    write(66) array(i)%vec(1:array(i)%used)
  endif
enddo

! End
close(unit=66)

if (array(0)%b .ne. 0.0d0) then
  print*,'DANGER THE FVAL B is NON-ZERO ', array(0)%b
  call abrt(__LINE__)
endif
call flush(6)
stop

!
!
!
! DONE
!
!
!

contains

subroutine abrt(err,err_str)
implicit none
integer err
character*(*), optional :: err_str
write(*,'(a,i)',advance='no') 'Died on line ', err
if (present(err_str)) write(*,'(a,a)',advance='no') ' because of ', err_str
write(*,*)
call flush(6)
stop
end subroutine abrt

subroutine add_element(inp_x,inp_b,inp_a)
implicit none
integer, intent(in) :: inp_x,inp_b
real*8, intent(in) :: inp_a
integer loc_int,loc_len,loc_len_old
real*8, allocatable, dimension(:),save :: temp_vec
integer, allocatable,dimension(:),save :: temp_ind
if (.not. allocated(temp_vec)) then
  allocate(temp_ind(ncolA),temp_vec(ncolA))
endif


if (array(inp_b)%length .eq. array(inp_b)%used) then ! GROWING!!!
  loc_len = array(inp_b)%length
  temp_vec(1:loc_len) = array(inp_b)%vec(1:loc_len)
  temp_ind(1:loc_len) = array(inp_b)%ind(1:loc_len)
  deallocate(array(inp_b)%vec,array(inp_b)%ind)
  loc_len_old = loc_len
  loc_len = loc_len*2
  if (loc_len .gt. ncolA) loc_len=ncolA
  allocate(array(inp_b)%vec(loc_len),array(inp_b)%ind(loc_len))
  array(inp_b)%vec(1:loc_len_old) = temp_vec(1:loc_len_old)
  array(inp_b)%ind(1:loc_len_old) = temp_ind(1:loc_len_old)
  array(inp_b)%length = loc_len
endif

loc_int = array(inp_b)%used + 1
if(loc_int .gt. ncolA) call abrt(__LINE__)
array(inp_b)%used = loc_int
array(inp_b)%vec(loc_int) = inp_a
array(inp_b)%ind(loc_int) = inp_x

! This code assumes sorted and no duplicates
if (array(inp_b)%used .gt. 1) then
  if(array(inp_b)%ind(loc_int) .le. array(inp_b)%ind(loc_int-1)) call abrt(__LINE__)
endif

end subroutine add_element


subroutine insert_element(inp_x,inp_b,inp_a)
implicit none
integer, intent(in) :: inp_x,inp_b
real*8, intent(in) :: inp_a
real*8 loc_dum2
integer loc_int,loc_len,loc_len_old
integer loc_i,loc_j,loc_dum1
real*8, allocatable, dimension(:),save :: temp_vec
integer, allocatable,dimension(:),save :: temp_ind
if (.not. allocated(temp_vec)) then
  allocate(temp_ind(ncolA),temp_vec(ncolA))
endif

if (array(inp_b)%length .eq. array(inp_b)%used) then ! GROWING!!!
  loc_len = array(inp_b)%length
  temp_vec(1:loc_len) = array(inp_b)%vec(1:loc_len)
  temp_ind(1:loc_len) = array(inp_b)%ind(1:loc_len)
  deallocate(array(inp_b)%vec,array(inp_b)%ind)
  loc_len_old = loc_len
  loc_len = loc_len*2
  if (loc_len .gt. ncolA) loc_len=ncolA
  allocate(array(inp_b)%vec(loc_len),array(inp_b)%ind(loc_len))
  array(inp_b)%vec(1:loc_len_old) = temp_vec(1:loc_len_old)
  array(inp_b)%ind(1:loc_len_old) = temp_ind(1:loc_len_old)
  array(inp_b)%length = loc_len
endif

loc_int = array(inp_b)%used + 1
if(loc_int .gt. ncolA) call abrt(__LINE__)
array(inp_b)%used = loc_int
array(inp_b)%vec(loc_int) = inp_a
array(inp_b)%ind(loc_int) = inp_x

! This code assumes sorted already, stay that way (Code is much faster if we add them in order)
if (array(inp_b)%used .gt. 1) then
  if(array(inp_b)%ind(loc_int) .le. array(inp_b)%ind(loc_int-1)) then ! not biggest one.  crud
    loc_j = loc_int
    do loc_i = loc_int-1,1,-1
      if (array(inp_b)%ind(loc_int) .lt. array(inp_b)%ind(i)) loc_j = loc_i ! first one that is bigger
      if (array(inp_b)%ind(loc_int) .eq. array(inp_b)%ind(i)) call abrt(__LINE__) ! This is insert, not replace
    enddo
    loc_dum1 = array(inp_b)%ind(loc_int)
    loc_dum2 = array(inp_b)%vec(loc_int)
    do loc_i = loc_int-1, loc_j, -1
       array(inp_b)%ind(loc_i+1) = array(inp_b)%ind(loc_i)
       array(inp_b)%vec(loc_i+1) = array(inp_b)%vec(loc_i)
    enddo
    array(inp_b)%ind(loc_j) = loc_dum1
    array(inp_b)%vec(loc_j) = loc_dum2
  endif
endif

end subroutine insert_element

end program convert
