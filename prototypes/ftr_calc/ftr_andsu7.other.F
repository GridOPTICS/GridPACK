!------------------------------------------

!           FTR simulation nonlinear dynamical system solver
!
!------------------------------------------
!
!      This program performs the FTR optimization after it reading an FTR setup file.
!      This version is based on Matlab code written by Karanjit Kalsi,
!      based on the paper by Khanh V. Nguyen, 'A Nonlinear Neural Network
!      for Solving Linear Programming Problems,' International Symposium on
!      Mathematical Programming 2000, Atlanta, GA, August 7-11,2000.
!      Fortran implementation by Stephen Elbert.
!
!-------------------------------------------
!      Revision log:
!      Steve Elbert 2011 01 12 create separate setup program
!      Steve Elbert 2011 01 17 remove explicit transpose
!      Steve Elbert 2011 01 21 change printout
!      Steve Elbert 2011 02 07 read convergence threshold from input; change dt up/down heuristic
!      Steve Elbert 2011 02 16 introduce dtscale, sp_err_max, sp_err_it, thresh
!      Steve Elbert 2011 02 17 enforce dtmin
!      Steve Elbert 2011 02 21 put floor on thresh (.2)
!      Steve Elbert 2011 02 23 converge for 2 iterations, raise thresh floor to .5
!      Steve Elbert 2011 02 24 fmax and fval must be within err_max, thresh floor to .4
!      Steve Elbert 2011 04 04 add ACCESS='STREAM' to OPEN
!      Steve Elbert 2011 07 07 evolve code to be compatable with dense matrices
!                           12 use A(transpose) to achieve unit stride
!                           20 introduce sp_err_min
!      Steve Elbert 2011 11 11 implement adaptive nonlinear dynamical system from unds(andsu)
!                           17 implement optional bids for single interval (andsu2)
!                           22 fix initial iteration, clean up, improve output (andsu3)
!                           23 read profile to guide convergence (andsu5)
!                           23 use DDOT (Goto 2% slower @50k) and DAXPY (19% faster)(andsu6)
!                           30 use DGEMV (andsu7)


      PROGRAM FTR_andsu  ! use only unique (transposed) elements of A

      IMPLICIT NONE

      include 'mpif.h'

      LOGICAL :: more_data, read_xx

      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
      INTEGER, PARAMETER :: lu_input = 50
      INTEGER, PARAMETER :: lu_setup = 10
      INTEGER, PARAMETER :: lu_xx    = 20
      INTEGER, PARAMETER :: BUFF_ROWS = 16
      INTEGER :: i,ii,ij,j,jj,k,num_iter,itfmax,iter,iihi,max_ii,iii
      INTEGER :: sp_err_it, sp_err_it0,iter_close,nrowA,ncolA,nrow0A,namelen,nprofile
      INTEGER, ALLOCATABLE, DIMENSION(:) :: itprofile, iiprofile
      INTEGER, ALLOCATABLE, DIMENSION(:) :: signs,mapper
      INTEGER*8 :: memwrds,bigint1,bigint2,count

      type growable
         integer length
         integer used
         real(kind=dp), dimension(:), allocatable :: vec
         integer,dimension(:), allocatable :: ind
      end type growable

      real(kind=dp), dimension(:), allocatable :: vec
      integer,dimension(:), allocatable :: ind
      real(kind=dp), dimension(:), allocatable :: temp_vec
      integer,dimension(:), allocatable :: temp_ind

      integer used
      integer loc_len,loc_len_old,loc_int

      CHARACTER*120:: filename, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9
      CHARACTER*150:: filename_input,filename_setup,filename_out
      CHARACTER*10  :: hour, date

      REAL(4) t0,t1,times(2), ETIME
      REAL(KIND=DP) dxx_max,dyy_max,dt,fval,fmax,kx,ky,dtmax,dt_ky,dt_kx
      REAL(KIND=DP) k1,k2,sp_err,sp_err_prev,delerr,err_max
      REAL(KIND=DP) sp_err_max, sp_err_min, thresh, memGB, gain, gain0, fval_scl,dval_shift
      REAL(KIND=DP) delta_fval, fval0, dval, delta_dval, dval0, fval_close, dval_close
      REAL(KIND=DP), ALLOCATABLE,DIMENSION(:) :: f,bnn,max_xx,min_xx,shifter
      type(growable), ALLOCATABLE,DIMENSION(:) :: A_loci, At_locj
      INTEGER, PARAMETER :: diis_stepx = 10, diis_stepy = 10
      INTEGER, PARAMETER :: diis_lwork =  diis_stepx * diis_stepx + diis_stepy * diis_stepy
      REAL(KIND=DP), ALLOCATABLE,DIMENSION(:,:) :: diis_dxx,diis_dyy,diis_xa,diis_ya,diis_x,diis_y,diis_tmpx,diis_tmpy
      REAL(KIND=DP), ALLOCATABLE,DIMENSION(:) :: diis_yw,diis_xw
      REAL(KIND=DP) :: diis_work(diis_lwork)
      REAL(KIND=DP), ALLOCATABLE,DIMENSION(:) :: xx,yy,dxx,dyy,aty,xp,yp
      REAL(KIND=DP), ALLOCATABLE,DIMENSION(:) :: gprofile,dtprofile
      REAL(KIND=DP) :: l2norm4
      REAL(KIND=DP) :: norm2 ! This is built-into Fortran standard now
      INTEGER :: DIIS_IPIV(diis_stepx+1+diis_stepy)
      LOGICAL, ALLOCATABLE,DIMENSION(:) :: yy_can_be_neg ! For equalities, YY can be any value (think of it is positive, but with all signs flipped)

      INTEGER iargc
      INTEGER gain_cntrl
      INTEGER ierr,rank,numtasks,nsize,nextra,istart,iend,jstart,jend,style,numnodes
      INTEGER, ALLOCATABLE,DIMENSION(:) :: istarter,ilength,jstarter,jlength

!
!     -------------- End of declarations  -------------
!

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
      call print_parallel_info()


! --- read input file name, number of iterations, and time step (dt) from command argument
!

      if (iargc() .ne. 5) then
         if(rank .eq. 0) PRINT *,'Usage: FTR_nds <filename><num_iter><dt><gain><max_ii>'
         stop
      end if

      date = ''
      hour = ''
      CALL DATE_AND_TIME (date,hour)
      if(rank .eq. 0) PRINT *,'running andsu7 on ',date,' at ', hour

      err_max = 1.0d0
      CALL getarg(1,filename)
      namelen = LEN(TRIM(filename)) ! Check for existing .bin
      if (namelen .ge. 4) then
         if (filename(namelen-3:namelen) .eq. '.bin') then
            filename(1:namelen) = filename(1:namelen-4)//'    '
            filename = trim(filename)
         endif
      endif
      filename_setup = TRIM(filename)//'.bin'


      if(rank .eq. 0) PRINT *, 'Input will be read from ',filename_setup
      CALL getarg(2,arg2)
      read (arg2,*)             num_iter
      CALL getarg(3,arg3)
      read (arg3,*)             dtmax
      CALL getarg(4,arg4)
      read (arg4,*)             gain0
      CALL getarg(5,arg5)
      read (arg5,*)             max_ii

      if(rank .eq. 0) then
       PRINT *, 'Maximum number of iterations is',num_iter
       PRINT *, 'Gain is',gain0
       PRINT *, 'Maximum number of inner iterations is',max_ii
       PRINT *, 'Step size is',dtmax
       PRINT *, 'err_max is',err_max
      endif
      dt = dtmax

!     Read in optional gain changes from *.profile file
      if(rank .eq. 0)call read_profile()
      nprofile = 1

!
! ---- read setup arrays from filename_setup
!
      OPEN(UNIT=lu_setup, FILE=filename_setup, STATUS='OLD', ERR=911, ACCESS='STREAM',FORM='UNFORMATTED',readonly)

      READ(lu_setup) nrow0A, ncolA
      nrowA = nrow0A
      if(rank .eq. 0) then
       PRINT *, 'nrow0A',nrow0A,'ncolA',ncolA
       PRINT *, 'the number of constraints is', nrowA
       PRINT *, 'the number of variables is  ',ncolA
       PRINT 747, ' Input was: ',TRIM(filename),num_iter,dtmax,gain0,max_ii
      endif
 747  FORMAT(A,A,1X,I,1X,F15.9,1X,F8.2,1X,I)

! Send out how to break up loops over 1 to nrow0A
      nsize = nrow0A/(numtasks)
      nextra= nrow0A - (nsize*(numtasks))
      allocate (istarter(0:numtasks-1),ilength(0:numtasks-1))
      istarter(:) = 1
      ilength(:)  = 0
      j = 0
      if (nextra .gt. j) then
         ilength(j) = nsize + 1
      else
         ilength(j) = nsize
      endif
      do j = 1,numtasks-1
         istarter(j) = istarter(j-1) + ilength(j-1)
         if (nextra .gt. j) then
            ilength(j) = nsize + 1
         else
            ilength(j) = nsize
         endif
      enddo
      istart = istarter(rank)
      iend   = istart + ilength(rank) - 1
! Send out how to break up loops over 1 to ncolA
      nsize = ncolA/(numtasks)
      nextra= ncolA - (nsize*(numtasks))
      allocate (jstarter(0:numtasks-1),jlength(0:numtasks-1))
      jstarter(:) = 1
      jlength(:)  = 0
      j = 0
      if (nextra .gt. j) then
         jlength(j) = nsize + 1
      else
         jlength(j) = nsize
      endif
      do j = 1,numtasks-1
         jstarter(j) = jstarter(j-1) + jlength(j-1)
         if (nextra .gt. j) then
            jlength(j) = nsize + 1
         else
            jlength(j) = nsize
         endif
      enddo
      jstart = jstarter(rank)
      jend = jstart + jlength(rank) - 1

!
! ---- allocate setup arrays
!
      ALLOCATE(f(ncolA))
      ALLOCATE(bnn(nrowA))
      ALLOCATE(max_xx(ncolA))
      ALLOCATE(min_xx(ncolA))
      ALLOCATE(aty(0:ncolA))
      ALLOCATE(xx(ncolA))
      ALLOCATE(xp(ncolA))
      ALLOCATE(yy(nrowA))
      ALLOCATE(dxx(ncolA))
      ALLOCATE(dyy(0:nrowA))
      ALLOCATE(yy_can_be_neg(nrowA))
      yy_can_be_neg(:) = .false.
      ALLOCATE(yp(nrowA))
      allocate(shifter(ncola))
      allocate(mapper(0:ncola))

      READ(lu_setup) bnn(:),max_xx(:),min_xx(:)
      allocate(signs(ncolA))
      READ(lu_setup) signs(:)
      READ(lu_setup) shifter(:)
      mapper(0) = 0
      READ(lu_setup) mapper(1:ncolA)
      READ(lu_setup) yy_can_be_neg(1:nrowA)
      READ(lu_setup) f

      allocate ( A_loci(istart:iend)) ! A           break up nrowa
      allocate (At_locj(jstart:jend)) ! A-transpose break up ncola

      do j = jstart,jend
        At_locj(j)%used = 0
        At_locj(j)%length = 16
        allocate(At_locj(j)%vec(At_locj(j)%length))
        allocate(At_locj(j)%ind(At_locj(j)%length))
      enddo
      allocate(ind(1:ncola))
      allocate(vec(1:ncola))
      allocate(temp_vec(1:nrowa))
      allocate(temp_ind(1:nrowa))

      count = 0
      do i=1,nrowA
        read(lu_setup) used
        if(used .ne. 0) then
         read(lu_setup) ind(1:used)
         read(lu_setup) vec(1:used)
        endif
        count = count+used
        if( i .ge. istart .and. i .le. iend) then
          A_loci(i)%used = used
          A_loci(i)%length = used
          allocate(A_loci(i)%vec(used))
          allocate(A_loci(i)%ind(used))
          A_loci(i)%vec(1:used) = vec(1:used)
          A_loci(i)%ind(1:used) = ind(1:used)
        endif
        if (used .eq. 0) goto 7221
        if(ind(1) .le. jend .and. ind(used) .ge. jstart) then
         do j = 1, used
           jj = ind(j)
           if (jj .ge. jstart .and. jj .le. jend) then !TODO binary search
            if (At_locj(jj)%length .eq. At_locj(jj)%used) then ! growing
              loc_len = At_locj(jj)%length
              temp_vec(1:loc_len) = At_locj(jj)%vec(1:loc_len)
              temp_ind(1:loc_len) = At_locj(jj)%ind(1:loc_len)
              deallocate(At_locj(jj)%vec,At_locj(jj)%ind)
              loc_len_old = loc_len
              loc_len = loc_len*2
              if (loc_len .gt. nrowA) loc_len=nrowA
              allocate(At_locj(jj)%vec(loc_len),At_locj(jj)%ind(loc_len))
              At_locj(jj)%vec(1:loc_len_old) = temp_vec(1:loc_len_old)
              At_locj(jj)%ind(1:loc_len_old) = temp_ind(1:loc_len_old)
              At_locj(jj)%length = loc_len
            endif
            loc_int = At_locj(jj)%used + 1
            At_locj(jj)%used = loc_int
            At_locj(jj)%vec(loc_int) = vec(j)
            At_locj(jj)%ind(loc_int) = i
           else if (jj .gt. jend) then
            goto 7221 ! done
           endif
         enddo
 7221   endif
      enddo

      CLOSE(UNIT=lu_setup)
      if (rank .ne. 0) call worker()

      print*,'number of non-zero elements:',count


! --- maximize fval = f(transpose)*x
! --- subjext to A*x <= b  and x >=0 and x < ub
! --- dual
! --- minimize dval = b(transpose)*y
! --- subject to A(transpose)*y >= f and y >=0
!
! --- initialize the solution and derivative vectors
!
      memwrds = 0
      memGB = (8.0d0*memwrds)/(1024**3)
      PRINT 85, memGB,'GB memory needed for array storage'
      bigint1 = ncolA
      bigint2 = nrow0A
      memwrds = bigint1*bigint2
      memGB = (8.0d0*memwrds)/(1024**3)
      print 85, memGB,'GB memory needed for entire A Matrix'
 85   format(F5.1,A)

! Read in optional XX input
      dyy_max = 0.0d0
      call initialize_xx()

      fmax = 0.0d0
      yp(:) = yy(:)
      xp(:) = xx(:)
      dyy = 0.0d0
      iter_close = -1
      iter = 0
      gain_cntrl = 0
      dxx_max = 0.0d0
      dxx = 0.0d0
      fval = 0.0d0
      dval = 0.0d0
      Aty = 0.0d0

      itfmax = 0
      sp_err = 0.0d0
      sp_err_max = 0.0d0
      sp_err_min = 0.0d0
      sp_err_it = 0
      sp_err_it0 = 0
      delerr = 0.0d0
      sp_err_prev = 0.0d0
      delta_fval = 0.0d0
      fval0 = 0.0d0
      delta_dval = 0.0d0
      dval0 = 0.0d0

      t0 = ETIME(times)

      allocate(diis_dxx(ncola,diis_stepx),diis_dyy(nrowa,diis_stepy),diis_x(ncola,diis_stepx),diis_y(nrowa,diis_stepy)) ! Old stuff
      allocate(diis_xa(diis_stepx+1,diis_stepx+1),diis_ya(diis_stepy+1,diis_stepy+1),diis_tmpx(diis_stepx+1,diis_stepx+1),
     &           diis_tmpy(diis_stepy+1,diis_stepy+1)) ! Dot products of old stuff
      allocate(diis_yw(diis_stepy+1),diis_xw(diis_stepx+1)) ! Solution vectors
      diis_yw(:) = 0
      diis_yw(diis_stepy) = 1
      diis_xw(:) = 0
      diis_xw(diis_stepx) = 1
      diis_dxx(:,:) = 0
      diis_dyy(:,:) = 0
      diis_x(:,:) = 0
      diis_y(:,:) = 0
      diis_xa(:,:) = 0
      diis_xa(diis_stepx+1,:) = 1
      diis_xa(:,diis_stepx+1) = 1
      diis_xa(diis_stepx+1,diis_stepx+1) = 0
      diis_ya(:,:) = 0
      diis_ya(diis_stepy+1,:) = 1
      diis_ya(:,diis_stepy+1) = 1
      diis_ya(diis_stepy+1,diis_stepy+1) = 0

! --- solve the nonlinear dynamical system ************************ main iteration loop *********

      PRINT *, 'optimization start'
      write(*,1)'iteration        fval           dval     sp_err         kx         ky        dt     dx_max     dy_max       gain1'
1     format(a)
      gain=1
      gain0=1

      DO ! This is the outer loop over iterations.  We use goto to break out

      IF (dyy_max .eq. 0.0d0) THEN
         iihi = 1                       ! don't run inner loop during initial phase when no y contribution
         gain = gain0*(2**(gain_cntrl))      ! increase gain while no y contribution
         gain_cntrl = gain_cntrl + 1
         if (gain_cntrl .gt. 20) gain_cntrl = 20
      ELSE
         iihi = max_ii
         gain = gain0
         gain_cntrl = 0
      END IF
      IF (iter + iihi .ge. itprofile(nprofile)) THEN
         max_ii = iiprofile(nprofile)
         iihi = max_ii
         gain0 = gprofile(nprofile)
         gain = gain0
         gain_cntrl = 0
         dt = dtprofile(nprofile)
         PRINT *, 'switching to profile',nprofile,'with gain',gain0,'and max_ii',max_ii,'and dt',dt
         nprofile = nprofile +1
      END IF
      DO ii = 1, iihi
         iter = iter + 1
! DIIS
         if (mod(iter,20000) .ge. 10000) then ! TODO make this input controlled
           kx = 0 ! Flag diis in output
           ky = 0
           dt_kx = (0.01)*kx
           dt_ky = (0.01)*ky ! Fake dt
           gain_cntrl = 0
         else
           kx = gain/ii
           ky = 1.d0/kx
           dt_kx = dt*kx
           dt_ky = dt*ky
         endif

! --- setup and solve primal solution xx

!        dx = dt   *(f-A'*(y+(k1/dt)dy))  original
!        dx = dt*k1*(f-A'*(y+dy))         adaptive

         IF (dyy_max .GT. 0.0d0) THEN
            style = 1
            call mpi_bcast(style,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

            call calculate_Aty()

! Pass in bogus zero element, so that jstarter is actually displacements
            call mpi_gatherv(MPI_IN_PLACE,-1,-1,Aty(0),jlength,jstarter,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
         ELSE
            Aty(:) = 0.0d0
         END IF


!        dx = max(x+dx,0)-x;   x = x+dx;

         dval_shift = 0.0d0
         DO j = 1, ncolA ! Make dxx
            dxx(j) = dt_kx*(f(j)-Aty(j))
            if ( max_xx(j) .lt. 1.0d15) then
            dval_shift = dval_shift + max_xx(j)*max(0.0d0,f(j)-Aty(j))
            endif
            if ( min_xx(j) .gt. -1.0d15) then
            dval_shift = dval_shift + min_xx(j)*min(0.0d0,f(j)-Aty(j))
            endif

! dval_shift is needed because the x<max_xx constraint does not contribute to dval directly
! In the old days it was in the A matrix and all was well with dval, but convergence was terrible
!           Here's the Math
! (A-transpose' * y ) >= f   ! this is full A matrix including max_xx constraints
! For all i: (At' * y)(i) >= f(i)
! For all i: (At * y)(i) + y(i+m) >= f(i) ! Split out max_xx constraints
! For all i: y(i+m) >= F(i) - (At * y)(i)
! But, our goal is to minimize b*y, so choose equality (subject to y>=0 constraint)
! For all i: y(i+m) = max(0,F(i) - (At * y)(i))
! Thus the addition to dval is (with b = max_xx for these y's):
! dval_shift = sum(i) [ max(0,F(i) - (At * y)(i)) * max_xx(i) ]
! Note, this ALWAYS obeys At*y >= f constraint, because we made it that way
! Note, dval on iteration one is max_xx*f, ie. the maximum profit, assume no constraints.
         ENDDO
         dval_shift = max(0.0d0,dval_shift) ! Can't shift down

         if (iter.eq.1 .and. read_xx) dxx(:) = 0

         dxx_max = 0.0d0 ! Apply dxx
         DO j = 1, ncolA
            IF (xx(j)+dxx(j) .ge. min_xx(j) .and. xx(j)+dxx(j) .le. max_xx(j) ) THEN
                !
            ELSE IF (xx(j)+dxx(j) .lt. min_xx(j)) then
                dxx(j) = min_xx(j) - xx(j)
            ELSE
                dxx(j) = max_xx(j) - xx(j)
            END If
            dxx_max = MAX(dxx_max, ABS(dxx(j)))       !ste
         END DO
         do i = 1,diis_stepx-1
            diis_dxx(:,i) = diis_dxx(:,i+1)
            diis_x(:,i)   = diis_x(:,i+1)
         enddo
         diis_dxx(:,diis_stepx) = dxx(:)
         diis_x(:,diis_stepx)   = xx(:)

! Do DIIS
         do i = 1, diis_stepx-1
            diis_xa(i,:) = diis_xa(i+1,:)
         enddo
         do i = 1, diis_stepx-1
            diis_xa(:,i) = diis_xa(:,i+1)
         enddo
         do i = 1, diis_stepx
            diis_xa(i,diis_stepx) = dot_product(diis_dxx(:,diis_stepx),diis_dxx(:,i))
            diis_xa(diis_stepx,i) = diis_xa(i,diis_stepx)
         enddo 

! NOT DIIS
         if (mod(iter,20000) .lt. 10000+diis_stepx) then
           diis_xw(:) = 0
           diis_xw(diis_stepx) = 1
! DIIS
         else
           diis_xw(1:diis_stepx) = 0
           diis_xw(diis_stepx+1) = 1
           diis_tmpx(:,:) = diis_xa(:,:)
           Call DSYSV('U', diis_stepx+1, 1, diis_tmpx,diis_stepx+1,diis_IPIV,diis_xw,diis_stepx+1,DIIS_WORK,diis_LWORK,ierr)
           if (ierr .ne. 0 .and. iter .ge. diis_stepx+1) then
             diis_xw(:) = 0
             diis_xw(diis_stepx) = 1
           else if ((maxval(diis_xw(:)) .gt. 10) .or. (minval(diis_xw(:)) .lt. -10) .or. (iter .le. diis_stepx+1)) then
             diis_xw(:) = 0
             diis_xw(diis_stepx) = 1
           endif
         endif

         dxx(:) = 0
         xx(:) = 0
         do i = 1, diis_stepx
           dxx(:) = dxx(:) + diis_dxx(:,i)*diis_xw(i)
           xx(:)  = xx(:)  + diis_x(:,i)  *diis_xw(i)
         enddo
         xx(:) = xx(:) + dxx(:)
         DO j = 1, ncolA
           if (xx(j) .gt. max_xx(j)) xx(j) = max_xx(j)
           if (xx(j) .lt. min_xx(j)) xx(j) = min_xx(j)
         enddo

         DO j = 1, ncolA
           xp(j) = xx(j) + dxx(j)
           if (xp(j) .lt. min_xx(j)) xp(j) = min_xx(j)
           if (xp(j) .gt. max_xx(j)) xp(j) = max_xx(j)
         END DO

!        dy = dt*   (-b+A*(x+(k2/dt)dx))                original
!        dy = dt*k2*(-b+A*(x+dx))                       adaptive

         style = 2
         call mpi_bcast(style,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call mpi_bcast(dt_ky,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         call mpi_bcast(xp,ncolA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

         call calculate_dyy(xp)

! Pass in bogus zero element, so that istarter is actually displacements
         call mpi_gatherv(MPI_IN_PLACE,-1,-1,dyy(0),ilength,istarter,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)

!        dy = max(y+dy,0)-y; y = y + dy

         dyy_max = 0.0d0                                        !ste
         DO j = 1, nrowA
            IF (yy(j)+dyy(j) .ge. 0.0d0 .or. yy_can_be_neg(j)) THEN
                !
            ELSE
                dyy(j) = -yy(j)
            END IF
            dyy_max = MAX(dyy_max, ABS(dyy(j)))       !ste
         END DO

         do i = 1,diis_stepy-1
            diis_dyy(:,i) = diis_dyy(:,i+1)
            diis_y(:,i)   = diis_y(:,i+1)
         enddo
         diis_dyy(:,diis_stepy) = dyy(1:nrowa)
         diis_y(:,diis_stepy)   = yy(:)

! Doing DIIS
         do i = 1, diis_stepy-1
            diis_ya(i,:) = diis_ya(i+1,:)
         enddo
         do i = 1, diis_stepy-1
            diis_ya(:,i) = diis_ya(:,i+1)
         enddo
         do i = 1, diis_stepy
            diis_ya(i,diis_stepy) = dot_product(diis_dyy(:,diis_stepy),diis_dyy(:,i))
            diis_ya(diis_stepy,i) = diis_ya(i,diis_stepy)
         enddo

         if (mod(iter,20000) .lt. 10000+diis_stepy) then
! Not DIIS
           diis_yw(:) = 0
           diis_yw(diis_stepy) = 1
         else
! DIIS
           diis_yw(1:diis_stepy) = 0
           diis_yw(diis_stepy+1) = 1
           diis_tmpy(:,:) = diis_ya(:,:)
           Call DSYSV('U', diis_stepy+1, 1, diis_tmpy,diis_stepy+1,diis_IPIV,diis_yw,diis_stepy+1,DIIS_WORK,diis_LWORK,ierr)
           if (ierr .ne. 0 .and. iter .ge. diis_stepy+1) then
             diis_yw(:) = 0
             diis_yw(diis_stepy) = 1
           else if ((maxval(diis_yw(:)) .gt. 10) .or. (minval(diis_yw(:)) .lt. -10) .or. (iter .le. diis_stepy+1)) then
             diis_yw(:) = 0
             diis_yw(diis_stepy) = 1
           endif
         endif

         dyy(:) = 0
         yy(:) = 0
         do i = 1, diis_stepy
           dyy(1:nrowa) = dyy(1:nrowa) + diis_dyy(:,i)*diis_yw(i)
           yy(:)  = yy(:)  + diis_y(:,i)  *diis_yw(i)
         enddo
         yy(:) = yy(:) + dyy(1:nrowa)
         DO j = 1, nrowA
           if (yy(j) .lt. 0 .and. .not. yy_can_be_neg(j)) yy(j) = 0
         enddo

         DO j = 1, nrowA
           yp(j) = yy(j) + dyy(j)
           if(yp(j) .lt. 0.0d0 .and. .not. yy_can_be_neg(j)) yp(j) = 0.0d0
         END DO

!        Calculate objective function
         fval = 0.0d0                                                !ste
         DO j = 1, ncolA                                     !ste
            fval = fval + f(j) * xx(j)                       !ste
         END DO                                              !ste

         dval = 0.0d0
         DO j = 1, nrowA
            dval = dval + abs(bnn(j) * yy(j)) ! Abs because YY can be negative for equalities
         END DO
         dval = dval + dval_shift

         sp_err = ABS(fval - dval)
         delerr = sp_err - sp_err_prev
         sp_err_prev = sp_err
         delta_fval = fval - fval0                           !ste
         fval0 = fval                                        !ste
         delta_dval = dval - dval0
         dval0 = dval
         IF ( sp_err .gt. sp_err_max) THEN
            sp_err_max = sp_err
            sp_err_it = iter
            sp_err_min = sp_err_max  !ste ignore sp_err_min before max error reached
         END IF
         IF ( sp_err .lt. sp_err_min) THEN
            sp_err_min = sp_err
            sp_err_it0 = iter
            fval_close = fval
            dval_close = dval
         END IF
         PRINT 900, iter,fval,dval,sp_err,kx,ky,dt,
     *        dxx_max,dyy_max, gain
  900    FORMAT(I7,2F15.2,1P8D11.2)                          !ste

        if( mod(iter,1000) .eq. 0 .or. iter .eq. 1 .or. iter .eq. num_iter) then
         i = 0
         j = 0
         k = 0
         do jj = 1,ncolA
           if(min_xx(jj) .eq. 0.0d0 .and. max_xx(jj) .eq. 0.0d0) then
             ! Not really active variable
           else if (xx(jj) .le. min_xx(jj)) then
             i = i+1
           else if (xx(jj) .ge. max_xx(jj)) then
             j = j+1
           endif
         enddo
         do jj = 1,nrowA
           if (yy(jj) .ne. 0.0d0) k = k +1
         enddo
         l2norm4 = 0
         do jj = 1,nrowA
           if (yy_can_be_neg(jj) .or. diis_dyy(jj,diis_stepy) .gt. 0.0d0 ) then
              l2norm4 = l2norm4 + diis_dyy(jj,diis_stepy)**2
           endif
         enddo
         l2norm4 = sqrt(l2norm4)/dt_ky
         print 49, iter,'   MIN_XX=',i,'   MAX_XX=',j,'   active YY=',k,'  total=',i+j+k,' norm2=',l2norm4

 49      format(i,a,i,a,i,a,i,a,i,a,f)

        if( iter .eq. num_iter) then
         write(filename_out,'(a,a,i0,a)')trim(filename),'.',iter,'.xout'
         open(unit=lu_xx,FILE=trim(filename_out),err=631,status='new')
         do i = 1, ncolA
           do j = mapper(i-1) + 1, mapper(i) - 1
             write(lu_xx,*) '0.0d0             dummy'
           enddo
           write(lu_xx,*) (xx(i)-shifter(i))*signs(i)
         enddo
         close(unit=lu_xx)


         write(filename_out,'(a,a,i0,a)')trim(filename),'.',iter,'.yout'
         open(unit=lu_xx,FILE=trim(filename_out),err=631,status='new')
         do i = 1, nrowA
           write(lu_xx,*) yy(i)
         enddo
         close(unit=lu_xx)
 631     continue
        ENDIF

        endif

            IF (fval .gt. fmax) THEN                            !ste
                fmax = fval                                     !ste
                itfmax = iter                                   !ste
            END IF                                              !ste
            IF (abs(fval) .gt. 1.d+15 .or. abs(dval-dval_shift) .gt.1.d+15) THEN                             !ste
                t1 = ETIME(times)
                PRINT *, 'Nonlinear Dynamical System solver time: ', t1-t0
                PRINT *, 'numerical instability encountered; lower dt'
                call abort()
            END IF                                              !ste
            IF (delerr .le. 0 .and. sp_err .le. err_max .and. delta_fval .ge. 0.0d0) THEN
               IF (iter-1 .eq. iter_close .and. fmax-fval .lt. err_max) THEN
                  PRINT *,'converged at iteration',iter,' to',err_max
                  num_iter = iter
                  GOTO 100
               END IF
               iter_close = iter
            END IF
          IF (iter .ge. num_iter ) GOTO 99
         END DO ! inner loop on ii
      END DO ! Outer infinite loop
   99 PRINT *,'NOT converged to tolerance',err_max,' at iteration',iter
  100 CONTINUE
      t1 = ETIME(times)
      PRINT *, 'Nonlinear Dynamical System solver time: ', t1-t0

      PRINT *
      PRINT 32, 'fval=', fval, '   dval=',dval,'   primal-dual error=',sp_err,' after ', iter, ' iterations'
 32   FORMAT (a,f,a,f,a,f,a,i,a)
      PRINT *, 'maximum fval=', fmax, ' at iteration',itfmax
      PRINT *, 'sp_err_max =',sp_err_max,' at iteration',sp_err_it
      IF (sp_err_max .ne. sp_err_min) PRINT *, 'sp_err_min ='
     *   ,sp_err_min,'at iteration',sp_err_it0,'fval=',fval_close
     *   ,'dval=',dval_close

      call abort()

  911 CONTINUE
      PRINT *, 'Error opening unit',lu_setup,'with file',filename_setup
      call abort()

      contains

      subroutine worker()
      implicit none
      integer i,j,ii,jj,k

 2010  continue

      call mpi_bcast(style,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if (style .eq. 1) then

         call calculate_Aty()
         call mpi_gatherv(Aty(jstart),jend-jstart+1,MPI_DOUBLE_PRECISION,-1.0d0,-1,-1,-1,0,mpi_comm_world,ierr)

      else if (style .eq. 2) then

         call mpi_bcast(dt_ky,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         call mpi_bcast(xp,ncolA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         call calculate_dyy(xp)
         call mpi_gatherv(dyy(istart),iend-istart+1,MPI_DOUBLE_PRECISION,-1.0d0,-1,-1,-1,0,mpi_comm_world,ierr)

      else if (style .eq. -1) then

         call mpi_finalize(ierr)
         stop

      else

         print*,'Rank ',rank ,' Got style ', style
         stop

      end if

      goto 2010

      end subroutine

      SUBROUTINE ABORT
      IMPLICIT NONE

      include 'mpif.h'
      INTEGER ierr, style
      style = -1
      call MPI_bcast(style,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_finalize(ierr)

      STOP

      END SUBROUTINE

      Subroutine Calculate_Aty()
      implicit none
      integer j,k,ierr
      real(kind=dp) temp2
      real(kind=dp) Aty_pos(jstart:jend),Aty_neg(jstart:jend)

      call mpi_bcast(yp,nrowA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! only send part

      Aty_pos(jstart:jend) =0.0d0
      Aty_neg(jstart:jend) =0.0d0

      DO k = jstart,jend
         DO j = 1, At_locj(k)%used
            temp2 = At_locj(k)%vec(j) * yp(At_locj(k)%ind(j))
            if(temp2 .ge. 0) then 
              Aty_pos(k)=Aty_pos(k)+temp2
            else
              Aty_neg(k)=Aty_neg(k)+temp2
            endif
          END DO
      END DO
      Aty(jstart:jend) = Aty_pos(jstart:jend) + Aty_neg(jstart:jend) ! More stable dot-product

      end subroutine


      subroutine calculate_dyy(xx_in)
      implicit none
      REAL(KIND=DP), DIMENSION(ncolA),intent(IN) :: xx_in
      REAL(KIND=DP) temp2
      real(kind=dp) dyy_neg(istart:iend),dyy_pos(istart:iend)
      integer j,k

      dyy_pos(istart:iend) = 0
      dyy_neg(istart:iend) = 0
      DO j = istart,iend
         DO k = 1, A_loci(j)%used
          temp2 = A_loci(j)%vec(k) * xx_in(A_loci(j)%ind(k))
          if(temp2 .ge. 0) then
             dyy_pos(j) = dyy_pos(j) + temp2
          else
             dyy_neg(j) = dyy_neg(j) + temp2
          endif
        END DO
      END DO
      dyy(istart:iend) = dt_ky*((dyy_neg(istart:iend)+dyy_pos(istart:iend))-bnn(istart:iend)) ! More stable dot-product - split pos/neg


      end subroutine

      subroutine print_parallel_info()

      character*(MPI_MAX_PROCESSOR_NAME) name1, name2
      integer j
#ifdef _OPENMP
      character*1024  evalue
#endif
 
      name1 = ''
      call MPI_GET_PROCESSOR_NAME(name1,namelen,ierr)
      if(rank .ne. 0) then
         call mpi_send(name1,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,0,987,MPI_COMM_WORLD,ierr)
         return
      else
         print*,'MPI spawned ', numtasks, ' tasks'
         numnodes = 1
         do j = 1, numtasks-1 ! Assumes in order assignment of cores
            name2=''
            call mpi_recv(name2,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,j,987,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
            if(name1(1:MPI_MAX_PROCESSOR_NAME) .ne. name2(1:MPI_MAX_PROCESSOR_NAME)) then
               name1=name2
               numnodes = numnodes+1
            endif
         enddo
         print*, 'Running on this many unique nodes: ',numnodes
      endif
#ifdef _OPENMP
! Only node zero of MPI gets here
#error No OpenMP code at this time
      evalue = ''
      call getenv( 'OMP_NUM_THREADS', evalue )
      if (evalue .eq. '') then
        print*, 'Using OpenMP with OMP_NUM_THREADS unset'
      else
        print*, 'Using OpenMP with OMP_NUM_THREADS = ', TRIM(evalue)
      endif
#endif

      end subroutine print_parallel_info

      subroutine read_profile()
      implicit none
      integer jj,j
      real(kind=DP) temp1,temp2

      nprofile = 1
      OPEN(UNIT=lu_input, FILE=trim(filename)//'.profile', STATUS='OLD', FORM='FORMATTED', ERR=009)
      print*,'Reading profile information from ',trim(filename)//'.profile'
      DO  ! First we could how many lines
         READ(lu_input,*,END=008,ERR=008) jj,temp1,j,temp2
         if (jj .eq. 0 .or. j .eq. 0 .or. temp1 .eq. 0.0d0 .or. temp2 .eq. 0) goto 008 ! Blank lines
         nprofile = nprofile + 1
      END DO
  008 continue
      allocate (itprofile(nprofile),gprofile(nprofile),iiprofile(nprofile),dtprofile(nprofile))
      REWIND(lu_input)
      DO j = 1,nprofile-1
         READ(lu_input,*) itprofile(j),dtprofile(j),gprofile(j),iiprofile(j)
         PRINT *, 'profile:',itprofile(j),dtprofile(j),gprofile(j),iiprofile(j)
      END DO
  009 CONTINUE
      IF (.not. allocated(itprofile) ) then
        allocate (itprofile(1))
      else
        close(lu_input)
      endif
      itprofile(nprofile) = num_iter * 2 ! never go off of end
      end subroutine

      subroutine initialize_xx()
      implicit none
      integer i,j

      read_xx = .false.
      yy = 0.0d0
      xx = 0.0d0
      OPEN(UNIT=lu_xx, FILE=trim(filename)//'.xin', STATUS='OLD', ERR=383)
      print*,'Reading initial guess from file '//trim(filename)//'.xin'
      do j = 1, ncolA
        do i = mapper(j-1) + 1, mapper(j) - 1
          read(lu_xx,*)
        enddo
        read(lu_xx,*) xx(j)
      enddo
      read_xx = .true.
      CLOSE(UNIT=lu_xx)

 383  continue

      do j = 1, ncolA ! Intenal signs might not match external data
        xx(j) = xx(j)*signs(j) + shifter(j)
      enddo

      do j = 1, ncolA
        xx(j) = max(min_xx(j),xx(j))
        xx(j) = min(max_xx(j),xx(j))
      enddo


      OPEN(UNIT=lu_xx, FILE=trim(filename)//'.yin', STATUS='OLD', ERR=384)
      print*,'Reading initial guess from file '//trim(filename)//'.yin'
      do j = 1, nrowA
        read(lu_xx,*) yy(j)
      enddo
      dyy_max = 1.0d0
      CLOSE(UNIT=lu_xx)
 384  continue

      end subroutine

      END PROGRAM
