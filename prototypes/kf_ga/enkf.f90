	  PROGRAM ENKF
	  implicit none
	  integer(4) kstep,istep,i,j,k,nsteps,kfault,kclear,in,out
      integer(4) nX,nZ,nE,nG,nB,ksize,krecv,info,ninner,i2
      complex(8), dimension(:,:),target,allocatable :: recv_f
      complex(8), dimension(:,:),target,allocatable :: recv_0
      complex(8), dimension(:,:),pointer :: recv
      complex(8), dimension(:,:),allocatable :: vth
      complex(8), dimension(:,:),allocatable :: egen
      real(8),allocatable :: xpri(:,:),xfwd(:,:),dxdt1(:,:),dxdt2(:,:)
      real(8),allocatable :: xpri_a(:,:)
      real(8),allocatable :: xpost(:,:),xvar(:,:)
      real(8),allocatable :: hmeas(:,:),hmeas_a(:,:),meas_err(:,:)
      real(8),allocatable :: ph_t(:,:)
      real(8),allocatable :: qmat(:,:),qinv(:,:)
      real(8),allocatable :: wmat(:,:)
      real(8),allocatable :: zmat(:,:)
      real(8),allocatable :: kgain(:,:)
	  real(8),allocatable :: Pm(:),invXdp(:),invHmac(:),emag(:),D(:)
      real(8),allocatable :: kalman(:,:)
	  real(8) timestep,timestep_internal,oonE,nfact,rfact,rnfact
      real(8) omega0,omegaB,std_dev_2,std_dev	
      real(8) re,im,sth,din,s
      real(8) tfault,tclear
      real(8) randnormal
      external randnormal 
!ccccccccccccccccccccccccccccccccccccccccccc
	  in=36
      out=37
	  open(unit=in,file="data1.txt",status="old")
	  read(in,*)nsteps,tfault,tclear
      read(in,*)nB,nG,nZ,nX,ksize
      allocate(recv_0(nG,nB))
      allocate(recv_f(nG,nB))
      allocate(Pm(nG))
      allocate(D(nG))
      allocate(invHmac(nG))
      allocate(invXdp(nG))
      allocate(emag(nG))
      allocate(xpri(nE,nX))
      allocate(xfwd(nE,nX))
      allocate(dxdt1(nE,nX))
      allocate(dxdt2(nE,nX))
      allocate(vth(nE,nG))
      allocate(xpost(nE,nsteps))
      allocate(xvar(nE,nsteps))
      allocate(hmeas(nE,nZ))
      allocate(hmeas_a(nE,nZ))
      allocate(meas_err(nE,nZ))
      allocate(ph_t(nX,nZ))
      allocate(qmat(nE,nE))
      allocate(qinv(nE,nE))
      allocate(zmat(nE,nZ))
      allocate(wmat(nZ,nZ))
      allocate(kgain(nX,nZ))
      allocate(kalman(nZ,nsteps))
          do i=1,nG
             do j=1,nB
                read(in,*)din
                recv_0(i,j)=dcmplx(din,0.0d0);
             end do
          end do  
          do i=1,nG
             do j=1,nB
                read(in,*)din
                re=dreal(recv_0(i,j))*dcos(din)
                im=dreal(recv_0(i,j))*dsin(din)
                recv_0(i,j)=dcmplx(re,im)
             end do
          end do  
          do i=1,nG
             do j=1,nB
                read(in,*)din
                recv_f(i,j)=dcmplx(din,0.0d0);
             end do
          end do  
          do i=1,nG
             do j=1,nB
                read(in,*)din
                re=dreal(recv_f(i,j))*dcos(din)
                im=dreal(recv_f(i,j))*dsin(din)
                recv_f(i,j)=dcmplx(re,im)
             end do
          end do  
      do i=1,nX
         do j=1,nE
            read(in,*) xpri(j,i)
         end do
      end do  
      do i=1,nG
            read(in,*) din
			invHmac(i)=1.0d0/din 
      end do  
      do i=1,nG
            read(in,*) Pm(i)
      end do  
      do i=1,nG
            read(in,*) D(i)
      end do  
      do i=1,nG
            read(in,*) Emag(i)
      end do  
      do i=1,nG
            read(in,*) din
			invXdp(i)=1.0d0/din 
      end do  
      do j=1,nsteps
         do i=1,nB
            read(in,*) kalman(i,j)
         end do
      end do
      do j=nsteps+1,ksize
         do i=1,nB
            read(in,*) din
         end do
      end do
      do j=1,ksize
         do i=1,nB
            read(in,*) kalman(i+nB,j)
         end do
      end do
      close(in)
	  std_dev=0.1d0
      std_dev_2=0.01d0
      timestep=0.01d0
      timestep_internal=0.01d0      
	  kfault=int(nint(tfault/timestep))
      kclear=int(nint(tclear/timestep))
      oonE=1.0d0/nE
      nfact=1.0d0/(dble(nE)-1.0d0)
      rfact=1.0d0/std_dev_2
      rnfact=rfact*nfact
      omega0=1.0d0
      omegaB=60.0d0*datan(1.0d0)*8.0d0
      ninner=int(nint(timestep/timestep_internal))
      if (ninner.lt.1) ninner=1
      write(*,*)"nsteps= ",nsteps
      write(*,*)"ninner= ",ninner
      write(*,*)"nB    = ",nB
      write(*,*)"nG    = ",nG
      write(*,*)"nE    = ",nE        
	  do kstep=1,nsteps
         nullify(recv)
         if (kstep.gt.kfault .and. kstep.le.kclear) then
            recv=>recv_f
         else 
			recv=>recv_0
         end if
! start inner loop
          do istep=1,ninner                       
			call find_deriv(omega0,omegaB,xpri,recv,emag,invHmac,Pm,D,invXdp,vth,dxdt1,nE,nG,nX)
            do i=1,nX
               do j=1,nE
                  xfwd(j,i)=xpri(j,i)+dxdt1(j,i)*timestep_internal
               end do
            end do
            call find_deriv(omega0,omegaB,xfwd,recv,emag,invHmac,Pm,D,invXdp,vth,dxdt2,nE,nG,nX)
            do i=1,nX
               do j=1,nE
                  xpri(j,i)=xpri(j,i)+0.5d0*(dxdt1(j,i)+dxdt2(j,i))*timestep_internal
               end do
            end do
            
         end do
         
! end inner loop
	     call generate_meas(hmeas,xpri,emag,recv,nE,nB,nG,nZ,nX)
         call matrix_mean(nE,nZ,hmeas,hmeas_a)
		 call matrix_mean(nE,nX,xpri,xpri_a)
         do i=1,nZ
           im=kalman(i,kstep)
           do j=1,nE
              re=randnormal(1.0)
              meas_err(j,i)=im*(1.0d0+std_dev*re)-hmeas(j,i)
           end do
         end do
 		 do j=1,nE
            do i=1,nE
               qmat(i,j)=0.0d0
            end do
            qmat(i,i)=1.0d0
         end do
         do i=1,nZ
            do j=1,nZ
               wmat(j,i)=0.0d0
            end do
            wmat(i,i)=1.0d0
         end do
!!!!!!!!!!!!
         call dgemm('t','n',nX,nZ,nE,nfact,xpri_a,nE,hmeas_a,nE,0.d0,ph_t,nX)
         call dgemm('n','t',nE,nE,nZ,rnfact,hmeas_a,nE,hmeas_a,nE,1.0d0,qmat,nE)
         call dpotrs('u',nE,nE,qmat,nE,qinv,nE,info)
         call dgemm('t','t',nZ,nE,nE,rnfact,hmeas_a,nE,qinv,nE,1.0d0,zmat,nZ)
         call dgemm('t','n',nZ,nZ,nE,-1,zmat,nZ,hmeas_a,nE,1.0d0,wmat,nZ)
         call dgemm('n','n',nX,nZ,nZ,rfact,ph_t,nX,wmat,nZ,0.0d0,kgain,nX)
         call dgemm('n','t',nE,nX,nZ,1.0d0,meas_err,nE,kgain,nX,1.0d0,xpri,nE)
		 do i=1,nX
            s=0.0d0
            im=0.0d0
            do j=1,nE
               s=s+xpri(j,i)
               im=im+xpri(j,i)*xpri(j,i)
            end do
		    s=s*oonE
            im=im*oonE-s*s
            xpost(i,kstep)=s
            xvar(i,kstep)=im
         end do
      end do 
      open(unit=out,file="xpost.dat",status='new')
      do j=1,nE
         do kstep=1,nsteps
           write(out,1500,advance='no')xpost(j,kstep)
         end do
		 write(out,*)
      end do
	  write(out,*)" "
      do j=1,nE
         do kstep=1,nsteps
            write(out,1500,advance='no')xvar(j,kstep)
         end do
         write(out,*)
      end do
      close(out)
	  stop
1500 format(d15.5,1x,'+')
      end program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! find time derivative of state
!ccccccccccccccccccccccccccccccccccccc
      subroutine find_deriv(omega0,omegaB,x,recv,emag,invHmac,Pm,D,invXdp,dxdt,nE,nG,nX)
	  implicit none
      integer(4) nE,nG,nX
      integer(4) i2,i,j
	  complex(8) vth(1:nE,1:nG)
      complex(8) egen(1:nG,1:nE)
      complex(8),intent(in)  :: recv(:,:)
	  real(8),intent(in) :: x(:,:)
      real(8),intent(out):: dxdt(:,:)
      real(8),intent(in) :: emag(:),invXdp(:),D(:),invHmac(:),Pm(:)
      real(8),intent(in) :: omega0,omegaB
      real(8) re,im,vij,theta,s1      
!ccccccccc begin statemnets
      do i=1,nG
         i2=i+i-1
         re=emag(i)
         do j=1,nE
            im=x(j,i2)
            egen(j,i)=dcmplx(re*dcos(im),re*dsin(im))
         end do
      end do
      call zgemm('n','n',nE,nG,nG,dcmplx(1.0d0),egen,nE,recv,nG,dcmplx(0.0d0),vth,nE)
      do i=1,nG
		 i2=i+i
         do j=1,nE
            re=dreal(vth(j,i))
            im=dimag(vth(j,i))
            vij=dsqrt(re*re+im*im)
            theta=datan2(im,re)
            s1=dsin(x(j,i2)-omega0)
            dxdt(j,(i2-1))= omegaB*s1
            re=Pm(i)-D(i)*s1-emag(i)*vij*invXdp(i)*dsin(x(j,i2-1)-theta)
            dxdt(j,(i2))=omega0*invHmac(i)*0.5d0*re
         end do
      end do
      return
      end subroutine
!ccccccccccccccccccccccccccccccccccccc
      function randnormal(sig)
      implicit none
      real(8) randnormal
  	  real(8) x,y,r,sig
      randnormal=0.0d0
      end function
!Ccccccccccccccccccccccc       
	  subroutine matrix_mean(nr,nc, mat, mat_a)
      implicit none
      integer(4),intent(in):: nr,nc
      real(8),intent(in) :: mat(:,:)
      real(8),intent(out):: mat_a(:,:)
      integer(4) i,j
      real(8) s
      do i=1,nc
		 s=0.0d0
         do j=1,nr
            s=s+mat(j,i) 
         end do
         s=s/dble(nr)
         do j=1,nr
            mat_a(j,i)=mat(j,i)-s
         end do
      end do
      return
      end subroutine
!cccccccccccccccccccccccccccc
      subroutine make_unit_matrix(nr,nc,mat)
      implicit none
      integer(4),intent(in) :: nr,nc
      real(8),intent(out) :: mat(:,:)
      integer(4) i,j
      do i=1,nc
         do j=1,nr
            mat(j,i)=0.0d0 
         end do
         mat(i,i)=1.0d0
      end do
      return
      end subroutine
!Cccccccccccccccccccccccccccccccccccc
	  subroutine generate_meas(hmeas,x,emag,recv,nE,nB,nG,nZ,nX)
	  implicit none
      integer(4),intent(in) :: nE,nG,nX,nB,nZ
      integer(4) i2,i,j
	  complex(8) vth(1:nE,1:nB)
      complex(8) egen(1:nG,1:nE)
      complex(8),intent(in)  :: recv(1:nG,1:nB)
	  real(8),intent(in) :: x(1:nE,1:nX)
      real(8),intent(out):: hmeas(1:nE,1:nZ)
      real(8),intent(in):: emag(1:nG)
      real(8) re,im 
!ccccccccc begin statemnets
      do i=1,nG
         i2=i+i-1
         re=emag(i)
         do j=1,nE
            im=x(j,i2)
            egen(j,i)=dcmplx(re*dcos(im),re*dsin(im))
         end do
      end do
      call zgemm('n','n',nE,nB,nG,dcmplx(1.0d0),egen,nE,recv,nG,dcmplx(0.0d0),vth,nE)
      do i=1,nB
		 i2=i+i
         do j=1,nE
            re=dreal(vth(j,i))
            im=dimag(vth(j,i))
			hmeas(j,i2-1)=dsqrt(re*re+im*im)
            hmeas(j,i2)=datan2(im,re)
         end do
      end do
      return
      end subroutine
!cccccccccccccccccccccccc
