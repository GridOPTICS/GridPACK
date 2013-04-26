!------------ SUBROUTINE ----------------------
      SUBROUTINE LIMIT
!
!     WRITTEN BY YOUSU CHEN, PNNL, ON 7/29/2007
!     LAST MODIFIED BY YOUSU CHEN, PNNL, ON 7/29/2007
!     $ID: LIMIT.F V1.01: 2007/07/29
!---------------------------------------------
!
!    FIND THE FACTORS FOR VIORATIONS, TO MAKE SURE THAT THERE IS 
!    NOT VIORATION FOR BASE CASE.
!
! ---------------------------------------------
      USE ALLMODULE

      IMPLICIT NONE

      INTEGER :: I

      VMHIGH=1.1
      VMLOW=0.9

!
!     UPDATE THE CRITERI TO MAKE SURE THERE IS NOT VIORATION FOR BASE CASE
!

      DO I = 1,NB
         IF (VM(I) .GT. 1.1) VMHIGH(I)=VM(I)*1.1
         IF (VM(I) .LT. 0.9) VMLOW(I)=VM(I)*0.9
      ENDDO
    
      MVAHIGH = RATE_A
      DO I = 1, NBRCH
          MVALIM=SQRT(3.0)*ABS(DCMPLX(PF(I),QF(I)))/VM(F_BUS(I)) 
     	  IF(MVALIM .GT. RATE_A(I)) MVAHIGH(I)=1.01*MVALIM
      END DO	  

      END SUBROUTINE LIMIT
