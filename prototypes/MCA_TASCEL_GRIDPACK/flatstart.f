!-----------------  SUBROUTINE -------------------
      SUBROUTINE FLATSTART
!      
! WRITTEN BY YOUSU CHEN, PNNL, 03/17/2007
! LAST MODIFIED BY YOUSU CHEN, PNNL, 04/05/2007
! $ID: FLATSTART.F V1.02 2007/04/05
! ------------------------------------------------
! REVISION LOG:
!   1. ADD VGG ARRAY TO CONSIDER GENERATOR STATUS
! 2007/04/05
!
! ----------------- DESCRIPTION ------------------
!
!     CREATES FLAT START V0
!
!-------------------------------------------------

      USE BUSMODULE, ONLY: NB,VA,VM
      USE GENMODULE, ONLY: NG,NGON,GEN_STATUS,VG,GBUS
      USE OTHERMODULE,ONLY: V0,V0SAVE
      USE CONSTANTS
      USE INTRIFUNS
      IMPLICIT NONE

      INTEGER :: I,IND
      REAL(KIND=DP)::THETA
      COMPLEX(KIND=DP),DIMENSION(NG)::VGG

      DO I = 1,NB
         THETA = VA(I) * PI/180.0
         V0(I) = DCMPLX (VM(I)*DCOS(THETA), VM(I)*DSIN(THETA))                  
      END DO
      IND = 0
      DO I = 1, NG
         IF (GEN_STATUS(I) .EQ.1 ) THEN
	    IND = IND + 1
	    VGG(IND) = VG(I)
	 END IF
      END DO

      DO I = 1,NGON
         V0(GBUS(I)) = VGG(I)* V0(GBUS(I))/ABS(V0(GBUS(I)))
      END DO
      V0SAVE=V0
      
      RETURN
      END SUBROUTINE FLATSTART
