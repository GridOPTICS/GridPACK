      SUBROUTINE REPORTVIOLATION
!      SUBROUTINE REPORTVIOLATION(IDX)

!
!     WRITTEN BY YOUSU CHEN, PNNL, 8/22/2007
!     LAST MODIFIED BY YOUSU CHEN, PNNL, 8/22/2007
!     $ID REPROTVIOLATION.F V1.01: 2007/8/22
!
!
! --------------- DESCRIPTION -----------------------
!
!     REPORT VIOLATIONS FOR BUS VOLTAGE AND BRANCH MVA
!
      USE BUSMODULE 
      USE GENMODULE 
      USE BRCHMODULE 
      USE YBUSMODULE 
      USE OTHERMODULE 
      USE CAMODULE 
      USE MPI

      IMPLICIT NONE
     
!      INTEGER, INTENT(IN):: IDX
      INTEGER :: J,FLAGB, FLAGBR
      LOGICAL :: HERE
      REAL(KIND=DP) :: MEAN, PERC, TEMP
      REAL(KIND=DP) :: SREAL,SIMAG 
      COMPLEX(KIND=DP) :: S_POWER



      FLAGB = 0
      FLAGBR = 0
      IF (CONVERGED .EQ. 1) THEN

!         IF THERE IS VIOLATION FOR BUS DATA, SET FLAGB = 1
!         IF THERE IS VIOLATION FOR BRANCH DATA, SET FLAGBR = 1
!         IF THERE IS NO VIOLATION, NO OUTPUT FILE GENERATED
!
!         FLAG = 1: BRCH
!         FLAG = 2: GENERATOR

      DO J = 1,NB
         IF (VM(J) .GT. VMHIGH(J) .OR. VM(J) .LT. VMLOW(J)) THEN
            FLAGB = 1
            EXIT
         ENDIF
      ENDDO
      DO J = 1,NBRCH
       !  MVALIM=SQRT(3.0)*ABS(DCMPLX(PF(J),QF(J)))/VM(F_BUS(J))
         MVALIM=ABS(DCMPLX(PF(J),QF(J)))
         IF(MVALIM.GT.MVAHIGH(J)) THEN 
!         IF((MVALIM.GT.MVAHIGH(J)).AND. ( J .NE. IDX ) ) THEN 
             FLAGBR = 1
             EXIT
         ENDIF
      ENDDO

      IF (FLAGB.EQ.0 .AND. FLAGBR .EQ.0 ) THEN
          NOK=NOK+1
          CASEFLAG='O'
          GOTO 100
      ELSE
          NVIO=NVIO+1
          CASEFLAG='V'
      ENDIF

      WRITE(1,FMT=105) 'BUS'

      DO J = 1,NB
!        MEAN = (VMHIGH(J)+VMLOW(J))/2
        IF (OUTFLAG .GE. 2) THEN ! REPORT OVER-VOLTAGE VIOLATION ONLY
!           TEMP = VMHIGH(J)
!        ELSE !REPORT ALL OUTPUTS
!           TEMP = 1.0  
!        ENDIF
!        IF (VM(J) >=TEMP) THEN
          IF ((VM(J) >=1.05) .OR. (VM(J) <= 0.95)) THEN
            WRITE(1,FMT=106) BUS_I(J),VM(J)
          ENDIF
        ELSE
          WRITE(1,FMT=106) BUS_I(J),VM(J)
!             PERC = (VM(J)-MEAN)/(VMHIGH(J)-MEAN)*100
        ENDIF
      ENDDO
      WRITE(1,FMT=108) 'BRCH'
      DO J = 1,NBRCH
!         IF (J .NE. IDX) THEN
        IF (BR_STATUS(J) .NE. 0) THEN
             MVALIM=ABS(DCMPLX(PF(J),QF(J)))
!             MVALIM=SQRT(3.0)*ABS(DCMPLX(PF(J),QF(J)))/VM(F_BUS(J))
        !     S_POWER=SQRT(3.0)*(DCMPLX(PF(J),QF(J)))/VM(F_BUS(J))
        !     SREAL = REAL(S_POWER)
        !     SIMAG = AIMAG(S_POWER)
          IF (OUTFLAG .GE. 2) THEN ! REPORT VIOLATION ONLY
             IF(MVALIM .GT. 0.9*MVAHIGH(J)) THEN !90% OVERLOADING
                WRITE(1,FMT=107) J,PF(J),QF(J)
             ENDIF
          ELSE
             WRITE(1,FMT=107) J,PF(J),QF(J)
          ENDIF   
        ENDIF
      ENDDO
      ENDIF

  105 FORMAT(A3) 
  106 FORMAT(I6,2X,F6.4) 
  107 FORMAT(I6,2X,F12.3,1X,F12.3)
  108 FORMAT(A4) 
  100 CONTINUE
      END SUBROUTINE REPORTVIOLATION
