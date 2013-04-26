      SUBROUTINE CHKINPUT


!     WRITTEN BY YOUSU CHEN, PNNL, 06/15/20R07
!
!
! --------------- DESCRIPTION ---------------
!     OUTPUT INPUT MODEL FOR CHECKING
!      
!-------------------------------------------------

      USE BUSMODULE
      USE GENMODULE
      USE BRCHMODULE
      IMPLICIT NONE
      
      INTEGER :: I

      OPEN(1,FILE='FILE.CHK',STATUS='UNKNOWN')
      WRITE (1,'(A12,F7.3)') 'BASEMVA =',BASEMVA
      WRITE(1,*)''
      WRITE (1,*)'%% BUS DATA'
      WRITE (1,*)' BUS_I  TYPE   PD     QD     GS     BS  AREA  VM 
     &   VA   BASEKV ZONE VMAX  VMIN'
      DO I = 1,NB 
          WRITE(1,FMT=101), BUS_I(I),BUS_TYPE(I),PD(I),
     &    QD(I),GS(I),BS(I),BUS_AREA(I),VM(I),
     &    VA(I),BASE_KV(I),ZONE(I),VMAX(I),
     &    VMIN(I)
      END DO
      WRITE(1,*)''
      WRITE(1,*)'%% GENERATOR DATA'
      WRITE(1,*)'  BUS     PG       QG     QMAX     QMIN      VG 
     &    MBASE  STATUS  PMAX    PMIN' 
      DO I = 1,NG
          WRITE(1,FMT=102), GEN_BUS(I),PG(I),QG(I),
     &    QMAX(I),QMIN(I),VG(I),MBASE(I),
     &    GEN_STATUS(I),PMAX(I),PMIN(I)
      END DO
      WRITE(1,*)''
      WRITE(1,*)'%% BRANCH DATA'
      WRITE(1,'(A97)')' FBUS  TBUS      R          X         B      
     &RATEA     RATEB    RATEC    RATIO   ANGLE   STATUS' 
      DO I = 1,NBRCH
          WRITE(1,FMT=103), F_BUS(I),T_BUS(I),
     &        BR_R(I),BR_X(I),BR_B(I),
     &    RATE_A(I),RATE_B(I),RATE_C(I),
     &    TAP(I),SHIFT(I),BR_STATUS(I)
      END DO
  101 FORMAT(I5,6X,I1,1X,F10.2,1X,3(F10.2,1X),1X,I4,3(F8.2,1X),I4,
     & 1X,2F8.2)  
  102 FORMAT(I5,6(1X,F8.2),3X,I4,4X,F7.2,1X,F7.2)   
  103 FORMAT(I5,1X,I5,3(1X,F10.5),5(1X,F8.3),3X,I4)  
      CLOSE (1)

      END SUBROUTINE CHKINPUT
