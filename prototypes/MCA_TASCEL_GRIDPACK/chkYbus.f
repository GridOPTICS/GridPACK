! ------------- SUBROUTIEN ------------------

      SUBROUTINE CHKYBUS
      
!     WRITTEN BY YOUSU CHEN, PNNL, 02/15/2013
! --------------- DESCRIPTION ---------------
!     OUTPUT YBUS MATRIX FOR CHECKING
!      
!-------------------------------------------------

      USE BUSMODULE, ONLY:NB
      USE YBUSMODULE, ONLY:IROWYBUS,ICOLYBUS,VALYBUS
      IMPLICIT NONE
      INTEGER I,J,ID

      OPEN(1,FILE='YBUS.CHK',STATUS='UNKNOWN')
      ID = 0
      DO I = 1, NB
        DO J = 1, IROWYBUS(I+1)-IROWYBUS(I)
          ID = ID+1
          WRITE (1,'(I5,1X,I5,2(1X,F12.4))') I, ICOLYBUS(ID),VALYBUS(ID)
        ENDDO
      ENDDO
      
      RETURN
      ENDSUBROUTINE

