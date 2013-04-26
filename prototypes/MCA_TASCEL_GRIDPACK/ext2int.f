! ------------- SUBROUTINE -------------------------------

      SUBROUTINE EXT2INT

!     WRITTEN BY YOUSU CHEN, PNNL, 03/23/2007
!     LAST MODIFIED BY YOUSU CHEN, PNNL, 03/29/2007
!     ID: EXT2INT.F V1.02 : 2007/03/29

!     ------------------- DESCRIPTION -----------------
!    CONVERT EXTERNAL NODE NUMBER TO CONSECUTIVE INTERNAL
!    BUS NUMBERS WHICH START AT 1
!
      USE BUSMODULE,ONLY: BUS_I,NB
      USE GENMODULE,ONLY: GEN_BUS,NG
      USE BRCHMODULE,ONLY: F_BUS,T_BUS,NBRCH
      USE OTHERMODULE,ONLY: E2I,I2E
      IMPLICIT NONE

      ! LOCAL ARGUMENTS
      INTEGER :: I

      DO I = 1, NB
         I2E(I) = BUS_I(I)
         E2I(I2E(I)) = I
         BUS_I(I) = E2I(BUS_I(I))
      END DO

      DO I = 1,NG
            GEN_BUS(I) = E2I(GEN_BUS(I))
      END DO

      DO I = 1,NBRCH
            F_BUS(I) = E2I(F_BUS(I))
            T_BUS(I) = E2I(T_BUS(I))
      END DO


      RETURN
      END
