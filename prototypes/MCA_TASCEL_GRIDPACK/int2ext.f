! ------------ SUBROUTINE -------------------      
      SUBROUTINE INT2EXT
!     WRITTEN BY YOUSU CHEN 03/23/2007      
!     LAST MODIFIED BY YOUSU CHEN 03/23/2007
!     $ID V1.10 2007/03/23
!
!------------  DESSCRIPTION --------------------
!      CONVERT INTERNAL BUS/BRANCH NUMBER TO ORIGINAL EXTERNAL NUMBER
!

      USE ALLMODULE
      IMPLICIT NONE
      ! LOCAL ARGUMENTS
      INTEGER :: I


      DO I = 1,NB
            BUS_I(I) = I2E(BUS_I(I))
      END DO

      DO I = 1,NG
            GEN_BUS(I) = I2E(GEN_BUS(I))
      END DO

      DO I = 1,NBRCH
            F_BUS(I) = I2E(F_BUS(I))
            T_BUS(I) = I2E(T_BUS(I))
      END DO
                  
      END
