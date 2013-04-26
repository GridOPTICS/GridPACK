!------------ SUBROUTINE ----------------------
      SUBROUTINE PFSOLN
!
!     WRITTEN BY YOUSU CHEN, PNNL, ON 3/16/2007
!     LAST MODIFIED BY YOUSU CHEN, PNNL, ON 6/15/2007
!     $ID: PFSOLN.F V1.04: 2007/06/15
!---------------------------------------------
!     REVISION LOG:
!     2007/03/28
!        1. UPDATED TO F90
!     2007/04/05
!        1. CONSIDERED MORE THAN ONE GENERATOR AT SWING BUS
!        2. ADDED SUMQGG AND SUMQG TO FIND THE NUMBER OF
!           GENERATORS AT SAME BUSES
!     2007/06/15
!        1. IMPLEMENTED DIVIDING QG IN PROPORTION TO THE REACTIVE
!           RANGE OF THE GENERATOR
!
! ------------ DESCRIPTION --------------------
!
!    UPDATES BUS, GENERATOR,BRANCH DATA 
!
! ---------------------------------------------
      USE ALLMODULE
      USE FLAGS
      USE CONSTANTS
      USE SUBS
      IMPLICIT NONE

      INTEGER :: I,J,ID,IND
      INTEGER, DIMENSION(:) :: IROWCG(NGON+1),ICOLCG(NGON),VALCG(NGON)
      INTEGER, DIMENSION(:) :: IROWCGT(NB+1),ICOLCGT(NGON),VALCGT(NGON)
      INTEGER, DIMENSION(:) :: SUMQGG(NGON),SUMQG(NGON)
      INTEGER, ALLOCATABLE,DIMENSION(:) :: RGEN
      REAL(KIND=DP) :: SUMREFGEN
      REAL(KIND=DP),DIMENSION(NB) :: QG_TOT,CMIN,CMAX,EPSTMP,TMP
      REAL(KIND=DP),DIMENSION(NGON) :: TMP2

      VM = ABS(V)
      VA = ATAN2(AIMAG(V),REAL(V))*180.0/PI
!
!     FIND THE NUMBER OF REFERENCE GENERATOR(S)
!     (COULD BE MORE THAN ONE GENERATOR AT SWING BUS)
!

      IND = 0
      DO I = 1,NGON
          IF (GBUS(I) .EQ. SLACK) IND = IND + 1
      END DO
      IF(ALLOCATED(RGEN)) DEALLOCATE(RGEN)
      ALLOCATE(RGEN(IND),STAT=ERROR)
      IF (ERROR /= 0) THEN
          PRINT *, " ------------------ ERROR ----------------------"
          PRINT *, " PROGRAM COULD NOT ALLOCATE SPACE FOR RGEN ARRAYS"
          STOP
      END IF
      IND = 0
      DO I = 1,NGON
          IF (GBUS(I) .EQ. SLACK) THEN
            IND = IND + 1
            RGEN(IND) = I
	    IF (GENOPT .EQ. 1) THEN
	       PMAXSLACK = PPMAX(I)
	    END IF   
          END IF
      END DO	  

!      
!     YBUS(GBUS,:) 
!     
      CALL YBUSGBUS
!     
!      YBUS(GBUS,:) * V
!
      
      CALL  MATMULVEC(IROWG,ICOLG,VALG,V,R,NGON,NB,NZG)
      
!
!      TOTAL INJECTED BUS POWERS :
!      SG = V(GBUS) .* CONJ (YBUS(GBUS,:) * V)
!      
      DO I = 1,NGON
            SG(I) = V(GBUS(I)) * CONJG(R(I))            
      END DO
!
!      UPDATE QG FOR ALL GENERATORS (GEN(ON,QG) = QGG )
!
      QGG = 0
      DO I = 1, NGON
            QGG(I) = AIMAG(SG(I)) * BASEMVA + QD(GBUS(I))
      END DO     
!
!     AT THIS POINT, ANY BUS WITH MORE THAN ONE GENERATOR WILL
!     HAVE THE TOTAL Q DISPATCH FOR THE BUS ASSIGNED TO EACH 
!     GENERATOR. THIS MUST BE SPLIT BETWEEN THEM. 
!     FIRST DO IT EQUALLY, THEN IN PROPORTION TO THE REACTIVE 
!     RANGE OF THE GENERATOR. 
!
!     NOTE: ASSUME IN GENERATOR DATA, SAME GENERATORS ARE POSITIONED TOGETHER
!
      
      IF (NGON .GT. 1) THEN
         IND = 1
         SUMQGG = 1
         DO I = 1, NGON-1
             IF (GBUS(I+1) .EQ. GBUS(I) ) THEN
                 SUMQGG(IND) = SUMQGG(IND) + 1
             ELSE
                 IND = IND + 1
             END IF
         END DO
         SUMQG = 0
         ID = 1
         DO I = 1, IND
             DO J = 1,SUMQGG(I)
                 SUMQG(ID) = SUMQGG(I)
                 ID = ID + 1
             END DO
         END DO     
!
!     DIVIDED QG BY THE NUMBER OF GENERATORS AT THE BUS 
!     TO DISTRIBUTE EQUALLY
!
         QGG = QGG/SUMQG


       
!     FIND QMIN AND QMAX AT THE GENERATOR BUS, (QQMIN, QQMAX)
!
         IF(ALLOCATED(QQMAX)) DEALLOCATE(QQMAX)
         IF(ALLOCATED(QQMIN)) DEALLOCATE(QQMIN)
         ALLOCATE(QQMAX(NGON),QQMIN(NGON),STAT=ERROR)
         IF (ERROR /= 0) THEN
           PRINT *, " ------------------ ERROR ------------------------"
           PRINT *, " PROGRAM COULD NOT ALLOCATE SPACE FOR QQMAX ARRAYS"
           STOP
         END IF
         
         ID = 0
         DO I = 1, NG
             IF (GEN_STATUS(I) .EQ.1 ) THEN
                  ID = ID + 1
                  QQMAX(ID) = QMAX(I)
                  QQMIN(ID) = QMIN(I)
             END IF
         END DO
!        
!     DIVIDE PROPORTIONALLY
!
!      
!        CG = SPARSE([1,NGON]',GBUS,ONES(NGON,1),NGON,NB)
!
         IROWCG(1) = 1
         DO I = 1, NGON
              IROWCG(I+1) = IROWCG(I) + 1
              ICOLCG(I) = GBUS(I)
              VALCG(I) = 1
         END DO
!
!       QG_TOT = CG' * GEN(ON,QG)
!
         CALL GETTRANSPOSEI (IROWCG,ICOLCG,VALCG,NGON,NB,NGON,IROWCGT,
     &       ICOLCGT,VALCGT)
         CALL MATMULVECR(IROWCGT,ICOLCGT,VALCGT,QGG,QG_TOT,NB,NGON,NGON)
!      
          CMIN = 0
	  CMAX = 0
	  DO I = 1,NGON
	     CMIN(GBUS(I)) = CMIN(GBUS(I))+QQMIN(I)
	     CMAX(GBUS(I)) = CMAX(GBUS(I))+QQMAX(I)
	  END DO
	  EPSTMP=EPS
	  TMP = (QG_TOT-CMIN) / (CMAX - CMIN + EPSTMP)
          CALL MATMULVECR(IROWCG,ICOLCG,VALCG,TMP,TMP2,NGON,NB,NGON)
          QGG = QQMIN + TMP2 * (QQMAX-QQMIN)

      ENDIF
                                           
!      
!      UPDATE PG FOR SWING BUS
!      
!      PGG(REFGEN) = REAL(SG(REFGEN)) * BASEMVA + PD(SLACK)
      PGG(RGEN(1)) = REAL(SG(RGEN(1))) * BASEMVA + PD(SLACK)
      
!
!     IF MORE THAN ONE GENERATOR AT THE SWING BUS,  SUBSTRACT OFF
!     WHAT IS GENERATED BY OTHER GENERATORS AT THIS BUS
!
      SUMREFGEN = 0.0
      IF (SIZE(RGEN) .GT. 1) THEN
        DO I = 2, SIZE(RGEN)
              SUMREFGEN = SUMREFGEN + PGG(RGEN(I))
        END DO   
          PGG(RGEN(1)) = PGG(RGEN(1)) - SUMREFGEN
      END IF
!
!      SF = V(F_BUS) .* CONJ(YF * V) * BASEMVA;
!      
      NZYF = IROWYF(NBRCH+1)-1      

      CALL  MATMULVEC(IROWYF,ICOLYF,VALYF,V,R1,NBRCH,NB,NZYF)
      
      DO I = 1,NBRCH
            SF(I) = V(F_BUS(I)) * CONJG(R1(I))*BASEMVA      
      END DO
!
!      ST = V(T_BUS) .* CONJ(YT * V) * BASEMVA;
!      
      NZYT = IROWYT(NBRCH+1)-1

      CALL  MATMULVEC(IROWYT,ICOLYT,VALYT,V,R2,NBRCH,NB,NZYT)
      
      DO I = 1,NBRCH
            ST(I) = V(T_BUS(I)) * CONJG(R2(I))*BASEMVA      
      END DO
      
!
      PF = REAL(SF)
      QF = AIMAG(SF)
      PT = REAL(ST)
      QT = AIMAG(ST)
      DO I = 1, NBRCH
         IF (BR_STATUS(I).EQ.0) THEN
            PF(I) = 0.0
            QF(I) = 0.0
            PT(I) = 0.0
            QT(I) = 0.0
         ENDIF
      ENDDO
      RETURN
      END 
!
!---------- SUBROUTINE YBUSGBUS -----------------
!
      SUBROUTINE YBUSGBUS
!
! ------------- DESCRIPTION ------------
!
!      CREATE SUBMATRIX: YBUS(GBUS,:) 
!

      USE ALLMODULE
      IMPLICIT NONE

      ! LOCAL ARGUMENTS
      INTEGER :: I,J,COUNT,IAA,IAB,TMPI1,TMPI2
!      
!      --- END OF DECLARATION ---
!

      IROWG(1) = 1
      COUNT = 1

      DO I = 1,NGON
            IAA = IROWYBUS(GBUS(I))
            IAB = IROWYBUS(GBUS(I)+1)
            TMPI1 = IAB - IAA   ! THE NUMBER OF NONZERO ELEMENTS IN ROW ROW(I)
            IROWG(I+1) = IROWG(I) + TMPI1 ! NEW IROWG
            
            DO J = COUNT, COUNT+TMPI1-1                  
                  TMPI2 = IROWYBUS(GBUS(I))            
                  ICOLG(J) = ICOLYBUS(TMPI2+J-COUNT) ! NEW ICOLG
                  VALG(J) = VALYBUS(TMPI2+J-COUNT)   ! NEW VALG
            END DO
            COUNT = COUNT + TMPI1    ! COUNT = THE NUMBER OF NONZERO IN NEW MATRIX
      END DO
      NZG = COUNT -1
      
      RETURN
      END SUBROUTINE YBUSGBUS

!
!      --- SUBROUTINE MATMULVECR ---
!
      SUBROUTINE MATMULVECR(IROW,ICOL,A,X,R,N,M,NZ)

!     THE PRODUCT OF INTEGER MATRIX AND A REAL VECTOR.
!-----------------------------------------------------      
!      R = A * X (ALL ARE INTEGERS)
!      IROW, ICOL,A : INPUT MATRIX A
!      X :  VECTOR
!      R : RESULT MATRIX
!      N : THE NUMBER OF ROWS
!      M : THE NUMBER OF COLUMNS
!      NZ: THE NUMBER OF NON-ZERO ELEMENTS

      USE DEFDP      
      IMPLICIT NONE
      
      ! DUMMY ARGUMENTS
      INTEGER, INTENT(IN) :: N,M,NZ
      INTEGER, INTENT(IN), DIMENSION(:) :: IROW(N+1),ICOL(NZ)
      INTEGER, INTENT(IN), DIMENSION(:) :: A(NZ)
      REAL(KIND=DP), INTENT(IN), DIMENSION(:) :: X(M)
      REAL(KIND=DP), INTENT(OUT), DIMENSION(:) :: R(N)

      ! LOCAL ARGUMENTS
      INTEGER :: I, J
      REAL(KIND=DP) :: TEMP
!      
!      --- END OF DECLARATION ---
!
      DO I = 1, N
         TEMP = 0.0
         DO J = IROW(I), IROW(I + 1) - 1
            TEMP = TEMP + A(J) * X(ICOL(J))
         END DO
         R(I) = TEMP
      END DO
      
      RETURN
      END SUBROUTINE MATMULVECR
      
