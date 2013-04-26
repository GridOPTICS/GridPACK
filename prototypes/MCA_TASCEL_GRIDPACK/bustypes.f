!---------------- SUBROUTINE -----------------
       SUBROUTINE BUSTYPES

!    WRITETEN BY YOUSU CHEN, PNNL, 28 MARCH 2007
!    LAST MODIFIED BY YOUSU CHEN, PNNL, 28 MARCH 2007
!    $ID: BUSTYPES.F, V1.02 2007/4/03
!  -------------------------------------------
!    REVISION LOG:
!    2007/4/03:
!        (1) ADD CG MATRIX AND BUS_GEN_STATUS ARRAY TO CONSIDER GENERATOR ARRAYS
!        (2) ALLOCATE MEMORY TO PV AND PQ ARRAYS BASED ON GENERATOR STATUS
!      
!----------- DESCRIPTION ---------------------
!     BUILDS LIST OF SWING BUS, PV BUS AND PQ BUS. (ONLY CONSIDER ONE SWING BUS 
!     AT THIS STEP)
!     IF GEN_STATUS .EQ. 0, TREATED PV BUSES AS PQ BUSES.
!----------------------------------------------------------

      USE BUSMODULE,ONLY: NB,NPV,NPQ,SLACK,BUS_I,BUS_TYPE
      USE GENMODULE,ONLY: NG,GEN_BUS,GEN_STATUS
      USE OTHERMODULE,ONLY:PV,PQ	
      USE FLAGS
      IMPLICIT NONE

      INTEGER :: I,IND
!      INTEGER,DIMENSION(:) :: IROWCGT(NG+1),ICOLCGT(NB),VALCGT(NG)
!      INTEGER,DIMENSION(:) :: IROWCG(NB+1),BUS_GEN_STATUS(NB),X(NG)
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IROWCGT,ICOLCGT,VALCGT
      INTEGER,ALLOCATABLE,DIMENSION(:) :: BUS_GEN_STATUS,X
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IROWCG,ICOLCG,VALCG
      INTEGER,ALLOCATABLE,DIMENSION(:) :: PQQ,PVV 
!
!     CG = SPARSE(GEN(:,GEN_BUS),[1:NG]',GEN(:,GEN_STATUS),NB,NG)
!     OBTAIN CG' FIRST, THEN FIND CG
!     THE NUMBER OF NONZERO ELEMENTS IN CG EQUALS TO IND - 1
!
      IF (ALLOCATED(IROWCGT)) DEALLOCATE(IROWCGT)
      IF (ALLOCATED(ICOLCGT)) DEALLOCATE(ICOLCGT)
      IF (ALLOCATED(VALCGT)) DEALLOCATE(VALCGT)
      ALLOCATE(IROWCGT(NG+1),ICOLCGT(NB),VALCGT(NG),STAT=ERROR)
      IF (ERROR /=0) THEN
           PRINT *, " ------------------ ERROR ------------------------"
           PRINT *, " PROGRAM COULD NOT ALLOCATE SPACE FOR CG' ARRAYS"
           STOP
      END IF
      IROWCGT(1) = 1
      IND = 1
      DO I = 1, NG
         IF (GEN_STATUS(I) .EQ. 1) THEN
             IROWCGT(I+1)=IROWCGT(I)+1
     	     ICOLCGT(IND) = GEN_BUS(I)
	     VALCGT(IND) = GEN_STATUS(I)
	     IND = IND + 1
	 ELSE     
             IROWCGT(I+1)=IROWCGT(I)
         END IF
      END DO
!
!     ALLOCATE MEMORY TO ICOLCG AND VALCG
!

      IF (ALLOCATED(IROWCG)) DEALLOCATE(IROWCG)
      IF (ALLOCATED(ICOLCG)) DEALLOCATE(ICOLCG)
      IF (ALLOCATED(VALCG)) DEALLOCATE(VALCG)
      IF (ALLOCATED(BUS_GEN_STATUS)) DEALLOCATE(BUS_GEN_STATUS)
      IF (ALLOCATED(X)) DEALLOCATE(X)
      ALLOCATE(IROWCG(NB+1),ICOLCG(IND-1),VALCG(IND-1),
     & BUS_GEN_STATUS(NB),X(NG),STAT=ERROR)
      IF (ERROR /=0) THEN
           PRINT *, " ------------------ ERROR ------------------------"
           PRINT *, " PROGRAM COULD NOT ALLOCATE SPACE FOR CG ARRAYS"
           STOP
      END IF
!				       
!     OBTAIN CG      
!
      IND = IND -1
      CALL GETTRANSPOSEI(IROWCGT,ICOLCGT,VALCGT,NG,NB,IND,IROWCG,
     & ICOLCG,VALCG)
!
!     BUS_GEN_STATUS = CG * ONES(NG,1)
! 
      X = 1 
      CALL MATMULVECI(IROWCG,ICOLCG,VALCG,X,BUS_GEN_STATUS,NB,NG,IND)
!
!     ALLOCATE MEMORY TO PQQ,PVV (TEMP ARRAY OF PQ AND PV)
!

      IF (ALLOCATED(PQQ)) DEALLOCATE(PQQ)
      IF (ALLOCATED(PVV)) DEALLOCATE(PVV)
      ALLOCATE(PQQ(NB),PVV(NB),STAT=ERROR)
      IF (ERROR /=0) THEN
          PRINT *," ------------------ ERROR ------------------------"
          PRINT *," PROGRAM COULD NOT ALLOCATE SPACE FOR PQQ/PVV ARRAYS"
          STOP
      END IF
!
!     FIND SLACK, PV, AND PQ
!
      NPV=0
      NPQ=0
      DO I=1,NB
	   IF (BUS_TYPE(I) .EQ. 3 .AND. BUS_GEN_STATUS(I) .NE. 0 ) THEN 
                SLACK = BUS_I(I) ! CONSIDER ONLY ONE SWING BUS
	   ELSEIF (BUS_TYPE(I) .EQ. 2 .AND. BUS_GEN_STATUS(I) 
     &             .NE. 0 ) THEN
                NPV = NPV + 1
	        PVV(NPV) = BUS_I(I) 
	   END IF
	   IF (BUS_TYPE(I) .EQ. 1 .OR. BUS_GEN_STATUS(I) .EQ.0 ) THEN
	        NPQ = NPQ + 1
	        PQQ(NPQ) = BUS_I(I) 
           END IF
      END DO
!      PRINT *,'THE NUMBER OF PV BUSES:',NPV
!      PRINT *,'THE NUMBER OF PQ BUSES:',NPQ
!      PRINT *,'THE SWING BUS NUMBER:',SLACK

!
!     ALLOCATE MEMORY FOR PV AN PQ
!
      IF (ALLOCATED(PV)) DEALLOCATE(PV)
      IF (ALLOCATED(PQ)) DEALLOCATE(PQ)
      ALLOCATE(PV(NPV),PQ(NPQ),STAT=ERROR)
      IF (ERROR /=0) THEN
          PRINT *, " ------------------ ERROR ------------------------"
          PRINT *, " PROGRAM COULD NOT ALLOCATE SPACE FOR PV/PQ ARRAYS"
          STOP
      END IF
!
!     GET PV,PQ
!
      DO I = 1, NPV
         PV(I) = PVV(I)
      END DO
      DO I = 1, NPQ
         PQ(I) = PQQ(I)
      END DO
!
!     DEALLOCATE MEMORY 
!
      DEALLOCATE(IROWCGT)
      DEALLOCATE(ICOLCGT)
      DEALLOCATE(VALCGT)
      DEALLOCATE(IROWCG)
      DEALLOCATE(ICOLCG)
      DEALLOCATE(VALCG)
      DEALLOCATE(BUS_GEN_STATUS)
      DEALLOCATE(X)
      DEALLOCATE(PVV)
      DEALLOCATE(PQQ)
     
      RETURN
      END SUBROUTINE BUSTYPES

!
!      --- SUBROUTINE MATMULVECI ---
!
      SUBROUTINE MATMULVECI(IROW,ICOL,A,X,R,N,M,NZ)

!     THE PRODUCT OF MATRIX AND A VECTOR, ALL ELEMENTS 
!     ARE INTEGERS
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
      INTEGER, INTENT(IN), DIMENSION(:) :: A(NZ),X(M)
      INTEGER, INTENT(OUT), DIMENSION(:) :: R(N)

      ! LOCAL ARGUMENTS
      INTEGER :: I, J, TEMP
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
      END SUBROUTINE MATMULVECI

!
!     -------- SUBROUTINE GETTRANSPOSEI ------------
!

      SUBROUTINE GETTRANSPOSEI (IROW,ICOL,VAL,N,M,NZ,IROWT,ICOLT,VALT)

!     GET TRANSPORSE OF A MATRIX (ELEMENTS ARE INTEGERS)
!-----------------------------------------------------------------------
!      IROW,ICOL,VAL : ORIGINAL MATRIX
!      N : THE NUMBER OF ROWS
!      M : THE NUMBER OF COLUMNS
!      NZ : THE NUMBER OF NONZERO ELEMENTS
!      IROWT,ICOLT,VALT : TRANSPORSE MATRIX
!
      USE DEFDP
      IMPLICIT NONE
      
!     DUMMY ARGUMENTS
      INTEGER, INTENT(IN)      :: N,M,NZ
      INTEGER, INTENT(IN),DIMENSION(:) :: IROW(N+1),ICOL(NZ), VAL(NZ)
      INTEGER, INTENT(OUT),DIMENSION(:) :: IROWT(M+1),ICOLT(NZ),VALT(NZ)

!      LOCAL ARGUMENTS      
      INTEGER :: I,J,K,JP,IAA,IAB
!      
!      --- END OF DECLARATION ---
!
      IROWT=0
      IROWT(1) = 1
      IROWT(2) = 1

      
      DO I = 1, IROW(N + 1) - 1
            J = ICOL(I) + 2
            IF ( J .LE. M + 1) IROWT(J) = IROWT(J) + 1
      END DO


      DO I = 3, M + 1
            IROWT(I) = IROWT(I) + IROWT(I-1)
      END DO

      DO I = 1,N
            IAA = IROW(I)
            IAB = IROW(I+1) - 1
            IF (IAB .GE. IAA) THEN 
                  DO JP = IAA, IAB
                        J = ICOL(JP) + 1
                        K = IROWT(J)
                        ICOLT(K) = I
                        VALT(K) = VAL(JP)
                        IROWT(J) = K + 1
                  END DO
            END IF         
      END DO
      
      RETURN
      END SUBROUTINE GETTRANSPOSEI

