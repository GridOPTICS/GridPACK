! ------------------- SUBROUTINE -------------------
!      SUBROUTINE NEWTONPF(TOL,MAX_IT,ITER,CONVERGED,
!     &ISOLATE,SOLOPT)
      SUBROUTINE NEWTONPF
!
!     WRITTEN BY YOUSU CHEN, PNNL, 3/15/2007
!     LAST MODIFIED BY YOUSU CHEN, PNNL, 4/20/2007
!     $ID NEWTONPF.F V1.06: 2007/4/20
!
!     REVISION LOG
!     2007/3/15
!     USE UDU.F AS SPARSE SOLVER (DISCARD)
!     2007/3/16
!     USE GAUSS-SIEDAL.F AS SPARSE SOLVER (DISCARD)
!     2007/3/20
!     USE SPLIB LIBRARY AS SPARSE SOLVER (NO GOOD FOR 14090 BUSES) 
!     2007/4/20
!     (1) USE SUPERLU LIBRARY AS SPARSE SOLVER (GOOD FOR ALL)
!     (2) MOVE SUBROUTINES TO ALLO.F MODULE SUBS
!
! --------------- DESCRIPTION -----------------------
!
!     SOLVES THE POWER FLOW UDSING NEWTON'S METHOD
!
! ---------------------------------------------------------
      USE INTRIFUNS
      USE BUSMODULE,ONLY:NB,NPV,NPQ,VA,VM
      USE YBUSMODULE
      USE OTHERMODULE,ONLY: PV,PQ,V,V0,SBUS,V0SAVE
      USE JACOBI
      USE SUPERLU
      USE SUBS
      USE DEFDP
     

!      USE ALLMODULE 
      IMPLICIT NONE

      ! DUMMY ARGUMENTS
!      INTEGER, INTENT(IN) :: MAX_IT,ISOLATE,SOLOPT
!      REAL(KIND=DP), INTENT(IN) :: TOL
!      INTEGER, INTENT(INOUT) :: ITER, CONVERGED

      ! LOCAL ARGUMENTS
      
      INTEGER :: I,NPVQ,J,NZJ11,NZJ12,NZJ21,NZJ22,NZJ,NBJ,K
      INTEGER :: COUNTER
      INTEGER, DIMENSION(:) :: PVQ(NPV+NPQ)
      INTEGER, DIMENSION(:) :: IROWJ11(NB+1), ICOLJ11(NZYBUS)
      INTEGER, DIMENSION(:) :: IROWJ12(NB+1), ICOLJ12(NZYBUS)
      INTEGER, DIMENSION(:) :: IROWJ21(NB+1), ICOLJ21(NZYBUS)
      INTEGER, DIMENSION(:) :: IROWJ22(NB+1), ICOLJ22(NZYBUS)
!      INTEGER, DIMENSION(:) :: IROW12(NPV+NPQ+1), ICOL12(2*NZYBUS)
!      INTEGER, DIMENSION(:) :: IROW34(NPQ+1), ICOL34(2*NZYBUS)
!      INTEGER, DIMENSION(:) :: IROWJ(NPV+2*NPQ+1), ICOLJ(4*NZYBUS)
      REAL(KIND=DP) :: NORMF, TEMPNORMF
      REAL(KIND=DP), DIMENSION(:) :: F(NPV+2*NPQ),DX(NPV+2*NPQ)
      REAL(KIND=DP), DIMENSION(NZYBUS) :: REALDVA,IMGDVA,REALDVM,IMGDVM
      REAL(KIND=DP), DIMENSION(NZYBUS) :: VALJ11
      REAL(KIND=DP), DIMENSION(NZYBUS) :: VALJ12
      REAL(KIND=DP), DIMENSION(NZYBUS) :: VALJ21
      REAL(KIND=DP), DIMENSION(NZYBUS) :: VALJ22
!      REAL(KIND=DP), DIMENSION(2*NZYBUS) :: VAL12
!      REAL(KIND=DP), DIMENSION(2*NZYBUS) :: VAL34
      REAL(KIND=DP), DIMENSION(4*NZYBUS) :: VALJ
      COMPLEX(KIND=DP), DIMENSION(NB) :: MIS,IBUS
      COMPLEX(KIND=DP), DIMENSION(NZYBUS) :: DSBUS_DVM,DSBUS_DVA
!
!----- SPLIB ARGUMENTS-----
!
      INTEGER :: IRELTV,PCMETH,SOLMETH,SCALING,IPARMS(2),IREDO,ERRFLG
      INTEGER :: IUNIT(16),WORK(1000)
      INTEGER :: INSTLVL,MAXIT
      REAL(KIND=DP) :: DL(NB),DR(NB),RPARMS(2),TOLE
!
!----- SUPERLU ARGUMENT----
!

      INTEGER:: IOPT,INFO,LDB,NRHS
      INTEGER*8 :: FACTORS
      
!      
!      --- END OF DECLARATION ---
!

!      CALL ALLOCJACOBI
!      CALL ALLOCSUPERLU
      
!      PRINT *, ''
!      PRINT *, ''

!
!--- SPLIB PARAMETERS
!
		  
      IRELTV = 0
      PCMETH = 1
      SOLMETH = 1 
      SCALING = 0
      IPARMS(1) = 1
      IPARMS(2) = 0
      RPARMS(1) = 0.D+0
      RPARMS(2) = 0.D+0
      IREDO = 0
      INSTLVL = -1
      MAXIT=1000
      TOLE = 1.0D-8


      DO I = 1,NPV
            PVQ(I) = PV(I)
      END DO
      DO I = 1,NPQ
            PVQ(NPV+I) = PQ(I)
      END DO
      NPVQ=NPV+NPQ

!      
!     INITIALIZE      
!
      V0=V0SAVE
      V = V0
      VA = DATAN2(AIMAG(V0),REAL(V0))
      VM = ABS(V0)

      CONVERGED = 0
      ITER = 0      
      COUNTER=0
      
      VALJ11=0
      VALJ12=0
      VALJ21=0
      VALJ22=0
      VALJ=0
!      
!       MIS = V .* CONJ(YBUS * V) - SBUS
!
      CALL MATMULVEC(IROWYBUS,ICOLYBUS,VALYBUS,V,IBUS,NB,NB,NZYBUS)
      
                  
      DO I = 1,NB
         MIS(I) = V(I)*CONJG(IBUS(I))-SBUS(I)
      END DO
      
      F = 0
      DO I = 1,NPV
         F(I) = REAL(MIS(PV(I)))
      END DO
      DO I = 1,NPQ
         F(NPV+I) = REAL(MIS(PQ(I)))
      
      END DO
      DO I = 1,NPQ
         F(NPV+NPQ+I) = AIMAG(MIS(PQ(I)))      
      END DO
!
!       CHECK TOLERANCE
!
      NORMF = MAXVAL(ABS(F))

      IF (NORMF .LT. TOL) CONVERGED = 1

      DO WHILE (CONVERGED .NE. 1 .AND. ITER .LE. MAX_IT)
            ITER = ITER + 1

!	    IF (ITER .EQ. 1) THEN
          
            CALL DSBUS_DV(IROWYBUS,ICOLYBUS,VALYBUS,V,IBUS
     &             ,DSBUS_DVM,DSBUS_DVA,NB,NZYBUS,ISOLATE)


            REALDVA = REAL(DSBUS_DVA)
            IMGDVA = AIMAG(DSBUS_DVA)
            REALDVM = REAL(DSBUS_DVM)
            IMGDVM = AIMAG(DSBUS_DVM)
      !
      !   CREATE J11 = REALDVA(PVQ,PVQ)
      !
            CALL JSUBSET(IROWYBUS,ICOLYBUS,REALDVA,PVQ,PVQ,NB
     &      ,NZYBUS,NPVQ,NPVQ,IROWJ11,ICOLJ11,VALJ11,NZJ11)
      
      !
      !   CREATE J12 = REALDVM(PVQ,PQ)
      !
            CALL JSUBSET(IROWYBUS,ICOLYBUS,REALDVM,PVQ,PQ,NB
     &      ,NZYBUS,NPVQ,NPQ,IROWJ12,ICOLJ12,VALJ12,NZJ12)

      !
      !   CREATE J21 = IMGDVA(PQ,PVQ)
      !
            CALL JSUBSET(IROWYBUS,ICOLYBUS,IMGDVA,PQ,PVQ,NB
     &      ,NZYBUS,NPQ,NPVQ,IROWJ21,ICOLJ21,VALJ21,NZJ21)

      !
      !   CREATE J22 = REALDVM(PQ,PQ)
      !
            CALL JSUBSET(IROWYBUS,ICOLYBUS,IMGDVM,PQ,PQ,NB
     &      ,NZYBUS,NPQ,NPQ,IROWJ22,ICOLJ22,VALJ22,NZJ22)      
 
      !      
      !      CREATE [ J11 J12]
      !

            CALL SIDECATENATE(IROWJ11,ICOLJ11,VALJ11,IROWJ12,ICOLJ12, 
     &            VALJ12,IROW12,ICOL12,VAL12,NPVQ,NPVQ,NZJ11,NZJ12)
      
      !      
      !      CREATE [ J21 J22]
      !

            CALL SIDECATENATE(IROWJ21,ICOLJ21,VALJ21,IROWJ22,ICOLJ22, 
     &               VALJ22,IROW34,ICOL34,VAL34,NPQ,NPVQ,NZJ21,NZJ22)

      !                 
      !      CREATE J = [ J11 J12
      !                   J21 J22]   

            CALL UPDOWNCATENATE(IROW12,ICOL12,VAL12,IROW34,ICOL34,VAL34,
     &      IROWJ,ICOLJ,VALJ,NPVQ,NPQ,NPVQ+NPQ,NZJ11+NZJ12,NZJ21+NZJ22)
!            END IF
            NBJ = NPVQ+NPQ
            NZJ = NZJ11+NZJ12+NZJ21+NZJ22


                        
!        CALL SOLVE(IROWJ, ICOLJ, -VALJ, DX, F, NBJ,NZJ)
!        CALL GAUSSSEIDEL(IROWJ, ICOLJ, -VALJ, DX, F, NBJ,NZJ)    
!
!        SLOVE AX = B BY (1) SPLIB (2) SUPERLU


            IF (SOLOPT .EQ. 1) THEN
!           SPLIB SOLVER	    

!                CALL SPLIB(-VALJ,ICOLJ,IROWJ,NBJ,DX,F,
!     &          WORK,1000,TOLE,IRELTV,MAXIT,
!     &          PCMETH,SOLMETH,SCALING,DL,DR,IPARMS,RPARMS,IREDO,ERRFLG,
!     &          IUNIT, INSTLVL)
!                PRINT *, 'ERRFLG =', ERRFLG

	    ELSE ! SUPERLU

!            CALL MYPRINTSPARSE(IROWJ,ICOLJ,VALJ,F,NBJ)


!         CONVERT TO SUPERLU STORAGE FORMAT (COMPRESS COLUMN)

	        CALL GETTRANSPOSE (IROWJ,ICOLJ,VALJ,NBJ,NBJ,NZJ,COLPTR,
     &	                       ROWIND,VALUES)

!
!--- SUPERLU PARAMETERS
!
                 NRHS = 1
                 LDB = NBJ
	         DX = F
                 VALUES = -VALUES

             
                 IOPT = 1
	         CALL C_FORTRAN_DGSSV(IOPT,NBJ,NZJ,NRHS,VALUES,ROWIND,
     &                  COLPTR,DX, LDB, FACTORS, INFO )
!                 IF (INFO .EQ. 0) THEN
!                     WRITE (*,*) 'FACTORIZATION SUCCEEDED'
!                 ELSE
!                     WRITE(*,*) 'INFO FROM FACTORIZATION = ', INFO
!                 END IF
					 
!                 IOPT = 2
!	         CALL C_FORTRAN_DGSSV(IOPT,NBJ,NZJ,NRHS,VALUES,ROWIND,
!     &                  COLPTR,DX, LDB, FACTORS, INFO )
!                 IF (INFO .EQ. 0) THEN
!                     WRITE (*,*) 'SOLVE SUCCEEDED'
!	         ELSE
!	             WRITE(*,*) 'INFO FROM TRIANGULAR SOLVE = ', INFO
!	         END IF
		 
!                 IOPT = 3
!	         CALL C_FORTRAN_DGSSV(IOPT,NBJ,NZJ,NRHS,VALUES,ROWIND,
!     &                 COLPTR, DX, LDB, FACTORS, INFO )

            END IF
	    
!            UPDATE ANGEL AND VM OF V
!
            DO I = 1,NPV      
                  VA(PV(I)) = VA(PV(I)) + DX(I)
            END DO
            DO I = 1,NPQ
                  VA(PQ(I)) = VA(PQ(I)) + DX(NPV+I)
                  VM(PQ(I)) = VM(PQ(I)) + DX(NPV+NPQ+I)
            END DO

            V = DCMPLX(VM*DCOS(VA),VM*DSIN(VA))
!
!            UPDATE VM AND VA AGAIN IN CASE WE WRAPPED AROUND
!            WITH A NEGATIVE VM
!
            VA = DATAN2(AIMAG(V),REAL(V))
            VM = ABS(V)
!
!            UPDATE MIS AND F
!
            CALL MATMULVEC(IROWYBUS,ICOLYBUS,VALYBUS,V,IBUS,NB,NB,
     &           NZYBUS)
                  
            DO I = 1,NB
                  MIS(I) = V(I)*CONJG(IBUS(I))-SBUS(I)
            END DO
      
            F = 0
            DO I = 1,NPV
                  F(I) = REAL(MIS(PV(I)))      
            END DO
            
            DO I = 1,NPQ
                  F(NPV+I) = REAL(MIS(PQ(I)))      
            END DO
            
            DO I = 1,NPQ
                  F(NPV+NPQ+I) = AIMAG(MIS(PQ(I)))      
            END DO
            TEMPNORMF = NORMF
            NORMF=0.0
            NORMF = MAXVAL(ABS(F))
            IF ((TEMPNORMF > NORMF) .AND. NORMF >1000) THEN
                COUNTER=COUNTER+1
            ENDIF

            ! CONVERGED
            IF ( ABS(NORMF) .LT. TOL ) THEN
	        CONVERGED = 1
	        EXIT
            ELSE ! DIVERGED
	        IF (( ABS(NORMF) .GE. 10000.0) .OR.
     &              (COUNTER .GE. 4))   THEN
                  CONVERGED = 2
                  EXIT
                ENDIF
	    ENDIF  
      ENDDO
!      PRINT *, ME, ITER, NORMF, CONVERGED

      RETURN
      END 


                   
                        
