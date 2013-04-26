! --------------- SUBROUTINE ------------------
      SUBROUTINE PRINTPF
       
!      WRITTEN BY YOUSU CHEN,PNNL, ON 03/22/2007/       
!      LAST MODIFIED BY YOUSU CHEN, PNNL, ON 03/29/2007
!      $ID: PRINTPF.F V1.02: 2007/03/39
!
! ------------- DESCRIPTION -------------------
!
!      WRITE RESULTS INTO 'PFRESULT.OUT' FILE
!
! ---------------------------------------------
       
      USE ALLMODULE
      IMPLICIT NONE
      INTEGER :: I

      OPEN(2,FILE='pfresult.out',STATUS='UNKNOWN')
      OPEN(3,FILE='basecase.out',STATUS='UNKNOWN')
      OPEN(5,FILE='basevio.out',STATUS='UNKNOWN')
      WRITE(3,'(A3,1X,I1)') 'CTG', 0
      WRITE(3,'(I1)'), 1
      WRITE(3,FMT='(A3)') 'BUS'
      WRITE(2,*) '=================================================='
      WRITE(2,*) ''
      WRITE(2,*) '       SYSTEM SUMMARY                         '
      WRITE(2,*) ''
      WRITE(2,*) '=================================================='
      WRITE(2,*) ' '
      WRITE(2,*) '--------------------------------------------------'
      WRITE(2,*) ' '
      WRITE(2,'(A, I5)') ' THE NUMBER OF BUSES:            ', NB
      WRITE(2,'(A, I5)') ' THE NUMBER OF BRANCHES:      ', NBRCH
      WRITE(2,'(A, I5)') ' THE NUMBER OF GENERATORS:      ', NG
      WRITE(2,'(A, I5)') ' THE NUMBER OF COMMITTED GENERATORS:   ', NGON
      WRITE(2,'(A, I5)') ' THE NUMBER OF SHUNTS:        ', NS
      WRITE(2,*) ' TOTAL GENERATOR CAPACITY '
      WRITE(2,*) '  P(MW)        Q(MVAR) '
      WRITE(2,FMT=105) SUM(PMAX),SUM(QMIN),' TO ',SUM(QMAX) 

      WRITE(2,*) ' '
      WRITE(2,*) '--------------------------------------------------'
      WRITE(2,*) ' '
      WRITE(2,*) '=================================================='
      WRITE(2,*) '       BUS DATA'
      WRITE(2,*) '=================================================='
      WRITE(2,*) ' '
      WRITE(2,FMT=106)  'BUS','VOLTAGE','GENERATION','LOAD'
      WRITE(2,*)'    #   MAG(PU)  ANG(DEG)     P (MW)   Q (MVAR)   
     &  P (MW)   Q (MVAR)   '      
      WRITE(2,*) '-----------------------------------------------------
     &--------------'
  105 FORMAT(F10.3, 2X, F10.3,A,F10.3)
  106 FORMAT(4X,A,6X,A,13X,A,15X,A)
!
      PGP = 0.0
      QGP = 0.0
      
      PGP(GBUS(1)) = PGG(1)
      QGP(GBUS(1)) = QGG(1)
      DO I = 2,NGON
            IF (GBUS(I) .EQ. GBUS(I-1)) THEN
	        PGP(GBUS(I)) = PGG(I) + PGP(GBUS(I-1))
	        QGP(GBUS(I)) = QGG(I) + QGP(GBUS(I-1))
	    ELSE
	        PGP(GBUS(I)) = PGG(I)
	        QGP(GBUS(I)) = QGG(I)
            END IF
      END DO
!      PRINT *, QGG
!      PRINT *, QGP

      DO I = 1,NB
            WRITE(3,FMT=126) BUS_I(I),VM(I)
            IF ((VM(I) .GT. 1.05) .OR. (VM(I) .LT. 0.95)) THEN
                WRITE(5,FMT=126) BUS_I(I),VM(I)
            ENDIF
            IF (PGP(I) .NE. 0 .OR. QGP(I) .NE. 0) THEN      
                  IF (PD(I) .NE.0 .OR. QD(I) .NE. 0) THEN            
                       WRITE(2,FMT=107) BUS_I(I),VM(I),VA(I),PGP(I)
     &                                    ,QGP(I),PD(I),QD(I)
                  ELSE
                       WRITE(2,FMT=108) BUS_I(I),VM(I),VA(I),PGP(I)
     &                            ,QGP(I), '        -         -   '
                  END IF
            ELSE
                  IF (PD(I) .NE.0 .OR. QD(I) .NE. 0) THEN            
                       WRITE(2,FMT=109) BUS_I(I),VM(I),VA(I), 
     &                   '        -         -   ',PD(I),QD(I)
                  ELSE
                       WRITE(2,FMT=110) BUS_I(I),VM(I),VA(I), 
     &              '        -         -           -           - '
                  END IF
            END IF
      END DO
  126 FORMAT(I6,2X,F6.4) 
  107 FORMAT(I5,1X,F7.3,2X,F9.3,4(2X,F10.2))
  108 FORMAT(I5,1X,F7.3,2X,F9.3,2(2X,F10.2),A)
  109 FORMAT(I5,1X,F7.3,2X,F9.3,A,2(2X,F10.2))
  110 FORMAT(I5,1X,F7.3,2X,F9.3,A)
      WRITE(2,'(25X,A)')'   --------  --------     --------  --------  '
      WRITE(2,FMT=111) 'TOTAL:',SUM(PGG),SUM(QGG),SUM(PD),SUM(QD)       
  111 FORMAT(18X,A,1X,2(F10.3,1X),2X,2(F10.3,3X))

!      CALCULATE LOSS IN BRANCH
!
!      LOSS = BASEMVA * ABS(V(F_BUS(I)) ./ TAP - V(T_BUS(I)) .^ 2 ./ ...
!          (BRANCH(:, BR_R) - J * BRANCH(:, BR_X));
      
!       ADD PHASE SHIFTERS

!      DO I = 1, NBRCH
!            IF (TAP(I) .EQ. 0.0) TAP(I) = 1.0
!      END DO

!     TAP_CPLX = CMPLX(TAP*(COS(-SHIFT)),TAP*SIN(-SHIFT))
      
      DO I = 1, NBRCH
         IF (BR_STATUS(I).EQ.0) THEN
            LOSS(I) = 0.0
         ELSE
!            LOSS(I) = BASEMVA*(ABS(V(F_BUS(I))/TAP_CPLX(I) 
!     &      - V(T_BUS(I))))**2/CMPLX(BR_R(I),-BR_X(I))            
            LOSS(I) = BASEMVA*(ABS(V(E2I(F_BUS(I)))/TAP_CPLX(I)
     &      - V(E2I(T_BUS(I)))))**2/DCMPLX(BR_R(I),-BR_X(I))            
         ENDIF
      END DO

      WRITE(2,*) ' '
      WRITE(2,*) '==================================================='
      WRITE(2,*) '|     BRANCH DATA                                '
      WRITE(2,*) '==================================================='
      WRITE(2,*) ' '
      WRITE(2,*) 'BRNCH   FROM   TO    FROM BUS INJECTION   TO BUS 
     &INJECTION     LOSS (I^2 * Z)  '
      WRITE(2,*) '  #     BUS    BUS    P (MW)   Q (MVAR)   P (MW)   Q 
     &(MVAR)   P (MW)   Q (MVAR)'
      WRITE(2,*) ' -----  -----  -----  --------  --------  --------  
     &--------  --------  -------'
      WRITE(3,FMT='(A4)') 'BRCH'
      DO I = 1, NBRCH
      !   MVALIM=SQRT(3.0)*ABS(DCMPLX(PF(I),QF(I)))/VM(F_BUS(I))
!         MVALIM=ABS(DCMPLX(PF(I),QF(I)))
         MVALIM=SQRT(PF(I)*PF(I)+QF(I)*QF(I)) 
         IF (MVALIM .GT. 0.9*MVAHIGH(I)) THEN
             WRITE(5,FMT=114) I,PF(I),QF(I),!MVALIM,MVAHIGH(I),
     & MVALIM/MVAHIGH(I)
         ENDIF
          
         WRITE(2,FMT=112) I, F_BUS(I),T_BUS(I),PF(I),QF(I),PT(I),QT(I),
     &                      REAL(LOSS(I)), AIMAG(LOSS(I))
         WRITE(3,FMT=113) I,PF(I),QF(I)
      END DO
      DO I = 1, NBRCH
         WRITE(3,FMT=115) I,F_BUS(I),T_BUS(I),RATE_A(I)
      ENDDO
      CLOSE(3)
      CLOSE(5)

  112 FORMAT(I5,I7,I7,4(F10.2),F10.3,F10.2)
  113 FORMAT(I6,2X,F12.3,1X,F12.3)
  114 FORMAT(I6,2X,F12.3,1X,F12.3,1X,F12.3)
  115 FORMAT(I6,1X,I6,1X,I6,1X,F12.3)
!  114 FORMAT(I6,2X,F9.2,3(1X,F9.2))
      WRITE(2,'(63X,A)') '--------  --------'
      WRITE(2,'(54X,A,F9.3,1X,F9.2)') 'TOTAL:',
     &      SUM(REAL(LOSS)), SUM(AIMAG(LOSS))
      CLOSE(2)
  
      WRITE(*,*) ' '
      WRITE(*,*) ' POWER FLOW SIMULATION RESULTS ARE SAVED IN FILE:
     & PFRESULT.OUT'
      WRITE(*,*) ' '

      RETURN
      END
                 


