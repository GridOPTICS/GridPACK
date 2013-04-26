      SUBROUTINE READSIZE(FILENAME)
      
!    WRITETEN BY YOUSU CHEN, PNNL, 14 MARCH, 2007
!    LAST MODIFIED BY YOUSU CHEN, PNNL, 09 APRIL, 2007
!    $ID: READSIZE.F, V1.03 2007/4/09
!
!    REVISION LOG
!    1. 2007/4/1  ADDED PTI FORMAT
!    2. 2007/4/9  USED '-999' AS A FLAG FOR IEEE FORMAT, SINCE FOR SOME IEEE
!       INPUT FILES, THE LOCATION OF THE NUMBER OF BUS/BRANCH IS DIFFERENT.
!
! ------ DESCRIPTION ---------
!      THIS SUBROUTINE READS INPUT FILE (IEEE OR PTI) AND FIND THE 
!      THE NUMBER OF BUSES, GENERATORS, AND BRANCHES
!


      USE DEFDP
      USE BUSMODULE, ONLY: NB,NPV,NPQ
      USE GENMODULE, ONLY: NG
      USE BRCHMODULE, ONLY: NBRCH

      IMPLICIT NONE
      
!     DUMMY ARGUMENT DECLARATION
      CHARACTER*100, INTENT(IN) :: FILENAME

!     LOCAL ARGUMENT DECLARATION
      INTEGER :: I,TEMP,TMP,IFLG,IBUS,IGEN,IBRCH
      CHARACTER*100::TXT
!      
!      --- END OF DECLARATION ---
!      
!
! --------INITIALIZATION--------
!
      NB = 0
      NG = 0
      NPV = 0
      NPQ = 0
      NBRCH = 0
      IBUS = 0
      IGEN = 0
      IBRCH =0
      NXFMR_ADJ = 0
      NDC = 0
      NAREA = 0
      NSHUNT = 0
      IFMT = 2


      OPEN(1,FILE=FILENAME)
!
! ----- IEEE COMMON DATA FORMAT ---
!
      IF (IFMT .EQ. 1) THEN
        IFLG = 0
        DO WHILE (IFLG .EQ.0)
           READ(1,'(A)') TXT
           IBUS = IBUS + 1
           IF (TXT(1:4) .EQ. '-999') IFLG = 1
        END DO
        NB = IBUS - 3 
        DO WHILE (IFLG .EQ. 1)
          READ(1,'(A)') TXT
          IBRCH = IBRCH + 1
          IF (TXT(1:4) .EQ. '-999') IFLG = 2
        END DO
        NBRCH = IBRCH-2
        CLOSE(1)
        OPEN(1,FILE=FILENAME)
          READ(1,*)
          READ(1,*)
        DO I=1,NB
          READ(1,'(24X,I2)'), TEMP
          IF (TEMP .GE. 2) IGEN=IGEN+1
        END DO
        NG = IGEN
      END IF

!
! ------PTI FORMAT ----
!

      IF (IFMT .EQ. 2) THEN
        READ(1,*)
        READ(1,*)
        READ(1,*)
        IFLG = 0
        DO WHILE (IFLG .EQ. 0)  
          READ(1,'(I7,2X,I1,A)') TEMP,TMP,TXT
          IF (TEMP .NE. 0) THEN 
            NB=NB + 1
          ELSE 
            IFLG = 1
          END IF
          IF (TEMP .EQ. 0) IFLG = 1
        END DO   
        DO WHILE (IFLG .EQ. 1) 
          READ(1,'(A1)') TXT
          IF (TXT .NE. '0') THEN 
            NG = NG + 1
          ELSE 
            IFLG = 2
          END IF    
        END DO
        DO WHILE (IFLG .EQ. 2) 
          READ(1,'(A1)') TXT
          IF (TXT .NE. '0') THEN 
            NBRCH = NBRCH + 1
          ELSE 
            IFLG = 3
          END IF
        END DO
        DO WHILE (IFLG .EQ. 3) 
          READ(1,'(A1)') TXT
          IF (TXT .NE. '0') THEN 
            NXFMR_ADJ = NXFMR_ADJ + 1
          ELSE 
            IFLG = 4
          END IF
        END DO
        DO WHILE (IFLG .EQ. 4) 
          READ(1,'(A1)') TXT
          IF (TXT .NE. '0') THEN 
            NAREA = NAREA + 1
          ELSE 
            IFLG = 5
          END IF
        END DO
        DO WHILE (IFLG .EQ. 5) 
          READ(1,'(A1)') TXT
          IF (TXT .NE. '0') THEN 
            NDC = NDC + 1 
          ELSE 
            IFLG = 6
          END IF
        END DO
        DO WHILE (IFLG .EQ. 6) 
          READ(1,'(A1)') TXT
          IF (TXT .NE. '0') THEN 
            NSHUNT = NSHUNT + 1 
          ELSE 
            IFLG = 7
          END IF
        END DO
      END IF

!
! -----------PRINT SYSTEM INFORMATION--------
!
      CLOSE(1)
!      WRITE(*,*) 'NUMBER OF BUSES:',NB
!      WRITE(*,*) 'NUMBER OF BRANCHES:',NBRCH
!      WRITE(*,*) 'NUMBER OF GENERATORS:',NG

      END SUBROUTINE READSIZE


