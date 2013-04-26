!---------------- SUBROUTINE -----------------
       SUBROUTINE ALLOCCA

!    WRITETEN BY YOUSU CHEN, PNNL, 6 JUNE 2007
!    LAST MODIFIED BY YOUSU CHEN, PNNL, 6 JUNE 2007
!    $ID: READCA.F V1.01 2007/6/6
!
!----------- DESCRIPTION ---------------------
!      THIS PROGRAM READ THE CONTINGENCY LIST FOR CAFILE
!
      USE CAMODULE
      USE BRCHMODULE
      USE GENMODULE,ONLY: NGON

      USE FLAGS
      IMPLICIT NONE

      INTEGER :: I,J,ID
!      
!      --- END OF DECLARATION ---
!

      IF (CAOPT.EQ.1) THEN
         OPEN(8,FILE=CAFILE)
         READ(8,*) NCA0
	 CLOSE(8)
         IF (M.EQ.2) THEN
            NCA=NCA0*(NCA0-1)/2
         ELSE
            NCA=NCA0
         ENDIF
         ! FIND THE NUMBER OF LINES IN CALIST FILE WHEN M=0
         IF (M.EQ.0) THEN
           NLINES = 0
           OPEN(9,FILE=CAFILE)
           DO 
             READ(9,*, END=10)
             NLINES = NLINES + 1
           ENDDO
10         CLOSE(9)
         ENDIF
         PRINT *, 'NLINES =',NLINES
      ELSE  !CONSIDER BRANCH CA ONLY 
         NCA0=NBRCH
         IF (M.EQ.1) THEN
            NCA = NBRCH
         ELSEIF (M.EQ.2) THEN
            NCA = NBRCH * (NBRCH -1)/2
         ELSE
            PRINT *, 'ERROR, WHEN CAOPT=0, M SHOULD BE 1 OR 2'
            STOP
         ENDIF
      ENDIF

      IF (ALLOCATED(SEL_BRID)) DEALLOCATE(SEL_BRID)
      IF (ALLOCATED(CA_BRID)) DEALLOCATE(CA_BRID)
      IF (ALLOCATED(CA_FROM)) DEALLOCATE(CA_FROM)
      IF (ALLOCATED(CA_TO)) DEALLOCATE(CA_TO)
      IF (ALLOCATED(CA_IROW)) DEALLOCATE(CA_IROW)
      IF (ALLOCATED(N2_IDX1)) DEALLOCATE(N2_IDX1)
      IF (ALLOCATED(N2_IDX2)) DEALLOCATE(N2_IDX2)

      IF (NCASE.EQ.0) NCAN=NCA
      IF (NCASE.LT.NCA) NCAN=NCASE
      IF (NCASE.GE.NCA) NCAN=NCA
!     DYNAMIC ALLOCATE MEMORY ARRAY      
      IF (M.EQ.1) THEN
        ALLOCATE(SEL_BRID(NCA0), STAT=ERROR)  
        IF (ERROR /=0) THEN
          PRINT *, " ----------------- ERROR ------------------------" 
          PRINT *, " PROGRAM COULD NOT ALLOCATE SPACE FOR CA ARRAYS" 
          STOP
        ENDIF
      ENDIF
      IF (M.EQ.2) THEN
        ALLOCATE(SEL_BRID(NCA0),N2_IDX1(NCAN),N2_IDX2(NCAN),
     & STAT=ERROR)  
        IF (ERROR /=0) THEN
          PRINT *, " ----------------- ERROR ------------------------" 
          PRINT *, " PROGRAM COULD NOT ALLOCATE SPACE FOR CA ARRAYS" 
          STOP
        ENDIF
      ENDIF
      IF (M.EQ.0) THEN
        ALLOCATE(CA_BRID(NLINES-NCA0-1),CA_NAME(NCA0),N_OF_CA(NCA0), 
     & CA_INDEX(NCA0),CA_FROM(NLINES-NCA0-1), CA_TO(NLINES-NCA0-1), 
     & CA_IROW(NCA0),STAT=ERROR)  
        IF (ERROR /=0) THEN
          PRINT *, " ------------------ ERROR ------------------------" 
          PRINT *, " COULD NOT ALLOCATE SPACE FOR CA (M=0) ARRAYS" 
          STOP
        END IF
        CA_BRID=0
        CA_NAME=''
        CA_INDEX=0
        CA_FROM=0
        CA_TO=0
        N_OF_CA=0
      ENDIF

!     ASSIGN N-2 INDEX
      ID=0
      IF (M.EQ.2) THEN
         DO I=1,NCA0-1
            DO J=I+1,NCA0
               ID = ID +1
               IF (ID.GT.NCAN) GOTO 100
               N2_IDX1(ID)=I
               N2_IDX2(ID)=J
            ENDDO
         ENDDO
      ENDIF
  100 CONTINUE
      RETURN
      END

