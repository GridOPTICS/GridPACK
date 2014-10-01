! ----------------------------------------------------------------
! file: vector_test.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created May 15, 2014 by William A. Perkins
! Last Change: 2014-10-01 07:33:05 d3g096
! ----------------------------------------------------------------
PROGRAM vector_test

  USE gridpack_parallel
  USE gridpack_math

  IMPLICIT NONE

  TYPE (communicator) :: comm
  INTEGER(c_int) :: ma_stack = 200000, ma_heap = 200000

  TYPE (vector) :: vec, vec2
  INTEGER :: me, lo, hi, i, p
  COMPLEX(c_double_complex) :: myx
  COMPLEX(c_double_complex) :: norm1, norm2, norminf
  COMPLEX(c_double_complex), ALLOCATABLE :: xall(:)

  CALL gridpack_initialize_parallel(ma_stack, ma_heap)
  CALL gridpack_initialize_math()
  CALL comm%initialize()

  me = comm%rank()

  CALL vec%initialize(comm, 5)

  WRITE(*, *) me, ": I own ", vec%local_size(), " of ", vec%size(), " elements"

  myx = cmplx( me, me )
  CALL vec%fill(myx)
  CALL vec%print()

  CALL comm%barrier()

  CALL vec%local_index_range(lo, hi)
  !WRITE (*,*) lo, hi
  DO i = lo, hi-1
     myx = cmplx(i, me)
     CALL vec%set_element(i, myx)
  END DO
  CALL vec%ready()
  CALL vec%print()

  DO i = lo, hi-1
     myx = cmplx(i, me)
     CALL vec%add_element(i, myx)
  END DO
  CALL vec%ready()
  CALL vec%print()

  IF (me .EQ. 0) THEN
     WRITE (*,*) 'Using get_element():'
  END IF
  DO p = 0, comm%size()
     IF (p .EQ. me) THEN
        DO i = lo, hi-1
           CALL vec%get_element(i, myx)
           WRITE (*,"(I2, ': ', I2, ': (', F5.1, ', ', F5.1, ')')"), me, i, myx
        END DO
     END IF
     CALL comm%barrier()
  END DO

  IF (me .EQ. 0) THEN
     WRITE (*,*) 'Using get_all_elements():'
  END IF
  ALLOCATE(xall(vec%size()))
  CALL vec%get_all_elements(xall)
  IF (me .EQ. 0) THEN
     DO i = 0, vec%size() - 1
        WRITE (*, "(I2, ': (', F5.1, ', ', F5.1, ')')") i, xall(i+1)
     END DO
  END IF
  DEALLOCATE(xall)

  norm1 = vec%norm1()
  norm2 = vec%norm2()
  norminf = vec%norm_infinity()
  IF (me .EQ. 0) THEN
     WRITE (*,*) norm1
     WRITE (*,*) norm2
     WRITE (*,*) norminf
  END IF

  vec2 = vec%clone()

  norm1 = vec2%norm1()
  norm2 = vec2%norm2()
  norminf = vec2%norm_infinity()
  IF (me .EQ. 0) THEN
     WRITE (*,*) norm1
     WRITE (*,*) norm2
     WRITE (*,*) norminf
  END IF

  CALL vec%zero()
  CALL vec%print()

  norm1 = vec%norm1()
  norm2 = vec%norm2()
  norminf = vec%norm_infinity()
  IF (me .EQ. 0) THEN
     WRITE (*,*) norm1
     WRITE (*,*) norm2
     WRITE (*,*) norminf
  END IF

  CALL vec%finalize()
  CALL vec2%finalize()
  CALL comm%finalize()
  CALL gridpack_finalize_math()
  CALL gridpack_finalize_parallel()
  
END PROGRAM vector_test
