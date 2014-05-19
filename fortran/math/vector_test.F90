  ! ----------------------------------------------------------------
  ! file: vector_test.F90
  ! ----------------------------------------------------------------
  ! ----------------------------------------------------------------
  ! Battelle Memorial Institute
  ! Pacific Northwest Laboratory
  ! ----------------------------------------------------------------
  ! ----------------------------------------------------------------
  ! Created May 15, 2014 by William A. Perkins
  ! Last Change: 2014-05-16 09:24:37 d3g096
  ! ----------------------------------------------------------------
PROGRAM vector_test

  USE gridpack_parallel
  USE gridpack_math

  IMPLICIT NONE

  TYPE (communicator) :: comm
  INTEGER(c_int) :: ma_stack = 200000, ma_heap = 200000

  TYPE (vector) :: vec
  INTEGER :: me, lo, hi, i
  COMPLEX(c_double_complex) :: myx
  COMPLEX(c_double_complex) :: norm1, norm2, norminf

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
  DO i = lo, hi
     myx = cmplx(i, me)
     CALL vec%set_element(i, myx)
  END DO
  CALL vec%ready()
  CALL vec%print()

  DO i = lo, hi
     myx = cmplx(i, me)
     CALL vec%add_element(i, myx)
  END DO
  CALL vec%ready()
  CALL vec%print()

  norm1 = vec%norm1()
  norm2 = vec%norm2()
  norminf = vec%norm_infinity()
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
  CALL comm%finalize()
  CALL gridpack_finalize_math()
  CALL gridpack_finalize_parallel()
  
END PROGRAM vector_test
