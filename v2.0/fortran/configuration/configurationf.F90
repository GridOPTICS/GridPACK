! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! MODULE gridpack_configuration
! ----------------------------------------------------------------
MODULE gridpack_configuration
  
  USE iso_c_binding
  USE gridpack_parallel

  IMPLICIT NONE
  
  PRIVATE

  TYPE, PUBLIC :: cursor
     TYPE (c_ptr) :: impl = C_NULL_PTR
   CONTAINS
     ! PROCEDURE :: initialize
     PROCEDURE :: initialize
     PROCEDURE :: open
     PROCEDURE :: ok
     PROCEDURE :: finalize
     PROCEDURE :: set_path
     PROCEDURE :: get_cursor
     PROCEDURE :: get_bool
     PROCEDURE :: get_int
     PROCEDURE :: get_double
     PROCEDURE :: get_string
  END type cursor

  INTERFACE
     TYPE (c_ptr) FUNCTION configuration_open(comm, file) BIND(c)
       USE iso_c_binding, ONLY: c_char, c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: comm
       CHARACTER (kind=c_char), INTENT(IN) :: file(*)
     END FUNCTION configuration_open

     TYPE (c_ptr) FUNCTION configuration_cursor(cur, path) BIND(c)
       USE iso_c_binding, ONLY: c_char, c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: cur
       CHARACTER (kind=c_char), INTENT(IN) :: path(*)
     END FUNCTION configuration_cursor

     LOGICAL(c_bool) FUNCTION configuration_get_bool(cur, key, flag) BIND(c)
       USE iso_c_binding, ONLY: c_char, c_ptr, c_bool
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: cur
       CHARACTER (kind=c_char), INTENT(IN) :: key(*)
       LOGICAL (c_bool), INTENT(OUT) :: flag
     END FUNCTION configuration_get_bool

     LOGICAL(c_bool) FUNCTION configuration_get_int(cur, key, i) BIND(c)
       USE iso_c_binding, ONLY: c_char, c_ptr, c_bool, c_int
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: cur
       CHARACTER (kind=c_char), INTENT(IN) :: key(*)
       INTEGER (c_int), INTENT(OUT) :: i
     END FUNCTION configuration_get_int

     LOGICAL(c_bool) FUNCTION configuration_get_double(cur, key, d) BIND(c)
       USE iso_c_binding, ONLY: c_char, c_ptr, c_bool, c_double
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: cur
       CHARACTER (kind=c_char), INTENT(IN) :: key(*)
       REAL (c_double), INTENT(OUT) :: d
     END FUNCTION configuration_get_double

     LOGICAL(c_bool) FUNCTION configuration_get_string(cur, key, s, slen) BIND(c)
       USE iso_c_binding, ONLY: c_char, c_ptr, c_bool, c_int
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: cur
       CHARACTER (kind=c_char), INTENT(IN) :: key(*)
       CHARACTER (kind=c_char), INTENT(OUT) :: s(*)
       INTEGER (c_int), INTENT(OUT) :: slen
     END FUNCTION configuration_get_string

     LOGICAL(C_BOOL) FUNCTION configuration_ok(cur) BIND(c)
       USE iso_c_binding, ONLY: c_ptr, c_bool
       IMPLICIT NONE
       TYPE (c_ptr), VALUE, INTENT(IN) :: cur
     END FUNCTION configuration_ok

     SUBROUTINE configuration_destroy(cur) BIND(c)
       USE iso_c_binding, ONLY: c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(INOUT) :: cur
     END SUBROUTINE configuration_destroy

     
  END INTERFACE

CONTAINS

  ! ----------------------------------------------------------------
  ! SUBROUTINE f2cstring
  ! ----------------------------------------------------------------
  SUBROUTINE f2cstring(fstr, cstr)
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: fstr
    CHARACTER(KIND=c_char), INTENT(OUT) :: cstr(*)
    INTEGER :: i, len
    len = len_trim(fstr)
    DO i = 1, len
       cstr(i) = fstr(i:i)
    END DO
    cstr(len+1) = C_NULL_CHAR
  END SUBROUTINE f2cstring

  ! ----------------------------------------------------------------
  ! LOGICAL FUNCTION open
  ! ----------------------------------------------------------------
  LOGICAL FUNCTION open(this, comm, file)
    IMPLICIT NONE
    CLASS (cursor), INTENT(INOUT) :: this
    CLASS (communicator), INTENT(IN) :: comm
    CHARACTER (LEN=*), INTENT(IN) :: file
    CHARACTER (KIND=C_CHAR), ALLOCATABLE :: cfile(:)
    INTEGER :: flen

    IF (this%ok()) THEN
       CALL configuration_destroy(this%impl)
       this%impl = C_NULL_PTR
    END IF
    flen = LEN_TRIM(file)
    ALLOCATE(cfile(flen+1))
    CALL f2cstring(file, cfile)
    this%impl = configuration_open(comm%comm, cfile)
    DEALLOCATE(cfile)
    open = this%ok()
  END FUNCTION open

  ! ----------------------------------------------------------------
  ! SUBROUTINE initialize
  ! ----------------------------------------------------------------
  SUBROUTINE initialize(this)
    IMPLICIT NONE
    CLASS (cursor), INTENT(INOUT) :: this
    this%impl = C_NULL_PTR
  END SUBROUTINE initialize

  ! ----------------------------------------------------------------
  ! LOGICAL FUNCTION ok
  ! ----------------------------------------------------------------
  LOGICAL FUNCTION ok(this) 
    IMPLICIT NONE
    CLASS (cursor), INTENT(INOUT) :: this
    ok = configuration_ok(this%impl)
  END FUNCTION ok

  ! ----------------------------------------------------------------
  ! TYPE(CURSOR) FUNCTION get_cursor
  ! ----------------------------------------------------------------
  TYPE(CURSOR) FUNCTION get_cursor(this, path)
    IMPLICIT NONE
    CLASS (cursor), INTENT(IN) :: this
    CHARACTER (LEN=*), INTENT(IN) :: path
    CHARACTER(c_char), ALLOCATABLE :: cpath(:)
    INTEGER :: flen
    flen = LEN_TRIM(path)
    ALLOCATE(cpath(flen+1))
    CALL f2cstring(path, cpath)
    get_cursor%impl = configuration_cursor(this%impl, cpath)
    DEALLOCATE(cpath)
  END FUNCTION get_cursor

  ! ----------------------------------------------------------------
  ! SUBROUTINE set_path
  ! ----------------------------------------------------------------
  SUBROUTINE set_path(this, path)
    IMPLICIT NONE
    CLASS (cursor), INTENT(INOUT) :: this
    CHARACTER (LEN=*), INTENT(IN) :: path
    CHARACTER(c_char), ALLOCATABLE :: cpath(:)
    INTEGER :: flen
    TYPE(c_ptr) :: cnew
    flen = LEN_TRIM(path)
    ALLOCATE(cpath(flen+1))
    CALL f2cstring(path, cpath)
    cnew = configuration_cursor(this%impl, cpath)
    CALL configuration_destroy(this%impl);
    this%impl = cnew
  END SUBROUTINE set_path


  ! ----------------------------------------------------------------
  ! LOGICAL FUNCTION get_bool
  ! ----------------------------------------------------------------
  LOGICAL FUNCTION get_bool(this, key, flag)
    IMPLICIT NONE
    CLASS (cursor), INTENT(IN) :: this
    CHARACTER (LEN=*), INTENT(IN) :: key
    LOGICAL, INTENT(OUT) :: flag
    CHARACTER (c_char), ALLOCATABLE :: ckey(:)
    LOGICAL(c_bool) :: cflag
    INTEGER :: flen
    flen = LEN_TRIM(key)
    ALLOCATE(ckey(flen+1))
    CALL f2cstring(key, ckey)
    get_bool = configuration_get_bool(this%impl, ckey, cflag)
    DEALLOCATE(ckey)
    IF (get_bool) flag = cflag
  END FUNCTION get_bool

  ! ----------------------------------------------------------------
  ! LOGICAL FUNCTION get_int
  ! ----------------------------------------------------------------
  LOGICAL FUNCTION get_int(this, key, i)
    IMPLICIT NONE
    CLASS (cursor), INTENT(IN) :: this
    CHARACTER (LEN=*), INTENT(IN) :: key
    INTEGER, INTENT(OUT) :: i
    CHARACTER (c_char), ALLOCATABLE :: ckey(:)
    INTEGER(c_int) :: ci
    INTEGER :: flen
    flen = LEN_TRIM(key)
    ALLOCATE(ckey(flen+1))
    CALL f2cstring(key, ckey)
    get_int = configuration_get_int(this%impl, ckey, ci)
    DEALLOCATE(ckey)
    IF (get_int) i = ci
  END FUNCTION get_int

  ! ----------------------------------------------------------------
  ! LOGICAL FUNCTION get_double
  ! ----------------------------------------------------------------
  LOGICAL FUNCTION get_double(this, key, d)
    IMPLICIT NONE
    CLASS (cursor), INTENT(IN) :: this
    CHARACTER (LEN=*), INTENT(IN) :: key
    DOUBLE PRECISION, INTENT(OUT) :: d
    CHARACTER (c_char), ALLOCATABLE :: ckey(:)
    REAL(c_double) :: cd
    INTEGER :: flen
    flen = LEN_TRIM(key)
    ALLOCATE(ckey(flen+1))
    CALL f2cstring(key, ckey)
    get_double = configuration_get_double(this%impl, ckey, cd)
    DEALLOCATE(ckey)
    IF (get_double) d = cd
  END FUNCTION get_double

  ! ----------------------------------------------------------------
  ! LOGICAL FUNCTION get_string
  ! ----------------------------------------------------------------
  LOGICAL FUNCTION get_string(this, key, s)
    IMPLICIT NONE
    CLASS (cursor), INTENT(IN) :: this
    CHARACTER (LEN=*), INTENT(IN) :: key
    CHARACTER(len=*), INTENT(OUT) :: s
    CHARACTER (c_char), ALLOCATABLE :: ckey(:)
    CHARACTER(512) :: stmp
    CHARACTER(c_char) :: cs(512)
    INTEGER(C_INT) slen
    INTEGER :: flen, i
    flen = LEN_TRIM(key)
    ALLOCATE(ckey(flen+1))
    CALL f2cstring(key, ckey)
    get_string = configuration_get_string(this%impl, ckey, cs, slen)
    DEALLOCATE(ckey)
  !
  !  Use the extra copy to put a Fortran string termination character
  !  at the end of the string
  !
    do i = 1, slen
      stmp(i:i) = cs(i)
    end do
    s = stmp(1:slen)
  END FUNCTION get_string

  ! ----------------------------------------------------------------
  ! SUBROUTINE finalize
  ! ----------------------------------------------------------------
  SUBROUTINE finalize(this)
    USE iso_c_binding
    IMPLICIT NONE
    CLASS (cursor), INTENT(INOUT) :: this
    IF (this%ok()) THEN
       CALL configuration_destroy(this%impl)
    END IF
  END SUBROUTINE finalize


END MODULE gridpack_configuration
