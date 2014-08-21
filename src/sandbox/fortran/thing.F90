! ----------------------------------------------------------------
! MODULE thing_message
! ----------------------------------------------------------------
MODULE thing_module
  USE iso_c_binding
  IMPLICIT NONE
  
  PRIVATE

  TYPE, PUBLIC :: thing
   CONTAINS
     PROCEDURE :: message => base_message
  END type thing

  TYPE, PUBLIC, EXTENDS(thing) :: thing1
   CONTAINS
     PROCEDURE :: message => message1
  END type thing1

  TYPE, PUBLIC, EXTENDS(thing) :: thing2
   CONTAINS
     PROCEDURE :: message => message2
  END type thing2

  TYPE, PUBLIC :: thingwrapper
     CLASS(thing), POINTER :: t
  END type thingwrapper

CONTAINS

  SUBROUTINE base_message(this)
    IMPLICIT NONE
    CLASS(thing), INTENT(IN) :: this
    WRITE (*, *) "base thing"
  END SUBROUTINE base_message

  SUBROUTINE message1(this)
    IMPLICIT NONE
    CLASS (thing1), INTENT(IN) :: this
    WRITE (*, *) "thing 1"
  END SUBROUTINE message1
  
  SUBROUTINE message2(this)
    IMPLICIT NONE
    CLASS (thing2), INTENT(IN) :: this
    WRITE (*, *) "thing 2"
  END SUBROUTINE message2

  TYPE(c_ptr) FUNCTION create_thing(itype) BIND(c)
    IMPLICIT NONE
    INTEGER(c_int), VALUE, INTENT(IN) :: itype
    TYPE(thingwrapper), POINTER :: w
    CLASS(thing), POINTER :: t
    
    SELECT CASE (itype)
    CASE (0)
       ALLOCATE(thing::t)
    CASE (1)
       ALLOCATE(thing1::t)
    CASE (2)
       ALLOCATE(thing2::t)
    CASE DEFAULT

    END SELECT
    IF (ASSOCIATED(t)) THEN
       ALLOCATE(w)
       w%t => t
    ELSE 
       NULLIFY(w)
    END IF
    create_thing = C_LOC(w)
  END FUNCTION create_thing

  SUBROUTINE thing_message(ptr) BIND(c)
    IMPLICIT NONE
    TYPE(c_ptr), INTENT(IN), VALUE :: ptr
    TYPE(thingwrapper), POINTER :: w
    CALL C_F_POINTER(ptr, w)
    IF (ASSOCIATED(w)) THEN
       CALL w%t%message()
    END IF
  END SUBROUTINE thing_message


  SUBROUTINE destroy_thing(ptr) BIND(c)
    IMPLICIT NONE
    TYPE(c_ptr), INTENT(IN), VALUE :: ptr
    TYPE(thingwrapper), POINTER :: w
    CALL C_F_POINTER(ptr, w)
    IF (ASSOCIATED(w)) THEN
       IF (ASSOCIATED(w%t)) DEALLOCATE(w%t)
       DEALLOCATE(w)
    END IF
  END SUBROUTINE destroy_thing
    
  
END MODULE thing_module

