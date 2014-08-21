PROGRAM thing_program
  USE thing_module
  IMPLICIT NONE
  TYPE (thingwrapper) :: w
  CLASS (thing), POINTER :: ptr
  ALLOCATE(thing1::ptr)
  w%t => ptr
  CALL w%t%message()
  DEALLOCATE(w%t)
  ALLOCATE(thing2::ptr)
  w%t => ptr
  CALL w%t%message()
  DEALLOCATE(w%t)
END PROGRAM thing_program
