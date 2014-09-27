! ----------------------------------------------------------------
! file: configuration_test.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created August 13, 2014 by William A. Perkins
! Last Change: 2014-08-22 08:22:18 d3g096
! ----------------------------------------------------------------

PROGRAM configuration_test
  USE gridpack_parallel
  USE gridpack_configuration

  IMPLICIT NONE

  TYPE (communicator) :: comm
  TYPE (cursor) :: conf, cur1, cur2
  INTEGER(c_int) :: ma_stack = 200000, ma_heap = 200000
  LOGICAL :: flag
  DOUBLE PRECISION :: rvalue
  INTEGER :: ivalue
  CHARACTER(128) :: string

  CALL gridpack_initialize_parallel(ma_stack, ma_heap)
  CALL comm%initialize()
  
  CALL conf%initialize()
  IF (conf%ok()) THEN
     WRITE(*,*) 'Huh'
  END IF
  IF (.NOT. conf%open(comm, "gridpack.xml")) THEN
     WRITE (*, *) 'Cannot open configuration'
     STOP
  END IF

  cur1 = conf%get_cursor("Bogus")
  IF (cur1%ok()) THEN 
     WRITE(*,*) 'ERROR: found Bogus node'
  END IF
  CALL cur1%finalize()

  cur1 = conf%get_cursor("GridPACK")
  IF (.NOT. cur1%ok()) THEN
     WRITE(*,*) 'ERROR: did not find GridPACK node'
  END IF

  IF (cur1%get_bool("Boolean", flag)) THEN
     IF (.NOT. flag) THEN
        WRITE (*,*) 'ERROR: wrong value for GridPACK.Boolean'
     END IF
  ELSE
     WRITE(*,*) 'ERROR: did not find GridPACK.Boolean'
  END IF

  IF (cur1%get_int("Integer", ivalue)) THEN
     IF (ivalue .NE. 5) THEN
        WRITE (*,*) 'ERROR: wrong value for GridPACK.Integer: ', ivalue
     END IF
  ELSE 
     WRITE(*,*) 'ERROR: did not find GridPACK.Integer'
  END IF

  IF (cur1%get_double("Real", rvalue)) THEN
     IF (ABS(rvalue - 1.234D05) .GT. 1.0E-06) THEN
        WRITE (*,*) 'ERROR: wrong value for GridPACK.Real: ', rvalue
     END IF
  ELSE 
     WRITE(*,*) 'ERROR: did not find GridPACK.Real'
  END IF

  IF (cur1%get_string("String", string)) THEN
     IF (trim(string).ne.'Dog') THEN
        WRITE (*,*) 'ERROR: wrong value for GridPACK.String: ', trim(string)
     END IF
  ELSE 
     WRITE(*,*) 'ERROR: did not find GridPACK.String'
  END IF

  cur2 = conf%get_cursor("GridPACK.Subpath")
  IF (.NOT. cur2%ok()) THEN
     WRITE(*,*) 'ERROR: cur2: did not find GridPACK.Subpath'
  END IF

  IF (cur2%get_bool("Boolean", flag)) THEN
     IF (flag) THEN
        WRITE (*,*) 'ERROR: wrong value for GridPACK.Subpath.Boolean'
     END IF
  ELSE
     WRITE(*,*) 'ERROR: did not find GridPACK.Subpath.Boolean'
  END IF

  IF (cur2%get_int("Integer", ivalue)) THEN
     IF (ivalue .NE. 3) THEN
        WRITE (*,*) 'ERROR: wrong value for GridPACK.Subpath.Integer: ', ivalue
     END IF
  ELSE 
     WRITE(*,*) 'ERROR: did not find GridPACK.Subpath.Integer'
  END IF

  IF (cur2%get_double("Real", rvalue)) THEN
     IF (ABS(rvalue - 1.23D00) .GT. 1.0E-06) THEN
        WRITE (*,*) 'ERROR: wrong value for GridPACK.Subpath.Real: ', rvalue
     END IF
  ELSE 
     WRITE(*,*) 'ERROR: did not find GridPACK.Subpath.Real'
  END IF

  IF (cur2%get_string("String", string)) THEN
     IF (trim(string).ne.'Hamster') THEN
        WRITE (*,*) 'ERROR: wrong value for GridPACK.String: ', trim(string)
     END IF
  ELSE 
     WRITE(*,*) 'ERROR: did not find GridPACK.Subpath.String'
  END IF

  CALL conf%set_path("GridPACK.Subpath")
  IF (.NOT. conf%ok()) THEN
     WRITE(*,*) 'ERROR: conf: did not find GridPACK.Subpath'
  END IF

  IF (conf%get_bool("Boolean", flag)) THEN
     IF (flag) THEN
        WRITE (*,*) 'ERROR: wrong value for GridPACK.Subpath.Boolean'
     END IF
  ELSE
     WRITE(*,*) 'ERROR: did not find GridPACK.Subpath.Boolean'
  END IF

  IF (conf%get_int("Integer", ivalue)) THEN
     IF (ivalue .NE. 3) THEN
        WRITE (*,*) 'ERROR: wrong value for GridPACK.Subpath.Integer: ', ivalue
     END IF
  ELSE 
     WRITE(*,*) 'ERROR: did not find GridPACK.Subpath.Integer'
  END IF

  IF (conf%get_double("Real", rvalue)) THEN
     IF (ABS(rvalue - 1.23D00) .GT. 1.0E-06) THEN
        WRITE (*,*) 'ERROR: wrong value for GridPACK.Subpath.Real: ', rvalue
     END IF
  ELSE 
     WRITE(*,*) 'ERROR: did not find GridPACK.Subpath.Real'
  END IF

  IF (conf%get_string("String", string)) THEN
     IF (trim(string).ne.'Hamster') THEN
        WRITE (*,*) 'ERROR: wrong value for GridPACK.String: ', trim(string)
     END IF
  ELSE 
     WRITE(*,*) 'ERROR: did not find GridPACK.Subpath.String'
  END IF

  CALL cur1%finalize()
  CALL cur2%finalize()

  CALL conf%finalize()
  CALL comm%finalize()
  CALL gridpack_finalize_parallel()

END PROGRAM configuration_test
