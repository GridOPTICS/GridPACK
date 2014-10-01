! ----------------------------------------------------------------
! file: mathf.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created May 15, 2014 by William A. Perkins
! Last Change: 2014-05-29 13:29:02 d3g096
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! MODULE gridpack_math
! ----------------------------------------------------------------
MODULE gridpack_math

  USE iso_c_binding
  USE gridpack_vector
  USE gridpack_matrix
  USE gridpack_linear_solver

  IMPLICIT NONE

  CHARACTER (LEN=80), PRIVATE, SAVE :: rcsid = "$Id$"

  INTERFACE

     SUBROUTINE gridpack_initialize_math() BIND(c)
     END SUBROUTINE gridpack_initialize_math

     LOGICAL(c_bool) FUNCTION gridpack_math_initialized() BIND(c)
       USE iso_c_binding, ONLY: c_bool
     END FUNCTION gridpack_math_initialized

     SUBROUTINE gridpack_finalize_math() BIND(c)
     END SUBROUTINE gridpack_finalize_math

  END INTERFACE

END MODULE gridpack_math
