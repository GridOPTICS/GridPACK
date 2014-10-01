  ! ----------------------------------------------------------------
  ! file: nonlinear_solver_f.F90
  ! ----------------------------------------------------------------
  ! ----------------------------------------------------------------
  ! ----------------------------------------------------------------
  ! ----------------------------------------------------------------
  ! Created August 22, 2014 by William A. Perkins
  ! Last Change: 2014-10-01 08:01:58 d3g096
  ! ----------------------------------------------------------------
! ----------------------------------------------------------------
! MODULE gridpack_nonlinear_solver
! ----------------------------------------------------------------
MODULE gridpack_nonlinear_solver

  USE iso_c_binding
  USE gridpack_configuration
  USE gridpack_communicator
  USE gridpack_vector
  USE gridpack_matrix

  IMPLICIT NONE

  CHARACTER (LEN=80), PRIVATE, SAVE :: rcsid = "$Id$"

  PRIVATE

  ! ----------------------------------------------------------------
  ! TYPE builder
  ! This is a base type for types that define functions to fill a
  ! Jacobian matrix and a RHS vector describing a system of nonliner
  ! equations
  ! ----------------------------------------------------------------
  TYPE, ABSTRACT, PUBLIC :: builder
   CONTAINS
     PROCEDURE(jacobian_builder_method), DEFERRED :: jacobian
     PROCEDURE(function_builder_method), DEFERRED :: function
  END type builder

  ! Definitions of TYPE builder methods
  ABSTRACT INTERFACE 
     SUBROUTINE jacobian_builder_method(this, x, J)
       USE gridpack_vector
       USE gridpack_matrix
       IMPORT :: builder
       IMPLICIT NONE
       CLASS(builder), INTENT(IN) :: this
       TYPE (vector), INTENT(IN) :: x
       TYPE (matrix), INTENT(INOUT) :: J
     END SUBROUTINE jacobian_builder_method
     SUBROUTINE function_builder_method(this, x, F)
       USE gridpack_vector
       IMPORT :: builder
       IMPLICIT NONE
       CLASS(builder), INTENT(IN) :: this
       TYPE (vector), INTENT(IN) :: x
       TYPE (vector), INTENT(INOUT) :: F
     END SUBROUTINE function_builder_method
  END INTERFACE

  ! ----------------------------------------------------------------
  ! TYPE funcbuilder
  !
  ! This is a simple extension of type builder. The jacobian and
  ! function methods simply call user defined functions, compliant
  ! with the interfaces above.  No information is stored in the type.
  ! ----------------------------------------------------------------
  TYPE, PUBLIC, EXTENDS(builder) :: funcbuilder
     PROCEDURE (jacobian_builder_func), POINTER, NOPASS :: jfunc
     PROCEDURE (function_builder_func), POINTER, NOPASS :: ffunc
   CONTAINS
     PROCEDURE :: initialize => funcbuilder_initialize
     PROCEDURE :: jacobian => funcbuilder_jacobian
     PROCEDURE :: function => funcbuilder_function
  END type funcbuilder

  ! This defines the functions used by TYPE funcbuilder
  INTERFACE
     SUBROUTINE jacobian_builder_func(x, J)
       USE gridpack_vector
       USE gridpack_matrix
       IMPLICIT NONE
       TYPE (vector), INTENT(IN) :: x
       TYPE (matrix), INTENT(INOUT) :: J
     END SUBROUTINE jacobian_builder_func
     SUBROUTINE function_builder_func(x, F)
       USE gridpack_vector
       IMPLICIT NONE
       TYPE (vector), INTENT(IN) :: x
       TYPE (vector), INTENT(INOUT) :: F
     END SUBROUTINE function_builder_func
  END INTERFACE

  ! ----------------------------------------------------------------
  ! TYPE builder_wrap
  !
  ! A pointer to a builder instance cannot be passed to C and back and
  ! be useful.  So, this is a wrapper for builder.  A pointer to which
  ! can be passed back and forth to C.
  ! ----------------------------------------------------------------
  TYPE :: builder_wrap
     CLASS(builder), POINTER :: bldr
  END type builder_wrap


  ! ----------------------------------------------------------------
  ! TYPE nonlinear_solver
  ! The Fortran representation of NonlinearSolver instance.  
  ! ----------------------------------------------------------------
  TYPE, PUBLIC :: nonlinear_solver
     TYPE (c_ptr), PRIVATE :: impl
     TYPE (builder_wrap), POINTER, PRIVATE :: bldr
   CONTAINS
     PROCEDURE :: initialize => nonlinear_initialize
     PROCEDURE :: finalize      ! should be FINAL
     PROCEDURE :: solve
  END type nonlinear_solver

  ! ----------------------------------------------------------------
  ! 
  ! ----------------------------------------------------------------
  TYPE, EXTENDS(nonlinear_solver), PUBLIC :: newton_raphson_solver
     CONTAINS
       PROCEDURE :: initialize => newton_initialize
  END type newton_raphson_solver
     

  ! ----------------------------------------------------------------
  ! External C routines for nonlinear solver instances.
  ! ----------------------------------------------------------------
  INTERFACE
     FUNCTION nonlinear_solver_create(comm, conf, lsize, bldr) BIND(c)
       USE iso_c_binding, ONLY : c_ptr, c_int
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(IN), VALUE :: comm
       TYPE (c_ptr), INTENT(IN), VALUE :: conf
       INTEGER(c_int), INTENT(IN), VALUE :: lsize
       TYPE (c_ptr), INTENT(IN), VALUE :: bldr
       TYPE (c_ptr) :: nonlinear_solver_create
     END FUNCTION nonlinear_solver_create

     FUNCTION newton_solver_create(comm, conf, lsize, bldr) BIND(c)
       USE iso_c_binding, ONLY : c_ptr, c_int
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(IN), VALUE :: comm
       TYPE (c_ptr), INTENT(IN), VALUE :: conf
       INTEGER(c_int), INTENT(IN), VALUE :: lsize
       TYPE (c_ptr), INTENT(IN), VALUE :: bldr
       TYPE (c_ptr) :: newton_solver_create
     END FUNCTION newton_solver_create

     FUNCTION nonlinear_solver_jacobian(slvr) BIND(c)
       USE iso_c_binding, ONLY : c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(IN), VALUE :: slvr
       TYPE (c_ptr) :: nonlinear_solver_jacobian
     END FUNCTION nonlinear_solver_jacobian

     FUNCTION nonlinear_solver_rhs(slvr) BIND(c)
       USE iso_c_binding, ONLY : c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(IN), VALUE :: slvr
       TYPE (c_ptr) :: nonlinear_solver_rhs
     END FUNCTION nonlinear_solver_rhs

     SUBROUTINE nonlinear_solver_destroy(slvr) BIND(c)
       USE iso_c_binding, ONLY : c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(INOUT) :: slvr
     END SUBROUTINE nonlinear_solver_destroy

     SUBROUTINE nonlinear_solver_solve(slvr, x) BIND(c)
       USE iso_c_binding, ONLY : c_ptr
       IMPLICIT NONE
       TYPE (c_ptr), INTENT(IN), VALUE :: slvr
       TYPE (c_ptr), INTENT(IN), VALUE :: x
     END SUBROUTINE nonlinear_solver_solve
       
  END INTERFACE


CONTAINS

  ! ----------------------------------------------------------------
  ! SUBROUTINE funcbuilder_initialize
  ! ----------------------------------------------------------------
  SUBROUTINE funcbuilder_initialize(this, jproc, fproc)

    IMPLICIT NONE
    CLASS(funcbuilder), INTENT(INOUT) :: this
    PROCEDURE (jacobian_builder_func), POINTER :: jproc
    PROCEDURE (function_builder_func), POINTER :: fproc
    this%jfunc => jproc
    this%ffunc => fproc

  END SUBROUTINE funcbuilder_initialize

  ! ----------------------------------------------------------------
  ! SUBROUTINE funcbuilder_jacobian
  ! ----------------------------------------------------------------
  SUBROUTINE funcbuilder_jacobian(this, x, J)
    IMPLICIT NONE
    CLASS (funcbuilder), INTENT(IN) :: this
    TYPE (vector), INTENT(IN) :: x
    TYPE (matrix), INTENT(INOUT) :: J
    CALL this%jfunc(x, J)
  END SUBROUTINE funcbuilder_jacobian

  ! ----------------------------------------------------------------
  ! SUBROUTINE funcbuilder_function
  ! ----------------------------------------------------------------
  SUBROUTINE funcbuilder_function(this, x, F)
    IMPLICIT NONE
    CLASS (funcbuilder), INTENT(IN) :: this
    TYPE (vector), INTENT(IN) :: x
    TYPE (vector), INTENT(INOUT) :: F
    CALL this%ffunc(x, F)
  END SUBROUTINE funcbuilder_function


  ! ----------------------------------------------------------------
  ! SUBROUTINE builder_build_jacobian
  !
  ! Called from C. Extracts the builder instance from its wrapper and
  ! calls its jacobian method.
  ! ----------------------------------------------------------------
  SUBROUTINE builder_build_jacobian(bptr, jptr, xptr) BIND(c)

    IMPLICIT NONE
    TYPE (c_ptr), INTENT(IN), VALUE :: bptr, xptr, jptr
    TYPE (builder_wrap), POINTER :: bwrap
    TYPE (matrix) :: J
    TYPE (vector) :: X

    CALL C_F_POINTER(bptr, bwrap)
    ! Note: matix/vector initalize/finalize need to be avoided here
    J%mat = jptr
    X%vec = xptr

    CALL bwrap%bldr%jacobian(X, J)

  END SUBROUTINE builder_build_jacobian

  ! ----------------------------------------------------------------
  ! SUBROUTINE builder_build_function
  !
  ! Called from C. Extracts the builder instance from its wrapper and
  ! calls its function method.
  ! ----------------------------------------------------------------
  SUBROUTINE builder_build_function(bptr, fptr, xptr) BIND(c)
    
    IMPLICIT NONE
    TYPE (c_ptr), INTENT(IN), VALUE :: bptr, xptr, fptr
    TYPE (builder_wrap), POINTER :: bwrap
    TYPE (vector) :: F, X

    CALL C_F_POINTER(bptr, bwrap)

    ! Note: vector initalize/finalize need to be avoided here
    F%vec = fptr
    X%vec = xptr

    CALL bwrap%bldr%function(X, F)

  END SUBROUTINE builder_build_function

  ! ----------------------------------------------------------------
  ! SUBROUTINE nonlinear_initialize
  ! ----------------------------------------------------------------
  SUBROUTINE nonlinear_initialize(this, comm, conf, local_size, bldr)
    IMPLICIT NONE
    CLASS(nonlinear_solver), INTENT(INOUT) :: this
    CLASS(communicator), INTENT(IN) :: comm
    CLASS(cursor), INTENT(IN) :: conf
    INTEGER, INTENT(IN) :: local_size
    CLASS(builder), POINTER, INTENT(IN) :: bldr
    INTEGER(c_int) :: csize
    TYPE (c_ptr) :: bldrloc
    csize = local_size
    ALLOCATE(this%bldr)
    this%bldr%bldr => bldr
    bldrloc = C_LOC(this%bldr)
    this%impl = nonlinear_solver_create(comm%comm, conf%impl, csize, bldrloc)
  END SUBROUTINE nonlinear_initialize


  SUBROUTINE finalize(this)
    IMPLICIT NONE
    CLASS(nonlinear_solver), INTENT(INOUT) :: this
    DEALLOCATE(this%bldr)
    CALL nonlinear_solver_destroy(this%impl)
  END SUBROUTINE finalize

  SUBROUTINE solve(this, x)
    IMPLICIT NONE
    CLASS(nonlinear_solver), INTENT(INOUT) :: this
    CLASS(vector), INTENT(IN) :: x
    CALL nonlinear_solver_solve(this%impl, x%vec)
  END SUBROUTINE solve

  ! ----------------------------------------------------------------
  ! SUBROUTINE newton_initialize
  ! ----------------------------------------------------------------
  SUBROUTINE newton_initialize(this, comm, conf, local_size, bldr)
    IMPLICIT NONE
    CLASS(newton_raphson_solver), INTENT(INOUT) :: this
    CLASS(communicator), INTENT(IN) :: comm
    CLASS(cursor), INTENT(IN) :: conf
    INTEGER, INTENT(IN) :: local_size
    CLASS(builder), POINTER, INTENT(IN) :: bldr
    INTEGER(c_int) :: csize
    TYPE (c_ptr) :: bldrloc
    csize = local_size
    ALLOCATE(this%bldr)
    this%bldr%bldr => bldr
    bldrloc = C_LOC(this%bldr)
    this%impl = newton_solver_create(comm%comm, conf%impl, csize, bldrloc)
  END SUBROUTINE newton_initialize

END MODULE gridpack_nonlinear_solver

