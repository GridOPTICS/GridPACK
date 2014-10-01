! ----------------------------------------------------------------
! file: serial_io_f.F90
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Copyright (c) 2013 Battelle Memorial Institute
! Licensed under modified BSD License. A copy of this license can be found
! in the LICENSE file in the top level directory of this distribution.
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! Created September 17, 2014 by Bruce Palmer
! ----------------------------------------------------------------
!
!  Fortran serial IO functions
!
module gridpack_serial_io
  use, intrinsic :: iso_c_binding
  use gridpack_network
  implicit none
!
!  Define serial IO types
!
  private
!
  type, public :: bus_serial_io
    type(C_PTR) :: p_writer
    contains
    procedure::bus_create
    procedure::bus_destroy
    procedure::bus_open
    procedure::bus_close
    procedure::bus_write
    procedure::bus_header
  end type
!
  type, public :: branch_serial_io
    type(C_PTR) :: p_writer
    contains
    procedure::branch_create
    procedure::branch_destroy
    procedure::branch_open
    procedure::branch_close
    procedure::branch_write
    procedure::branch_header
  end type
  interface
!
! Create a bus serial IO object
! @param writer pointer to Fortran bus serial IO object
! @param max_len the maximum string length written by any bus
! @param network pointer to Fortran network object
!
    subroutine bus_serial_io_create(writer, max_len, network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: writer
      integer(C_INT), value, intent(in) :: max_len
      type(C_PTR), value, intent(in) :: network
    end subroutine bus_serial_io_create
!
! Destroy a bus serial IO object
! @param writer pointer to Fortran bus serial IO object
!
    subroutine bus_serial_io_destroy(writer) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: writer
    end subroutine bus_serial_io_destroy
!
! Open an external file and redirect output to it
! @param writer pointer to Fortran bus serial IO object
! @param filename name of new file
!
    subroutine bus_serial_io_open(writer, filename) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: writer
      character(C_CHAR), intent(in) :: filename(*)
    end subroutine bus_serial_io_open
!
! Close an external file and redirect output to standard out
! @param writer pointer to Fortran bus serial IO object
!
    subroutine bus_serial_io_close(writer) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: writer
    end subroutine bus_serial_io_close
!
! Write output from network to output
! @param writer pointer to Fortran bus serial IO object
!
    subroutine bus_serial_io_write(writer, signal) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: writer
      character(C_CHAR), intent(in) :: signal(*)
    end subroutine bus_serial_io_write
!
! Write single string to standard output. This is used to write headers for a
! data listing. It is mostly a convenience function so that users do not have
! to identify the head node
! @param writer pointer to Fortran bus serial IO object
! @param str character string containing the header
!
    subroutine bus_serial_io_header(writer, str) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: writer
      character(C_CHAR), intent(in) :: str(*)
    end subroutine bus_serial_io_header
!
! Create a branch serial IO object
! @param writer pointer to Fortran branch serial IO object
! @param max_len the maximum string length written by any branch
! @param network pointer to Fortran network object
!
    subroutine branch_serial_io_create(writer, max_len, network) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: writer
      integer(C_INT), value, intent(in) :: max_len
      type(C_PTR), value, intent(in) :: network
    end subroutine branch_serial_io_create
!
! Destroy a branch serial IO object
! @param writer pointer to Fortran branch serial IO object
!
    subroutine branch_serial_io_destroy(writer) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(inout) :: writer
    end subroutine branch_serial_io_destroy
!
! Open an external file and redirect output to it
! @param writer pointer to Fortran branch serial IO object
! @param filename name of new file
!
    subroutine branch_serial_io_open(writer, filename) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: writer
      character(C_CHAR), intent(in) :: filename(*)
    end subroutine branch_serial_io_open
!
! Close an external file and redirect output to standard out
! @param writer pointer to Fortran branch serial IO object
!
    subroutine branch_serial_io_close(writer) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: writer
    end subroutine branch_serial_io_close
!
! Write output from network to output
! @param writer pointer to Fortran branch serial IO object
!
    subroutine branch_serial_io_write(writer, signal) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: writer
      character(C_CHAR), intent(in) :: signal(*)
    end subroutine branch_serial_io_write
!
! Write single string to standard output. This is used to write headers for a
! data listing. It is mostly a convenience function so that users do not have
! to identify the head node
! @param writer pointer to Fortran branch serial IO object
! @param str character string containing the header
!
    subroutine branch_serial_io_header(writer, str) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value, intent(in) :: writer
      character(C_CHAR), intent(in) :: str(*)
    end subroutine branch_serial_io_header
  end interface
  contains
!
! Create a bus serial IO object
! @param p_writer pointer to Fortran bus serial IO object
! @param max_len the maximum string length written by any bus
! @param p_network pointer to Fortran network object
!
  subroutine bus_create(p_writer, max_len, p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_serial_io), intent(inout) :: p_writer
    integer, value, intent(in) :: max_len
    class(network), intent(in) :: p_network
    integer(C_INT) c_len
    c_len = max_len
    call bus_serial_io_create(p_writer%p_writer, c_len, p_network%p_network)
    return
  end subroutine bus_create
!
! Destroy a bus serial IO object
! @param p_writer pointer to Fortran bus serial IO object
!
  subroutine bus_destroy(p_writer)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_serial_io), intent(inout) :: p_writer
    call bus_serial_io_destroy(p_writer%p_writer)
    return
  end subroutine bus_destroy
!
! Open an external file and redirect output to it
! @param p_writer pointer to Fortran bus serial IO object
! @param filename name of new file
!
  subroutine bus_open(p_writer, filename)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_serial_io), intent(in) :: p_writer
    character, intent(in) :: filename(*)
    call bus_serial_io_open(p_writer%p_writer,filename)
    return
  end subroutine bus_open
!
! Close an external file and redirect output to standard out
! @param p_writer pointer to Fortran bus serial IO object
!
  subroutine bus_close(p_writer)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_serial_io), intent(in) :: p_writer
    call bus_serial_io_close(p_writer%p_writer)
    return
  end subroutine bus_close
!
! Write output from network to output
! @param p_writer pointer to Fortran bus serial IO object
!
  subroutine bus_write(p_writer, signal)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_serial_io), intent(in) :: p_writer
    character(len=*), intent(in) :: signal
    character(C_CHAR) c_string(512)
    integer slen, i
    slen = len(trim(signal))
    slen = min(slen,511)
    do i=1, slen
      c_string(i) = signal(i:i)
    end do
    c_string(slen+1) = C_NULL_CHAR
    call bus_serial_io_write(p_writer%p_writer, c_string)
    return
  end subroutine bus_write
!
! Write single string to standard output. This is used to write headers for a
! data listing. It is mostly a convenience function so that users do not have
! to identify the head node
! @param p_writer pointer to Fortran bus serial IO object
! @param str character string containing the header
!
  subroutine bus_header(p_writer, str)
    use, intrinsic :: iso_c_binding
    implicit none
    class(bus_serial_io), intent(in) :: p_writer
    character(len=*), intent(in) :: str
    character(C_CHAR) :: c_string(512)
    integer slen, i
    slen = len(trim(str))
    slen = min(slen,511)
    do i = 1, slen
      c_string(i) = str(i:i)
    end do
    c_string(slen+1) = C_NULL_CHAR
    call bus_serial_io_header(p_writer%p_writer,c_string)
    return
  end subroutine bus_header
!
! Create a branch serial IO object
! @param p_writer pointer to Fortran branch serial IO object
! @param max_len the maximum string length written by any bus
! @param p_network pointer to Fortran network object
!
  subroutine branch_create(p_writer, max_len, p_network)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_serial_io), intent(inout) :: p_writer
    integer, value, intent(in) :: max_len
    class(network), intent(in) :: p_network
    integer(C_INT) c_len
    c_len = max_len
    call branch_serial_io_create(p_writer%p_writer, max_len, p_network%p_network)
    return
  end subroutine branch_create
!
! Destroy a branch serial IO object
! @param p_writer pointer to Fortran branch serial IO object
!
  subroutine branch_destroy(p_writer)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_serial_io), intent(inout) :: p_writer
    call branch_serial_io_destroy(p_writer%p_writer)
    return
  end subroutine branch_destroy
!
! Open an external file and redirect output to it
! @param p_writer pointer to Fortran branch serial IO object
! @param filename name of new file
!
  subroutine branch_open(p_writer, filename)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_serial_io), intent(in) :: p_writer
    character, intent(in) :: filename(*)
    call branch_serial_io_open(p_writer%p_writer,filename)
    return
  end subroutine branch_open
!
! Close an external file and redirect output to standard out
! @param p_writer pointer to Fortran branch serial IO object
!
  subroutine branch_close(p_writer)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_serial_io), intent(in) :: p_writer
    call branch_serial_io_close(p_writer%p_writer)
    return
  end subroutine branch_close
!
! Write output from network to output
! @param p_writer pointer to Fortran branch serial IO object
!
  subroutine branch_write(p_writer, signal)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_serial_io), intent(in) :: p_writer
    character(len=*), intent(in) :: signal
    character(C_CHAR) c_string(512)
    integer slen, i
    slen = len(trim(signal))
    slen = min(slen,511)
    do i=1, slen
      c_string(i) = signal(i:i)
    end do
    c_string(slen+1) = C_NULL_CHAR
    call branch_serial_io_write(p_writer%p_writer, c_string)
    return
  end subroutine branch_write
!
! Write single string to standard output. This is used to write headers for a
! data listing. It is mostly a convenience function so that users do not have
! to identify the head node
! @param p_writer pointer to Fortran branch serial IO object
! @param str character string containing the header
!
  subroutine branch_header(p_writer, str)
    use, intrinsic :: iso_c_binding
    implicit none
    class(branch_serial_io), intent(in) :: p_writer
    character(len=*), intent(in) :: str
    character(C_CHAR) :: c_string(512)
    integer slen, i
    slen = len(trim(str))
    slen = min(slen,511)
    do i = 1, slen
      c_string(i) = str(i:i)
    end do
    c_string(slen+1) = C_NULL_CHAR
    call branch_serial_io_header(p_writer%p_writer,c_string)
    return
  end subroutine branch_header
end module gridpack_serial_io
