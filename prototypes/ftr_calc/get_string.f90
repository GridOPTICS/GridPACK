! Pass in string and break into tokens separted with spaces
! Count on caller giving us entries in the table
subroutine get_strings ( string, tokens, sublen)
    implicit none
    character(len=*), intent(in)                  :: string
    character(len=80), dimension(*), intent(out)  :: tokens
    integer, dimension(*), intent(out)            :: sublen
    integer posit
    integer num_strings    
    integer length
    integer substr

    length = len(string)
    posit  = 1
    substr = 1

    do while (posit .le. length)
       tokens(substr) = get_string(posit, string, sublen(substr))
       substr = substr + 1
    enddo
    sublen(substr) =  0

    return

contains

! First time it is called on a string, call with posit = 1
function get_string( posit, string, length )
    implicit none
    integer, intent(inout)         :: posit
    character(len=*), intent(in)   :: string
    integer, intent(out)           :: length

    character(len=len(string))     :: get_string

    integer                        :: lenstr
    integer                        :: str_start
    integer                        :: str_end

    lenstr      = len(string)

    do while ( posit .le. lenstr )
       if ( ' ' .eq. string(posit:posit) ) then
          posit = posit + 1
       else
          str_start = posit
          goto 99
       endif
    enddo
    get_string = ' '
    length     = 0 
    return ! Found nothing

 99 continue ! We are now passed all the spaces

    str_end = posit ! Just in case single character at end
    do while ( posit .le. lenstr )
       if ( ' ' .ne. string(posit:posit) ) then
          posit = posit + 1
       else
          str_end = posit - 1
          goto 88
       endif
    enddo
 88 continue
    length     = str_end-str_start+1
    get_string = string(str_start:str_end)
    return

end function get_string

end subroutine get_strings
