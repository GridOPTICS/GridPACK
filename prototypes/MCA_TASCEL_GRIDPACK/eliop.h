#ifndef  ELIOP_H
#define  ELIOP_H

#ifdef WIN32
#include <io.h>
#include "winutil.h"
#define F_OK 00
#endif

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>


#if 1

#define PRINT_AND_ABORT(msg, val) while(1){ fprintf(stderr," exiting in elio:%s %d\n",(msg),(val)); exit(2); }
#else

#define PRINT_AND_ABORT(msg, val) ga_error(msg, (long)val)
#ifndef GLOBAL_H
extern void ga_error(char*, long);
#endif

#endif


#if (defined(SP) || defined(SP1))
#define PIOFS 1
#endif


#if (defined(CRAY) && !defined(__crayx1)) || defined(NEC)
#        include <sys/statfs.h>
#        define  STATVFS statfs
#elif defined(KSR) || defined(__FreeBSD__) || defined(MACX)
#        include <sys/param.h>
#        include <sys/mount.h>
#        define  STATVFS statfs
#        define NO_F_FRSIZE 
#elif defined(WIN32)
#        define  STATVFS _stat 
#        define  S_ISDIR(mode) ((mode&S_IFMT) == S_IFDIR)
#        define  S_ISREG(mode) ((mode&S_IFMT) == S_IFREG)
#elif defined(CYGNUS) ||  defined(LINUX)  ||  defined(CYGWIN)
#        include <sys/vfs.h>
#        define  STATVFS statfs
#        define NO_F_FRSIZE 
#elif !defined(PARAGON)
#        include <sys/statvfs.h>
#        define  STATVFS statvfs
#endif

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <fcntl.h>
#if defined(PARAGON)
#  include <sys/mount.h>
#  include <nx.h>
#endif

#if (defined(CRAY) && defined(FFIO))
#        include <ffio.h>
#        include <sys/fstyp.h>
#        include <sys/fsid.h>
#endif


#include "elio.h"

extern int                   _elio_Errors_Fatal;
extern void                  elio_init(void);
extern int                   elio_pending_error;


#if !defined(PRINT_AND_ABORT)
#   if defined(SUN) && !defined(SOLARIS)
      extern int fprintf();
      extern void fflush();
#   endif
#   define PRINT_AND_ABORT(msg, val){\
     fprintf(stderr, "ELIO fatal error: %s %ld\n", msg,  val);\
     fprintf(stdout, "ELIO fatal error: %s %ld\n", msg,  val);\
     fflush(stdout);\
     exit(val);\
   }
#endif

/**************************** Error Macro ******************************/
/* ELIO defines error macro called in case of error
 * the macro can also use user-provided error routine PRINT_AND_ABORT
 * defined as macro to do some cleanup in the application before
 * aborting
 * The requirement is that PRINT_AND_ABORT is defined before
 * including ELIO header file - this file
 */

#define ELIO_ERROR_NULL(code, val){\
 PABLO_end(pablo_code);\
 if(! _elio_Errors_Fatal){\
     elio_pending_error= code;\
     return NULL;\
 }\
 if( _elio_Errors_Fatal)\
     PRINT_AND_ABORT(errtable[code-OFFSET], val);\
}

#define ELIO_ERROR(code, val) { \
 PABLO_end(pablo_code);\
 if(! _elio_Errors_Fatal) return(code);\
 else PRINT_AND_ABORT(errtable[code-OFFSET], (int)(val));\
}


/* error codes and messages */

#define ERRLEN 26
#define OFFSET    (-2000)
#define SEEKFAIL  (OFFSET + 0)
#define WRITFAIL  (OFFSET + 1)
#define AWRITFAIL (OFFSET + 2)
#define READFAIL  (OFFSET + 3)
#define AREADFAIL (OFFSET + 4)
#define SUSPFAIL  (OFFSET + 5)
#define HANDFAIL  (OFFSET + 6)
#define MODEFAIL  (OFFSET + 7)
#define DIRFAIL   (OFFSET + 8)
#define STATFAIL  (OFFSET + 9)
#define OPENFAIL  (OFFSET + 10)
#define ALOCFAIL  (OFFSET + 11)
#define UNSUPFAIL (OFFSET + 12)
#define DELFAIL   (OFFSET + 13)
#define CLOSFAIL  (OFFSET + 14)
#define INTRFAIL  (OFFSET + 15)
#define RETUFAIL  (OFFSET + 16)
#define LONGFAIL  (OFFSET + 17)
#define FTYPFAIL  (OFFSET + 18)
#define CONVFAIL  (OFFSET + 19)
#define TYPEFAIL  (OFFSET + 20)
#define PROBFAIL  (OFFSET + 21)
#define TRUNFAIL  (OFFSET + 22)
#define EOFFAIL   (OFFSET + 23)
#define FSYNCFAIL (OFFSET + 24)
#define UNKNFAIL  (OFFSET + 25)

extern  char *errtable[ERRLEN];

#define ELIO_FILENAME_MAX 1024
#define SDIRS_INIT_SIZE 1024

#  define PABLO_elio_write      710000
#  define PABLO_elio_awrite     710001
#  define PABLO_elio_read       710002
#  define PABLO_elio_aread      710003
#  define PABLO_elio_wait       710004
#  define PABLO_elio_probe      710005
#  define PABLO_elio_stat       710006
#  define PABLO_elio_open       710007
#  define PABLO_elio_gopen      710008
#  define PABLO_elio_close      710009
#  define PABLO_elio_set_cb     710010
#  define PABLO_elio_delete     710011
#  define PABLO_elio_truncate   710012
#  define PABLO_elio_length     710014

#  define PABLO_init 
#  define PABLO_start(_id) 
#  define PABLO_end( _id ) 
#  define PABLO_terminate  

#endif
