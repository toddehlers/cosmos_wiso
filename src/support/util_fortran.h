#include <config.h>

/* Should provide the information that this functions are f90 interfaced */
 
#define FORTRAN_CALL

#if defined(HAVE_FORTRAN_H) && defined(CRAY)
#include <fortran.h>
#else
#define _fcd        char *
#define _fcdtocp(x) x
#endif

#if SIZEOF_LONG == 8
typedef unsigned long ULONG ;
typedef long LONG; 
#elif SIZEOF_LONG_LONG == 8 
typedef unsigned long long ULONG ;
typedef long long LONG;
#endif

/* INT should be native Fortran INTEGER */

#if FORT_INT == 6
typedef int INT;
#elif FORT_INT == 8
typedef long INT;
#elif FORT_INT == 13
typedef long long INT;
#endif

#if FORT_REAL == 10
typedef float REAL;
#elif FORT_REAL == 11
typedef double REAL; 
#endif

