
/* Generated automatically from m214003 on Mon Aug 16 13:55:18 CEST 2004 */

#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <sys/types.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define HAVE_CF_INTERFACE 1

#define RENAME_PB         1

#ifndef _GRIB_H
#define _GRIB_H

#include <stdio.h>
#include <sys/types.h>

#define  GRIB_MISSVAL  -9.E33

/* Level Types */
#define  LTYPE_SURFACE               1
#define  LTYPE_99                   99
#define  LTYPE_ISOBARIC            100
#define  LTYPE_MEANSEA             102
#define  LTYPE_ALTITUDE            103
#define  LTYPE_HEIGHT              105
#define  LTYPE_SIGMA               107
#define  LTYPE_HYBRID              109
#define  LTYPE_DOWN                111
#define  LTYPE_GROUND              112
#define  LTYPE_DEPTH_BELOW_SEA     160
#define  LTYPE_99_MARGIN          1000

/*
 *  Data representation type (Grid Type) [Table 6]
 */
#define  GTYPE_LATLON                0  /*  latitude/longitude                       */
#define  GTYPE_LATLON_ROT           10  /*  rotated latitude/longitude               */
#define  GTYPE_LATLON_STR           20  /*  stretched latitude/longitude             */
#define  GTYPE_LATLON_ROTSTR        30  /*  rotated and stretched latitude/longitude */
#define  GTYPE_GAUSSIAN              4  /*  gaussian grid                            */
#define  GTYPE_GAUSSIAN_ROT         14  /*  rotated gaussian grid                    */
#define  GTYPE_GAUSSIAN_STR         24  /*  stretched gaussian grid                  */
#define  GTYPE_GAUSSIAN_ROTSTR      34  /*  rotated and stretched gaussian grid      */
#define  GTYPE_SPECTRAL             50  /*  spherical harmonics                      */
#define  GTYPE_TRIANGULAR          192  /*  triangular grid                          */

/*
 *  Macros for the indicator section ( Section 0 )
 */
#define  ISEC0_GRIB_Len             (isec0[ 0])  /*  Number of octets in the GRIB message         */
#define  ISEC0_GRIB_Version         (isec0[ 1])  /*  GRIB edition number                          */


/*
 *  Macros for the product definition section ( Section 1 )
 */
#define  ISEC1_TABLE4_MINUTE  0
#define  ISEC1_TABLE4_HOUR    1
#define  ISEC1_TABLE4_DAY     2


#define  ISEC1_CodeTable            (isec1[ 0])  /*  Version number of code table                 */
#define  ISEC1_CenterID             (isec1[ 1])  /*  Identification of centre                     */
#define  ISEC1_ModelID              (isec1[ 2])  /*  Identification of model                      */
#define  ISEC1_GridDefinition       (isec1[ 3])  /*  Grid definition                              */
#define  ISEC1_Sec2Or3Flag          (isec1[ 4])  /*  Section 2 or 3 included                      */
#define  ISEC1_Parameter            (isec1[ 5])  /*  Parameter indicator                          */
#define  ISEC1_LevelType            (isec1[ 6])  /*  Type of level indicator                      */
#define  ISEC1_Level1               (isec1[ 7])  /*  Level 1                                      */
#define  ISEC1_Level2               (isec1[ 8])  /*  Level 2                                      */
#define  ISEC1_Year                 (isec1[ 9])  /*  Year of century (YY)                         */
#define  ISEC1_Month                (isec1[10])  /*  Month (MM)                                   */
#define  ISEC1_Day                  (isec1[11])  /*  Day (DD)                                     */
#define  ISEC1_Hour                 (isec1[12])  /*  Hour (HH)                                    */
#define  ISEC1_Minute               (isec1[13])  /*  Minute (MM)                                  */
#define  ISEC1_TimeUnit             (isec1[14])  /*  Time unit indicator                          */
#define  ISEC1_TimePeriod1          (isec1[15])  /*  P1 Time period                               */
#define  ISEC1_TimePeriod2          (isec1[16])  /*  P2 Time period                               */
#define  ISEC1_TimeRange            (isec1[17])  /*  Time range indicator                         */
#define  ISEC1_AvgNum               (isec1[18])  /*  Number of products included in an average    */
#define  ISEC1_AvgMiss              (isec1[19])  /*  Number of products missing form an average   */
#define  ISEC1_Century              (isec1[20])  /*  Century                                      */
#define  ISEC1_SubCenterID          (isec1[21])  /*  Subcenter identifier                         */
#define  ISEC1_DecScaleFactor       (isec1[22])  /*  Decimal scale factor                         */
#define  ISEC1_LocalFLag            (isec1[23])  /*  Flag field to indicate local use in isec1    */

#define  ISEC1_ECMWF_LocalExtention (isec1[36])
#define  ISEC1_ECMWF_Class          (isec1[37])


/*
 *  Macros for the grid definition section ( Section 2 )
 */
#define  ISEC2_GridType             (isec2[ 0])  /* Data representation type */

/* Triangular grids */

#define  ISEC2_TRI_NI2              (isec2[ 1])  /*  Number of factor 2 in factorisation of Ni    */
#define  ISEC2_TRI_NI3              (isec2[ 2])  /*  Number of factor 3 in factorisation of Ni    */
#define  ISEC2_TRI_ND               (isec2[ 3])  /*  Nubmer of diamonds                           */
#define  ISEC2_TRI_NI               (isec2[ 4])  /*  Number of tri. subdiv. of the icosahedron    */
#define  ISEC2_TRI_AFlag            (isec2[ 5])  /*  Flag for orientation of diamonds (Table A)   */
#define  ISEC2_TRI_LatPP            (isec2[ 6])  /*  Latitude of pole point                       */
#define  ISEC2_TRI_LonPP            (isec2[ 7])  /*  Longitude of pole point                      */
#define  ISEC2_TRI_LonMPL           (isec2[ 8])  /*  Longitude of the first diamond               */
#define  ISEC2_TRI_BFlag            (isec2[ 9])  /*  Flag for storage sequence (Table B)          */

/* Spherical harmonic coeficients */

#define  ISEC2_PentaJ               (isec2[ 1])  /*  J pentagonal resolution parameter            */
#define  ISEC2_PentaK               (isec2[ 2])  /*  K pentagonal resolution parameter            */
#define  ISEC2_PentaM               (isec2[ 3])  /*  M pentagonal resolution parameter            */
#define  ISEC2_RepType              (isec2[ 4])  /*  Representation type                          */
#define  ISEC2_RepMode              (isec2[ 5])  /*  Representation mode                          */

/* Gaussian grids */

#define  ISEC2_NumLon               (isec2[ 1])  /*  Number of points along a parallel (Ni)       */
#define  ISEC2_NumLat               (isec2[ 2])  /*  Number of points along a meridian (Nj)       */
#define  ISEC2_FirstLat             (isec2[ 3])  /*  Latitude of the first grid point             */
#define  ISEC2_FirstLon             (isec2[ 4])  /*  Longitude of the first grid point            */
#define  ISEC2_ResFlag              (isec2[ 5])  /*  Resolution flag: 128 regular grid            */
#define  ISEC2_LastLat              (isec2[ 6])  /*  Latitude of the last grid point              */
#define  ISEC2_LastLon              (isec2[ 7])  /*  Longitude of the last grid point             */
#define  ISEC2_LonIncr              (isec2[ 8])  /*  i direction increment                        */
#define  ISEC2_LatIncr              (isec2[ 9])  /*  j direction increment                        */
#define  ISEC2_NumPar               (isec2[ 9])  /*  Number of parallels between a pole and the E.*/
#define  ISEC2_ScanFlag             (isec2[10])  /*  Scanning mode flags                          */
#define  ISEC2_NumVCP               (isec2[11])  /*  Number of vertical coordinate parameters     */

#define  ISEC2_Reduced              (isec2[16])  /* 0: regular, 1: reduced grid                   */

#define  ISEC2_RowLonPtr            (&isec2[22])
#define  ISEC2_RowLon(i)            (isec2[22+i]) /* Number of points along each parallel         */

/* */

#define  ISEC2_LatSP                (isec2[12])  /* Latitude of the southern pole of rotation     */
#define  ISEC2_LonSP                (isec2[13])  /* Longitude of the southern pole of rotation    */

#define  FSEC2_RotAngle             (fsec2[ 0])  /* Angle of rotation                             */
#define  FSEC2_StrFact              (fsec2[ 1])  /* Stretching factor                             */

/*
 *  Macros for the bit map section ( Section 3 )
 */
#define  ISEC3_PredefBitmap         (isec3[ 0])  /* Predefined bitmap                             */
#define  ISEC3_MissVal              (isec3[ 1])  /* Missing data value for integers               */
#define  FSEC3_MissVal              (fsec3[ 1])  /* Missing data value for floats                 */

/*
 *  Macros for the binary data section ( Section 4 )
 */
#define  ISEC4_NumValues            (isec4[ 0])  /* Number of data values for encode/decode       */
#define  ISEC4_NumBits              (isec4[ 1])  /* Number of bits used for each encoded value    */
#define  ISEC4_NumNonMissValues     (isec4[20])  /* Number of non-missing values                  */




void  gribSetDebug(int debug);
void  gribSetRound(int round);
void  gribSetRefDP(double refval);
void  gribSetRefSP(float  refval);
void  gribSetValueCheck(int vcheck);


void  gribExSP(int *isec0, int *isec1, int *isec2, float *fsec2, int *isec3,
               float *fsec3, int *isec4, float *fsec4, int klenp, int *kgrib,
               int kleng, int *kword, char *hoper, int *kret);

void  gribExDP(int *isec0, int *isec1, int *isec2, double *fsec2, int *isec3,
               double *fsec3, int *isec4, double *fsec4, int klenp, int *kgrib,
               int kleng, int *kword, char *hoper, int *kret);


const char *gribLibraryVersion(void);

void  gribDebug(int debug);


void  gribDateTime(int *isec1, int *date, int *time);
int   gribRefDate(int *isec1);
int   gribRefTime(int *isec1);
int   gribTimeIsFC(int *isec1);

void  gribPrintSec0(int *isec0);
void  gribPrintSec1(int *isec0, int *isec1);
void  gribPrintSec2DP(int *isec0, int *isec2, double *fsec2);
void  gribPrintSec2SP(int *isec0, int *isec2, float  *fsec2);
void  gribPrintSec3DP(int *isec0, int *isec3, double *fsec3);
void  gribPrintSec3SP(int *isec0, int *isec3, float  *fsec3);
void  gribPrintSec4DP(int *isec0, int *isec4, double *fsec4);
void  gribPrintSec4SP(int *isec0, int *isec4, float  *fsec4);
void  gribPrintSec4Wave(int *isec4);

void  gribPrintALL(int nrec, long offset, long recpos, long recsize, unsigned char *gribbuffer);
void  gribPrintPDS(int nrec, long recpos, long recsize, unsigned char *gribbuffer);
void  gribPrintGDS(int nrec, long recpos, long recsize, unsigned char *gribbuffer);
void  gribPrintBMS(int nrec, long recpos, long recsize, unsigned char *gribbuffer);
void  gribPrintBDS(int nrec, long recpos, long recsize, unsigned char *gribbuffer);


int   gribOpen(const char *filename, const char *mode);
void  gribClose(int fileID);

int   gribRead(int fileID, unsigned char *buffer, size_t *buffersize);
int   gribWrite(int fileID, unsigned char *buffer, size_t buffersize);
off_t gribGetPos(int fileID);
int   gribGetSize(int fileID);
int   gribFileSeek(int fileID, long *offset);
int   gribReadSize(int fileID);
int   gribVersion(unsigned char *buffer, size_t buffersize);


#endif  /* _GRIB_H */
#ifndef _DTYPES_H
#define _DTYPES_H

#include <stdio.h>
#include <limits.h>

/* INT32 */

#if ! defined (INT_MAX)
#  error INT_MAX undefined
#endif

#undef  INT32
#if  INT_MAX == 2147483647L
#  define  INT32  int
#elif LONG_MAX == 2147483647L
#  define  INT32  long
#endif

/* INT64 */

#if ! defined (LONG_MAX)
#  error LONG_MAX undefined
#endif

#undef  INT64
#if  LONG_MAX > 2147483647L
#  define  INT64  long
#else
#  define  INT64  long long
#endif

/* REAL4 */

#undef   REAL4
#define  REAL4  float

/* Real8 */

#undef   REAL8
#define  REAL8  double

#endif  /* _DTYPES_H */
#ifndef _ERROR_H
#define _ERROR_H

#define  _FATAL     1     /* Error flag: exit on error  */
#define  _VERBOSE   2     /* Error flag: report errors  */
#define  _DEBUG     4     /* Error flag: debug          */

extern int _ExitOnError;  /* If set to 1, exit on error (default 1)       */
extern int _Verbose;      /* If set to 1, errors are reported (default 1) */
extern int _Debug;        /* If set to 1, debuggig (default 0)            */

void SysError(const char *caller, const char *fmt, ...);
void    Error(const char *caller, const char *fmt, ...);
void  Warning(const char *caller, const char *fmt, ...);
void  Message(const char *caller, const char *fmt, ...);

#endif  /* _ERROR_H */
#ifndef _CALENDAR_H
#define _CALENDAR_H

typedef struct {
  double  origin;  /* origin datatype */
  double  factor;  /* scaling factor  */
}
sUnit;

void Calendar(double value, sUnit unit, int *year, int *month,
	      int *day, int *hour, int *minute, double *second);
void InvCalendar(int year, int month, int day, int hour,
		 int minute, double second, sUnit unit, double *value);

#endif  /* _CALENDAR_H */
#ifndef _GRIB_INT_H
#define _GRIB_INT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


#if ! defined   (_GRIB_H)
#  include "grib.h"
#endif
#if ! defined   (_ERROR_H)
#  include "error.h"
#endif
#if ! defined   (_DTYPES_H)
#  include "dtypes.h"
#endif

#if ! defined   (FALSE)
#  define  FALSE  0
#endif

#if ! defined   (TRUE)
#  define  TRUE  1
#endif

#if ! defined   (UCHAR)
#  define  UCHAR  unsigned char
#endif

#if  defined  (INT32)
#  define  GRIBPACK     INT32
#  define  PACK_GRIB    packInt32
#  define  UNPACK_GRIB  unpackInt32
#else
#  define  GRIBPACK     INT64
#  define  PACK_GRIB    packInt64
#  define  UNPACK_GRIB  unpackInt64
#endif


extern FILE *grprsm;

extern int  GRB_Debug;

void   gprintf(const char *caller, const char *fmt, ...);

void   grsdef(void);

void   prtbin(int kin, int knbit, int *kout, int *kerr);
void   confp3(double pval, int *kexp, int *kmant, int kbits, int kround);
double decfp2(int kexp, int kmant);
void   ref2ibm(double *pref, int kbits);

void   scaleComplex(double *fpdata, int pcStart, int pcScale, int truncation);
void   scatterComplex(double *fpdata, int pcStart, int truncation, int dimSP);
void   scm0(double *pdl, double *pdr, double *pfl, double *pfr, int klg);
int    rowina2(double *p, int ko, int ki, double *pw,
	       int kcode, double msval, int *kret);
int    qu2reg2(double *pfield, int *kpoint, int klat, int klon,
	       double *ztemp, double msval, int *kret);

#if  defined  (INT32)
long   packInt32(INT32 *up, char *cp, long bc, long tc);
#endif
long   packInt64(INT64 *up, char *cp, long bc, long tc);
#if  defined  (INT32)
long   unpackInt32(unsigned char *cp, INT32 *up, long bc, long tc);
#endif
long   unpackInt64(unsigned char *cp, INT64 *up, long bc, long tc);

void  gribEncode(int *isec0, int *isec1, int *isec2, double *fsec2, int *isec3,
		 double *fsec3, int *isec4, double *fsec4, int klenp, int *kgrib,
		 int kleng, int *kword, int efunc, int *kret);

void  gribDecode(int *isec0, int *isec1, int *isec2, double *fsec2, int *isec3,
		 double *fsec3, int *isec4, double *fsec4, int klenp, int *kgrib,
		 int kleng, int *kword, int dfunc, int *kret);

#endif  /* _GRIB_INT_H */
#ifndef _GRIBDECODE_H
#define _GRIBDECODE_H

#define  UNDEFINED          9.999e20


#define  GET_INT3(a,b,c)    ((1-(int) ((unsigned) (a & 128) >> 6)) * (int) (((a & 127) << 16)+(b<<8)+c))
#define  GET_INT2(a,b)      ((1-(int) ((unsigned) (a & 128) >> 6)) * (int) (((a & 127) << 8) + b))

/* this requires a 32-bit default integer machine */
#define  GET_UINT4(a,b,c,d) ((int) ((a << 24) + (b << 16) + (c << 8) + (d)))
#define  GET_UINT3(a,b,c)   ((int) ((a << 16) + (b << 8)  + (c)))
#define  GET_UINT2(a,b)     ((int) ((a << 8)  + (b)))
#define  GET_UINT1(a)       ((int)  (a))


/* Section 0: Indicator Section (IS) */

#define  IS_GRIB_Len         GET_INT3(is[ 4], is[ 5], is[ 6])
#define  IS_GRIB_Version     GET_UINT1(is[ 7])

/* Section 1: Product Definition Section (PDS) */

#define  PDS_Len             GET_UINT3(pds[ 0], pds[ 1], pds[ 2])
#define  PDS_CodeTable       GET_UINT1(pds[ 3])
#define  PDS_CenterID        GET_UINT1(pds[ 4])
#define  PDS_ModelID         GET_UINT1(pds[ 5])
#define  PDS_GridDefinition  GET_UINT1(pds[ 6])
#define  PDS_Sec2Or3Flag     GET_UINT1(pds[ 7])
#define  PDS_HAS_GDS         ((pds[7] & 128) != 0)
#define  PDS_HAS_BMS         ((pds[7] &  64) != 0)
#define  PDS_Parameter       GET_UINT1(pds[ 8])
#define  PDS_LevelType       GET_UINT1(pds[ 9])
#define  PDS_Level1          (pds[10])
#define  PDS_Level2	     (pds[11])
#define  PDS_Level	     GET_UINT2(pds[10], pds[11])
#define  PDS_Year            GET_UINT1(pds[12])
#define  PDS_Month           GET_UINT1(pds[13])
#define  PDS_Day             GET_UINT1(pds[14])
#define  PDS_Hour            GET_UINT1(pds[15])
#define  PDS_Minute          GET_UINT1(pds[16])
#define  PDS_Date            (PDS_Year*10000+PDS_Month*100+PDS_Day)
#define  PDS_Time            (PDS_Hour*100+PDS_Minute)
#define  PDS_TimeUnit        GET_UINT1(pds[17])
#define  PDS_TimePeriod1     GET_UINT1(pds[18])
#define  PDS_TimePeriod2     GET_UINT1(pds[19])
#define  PDS_TimeRange       GET_UINT1(pds[20])
#define  PDS_AvgNum          GET_UINT2(pds[21], pds[22])
#define  PDS_AvgMiss         GET_UINT1(pds[23])
#define  PDS_Century         GET_UINT1(pds[24])
#define  PDS_Subcenter       GET_UINT1(pds[25])
#define  PDS_DecimalScale    GET_INT2(pds[26],pds[27])


/* Section 2: Grid Description Section (GDS) */

#define  GDS_Len             ((gds) == NULL ? 0 : GET_UINT3(gds[ 0], gds[ 1], gds[ 2]))
#define  GDS_NV              GET_UINT1(gds[ 3])
#define  GDS_PVPL            GET_UINT1(gds[ 4])
#define  GDS_PV	             ((gds[3] == 0) ? -1 : (int) gds[4] - 1)
#define  GDS_PL	             ((gds[4] == 255) ? -1 : (int) gds[3] * 4 + (int) gds[4] - 1)
#define  GDS_GridType        GET_UINT1(gds[ 5])


/* Triangular grid of DWD */
#define  GDS_TRI_NI2         GET_UINT2(gds[ 6], gds[ 7])
#define  GDS_TRI_NI3         GET_UINT2(gds[ 8], gds[ 9])
#define  GDS_TRI_ND          GET_UINT3(gds[10], gds[11], gds[12])
#define  GDS_TRI_NI          GET_UINT3(gds[13], gds[14], gds[15])
#define  GDS_TRI_AFlag       GET_UINT1(gds[16])
#define  GDS_TRI_LatPP       GET_INT3(gds[17], gds[18], gds[19])
#define  GDS_TRI_LonPP       GET_INT3(gds[20], gds[21], gds[22])
#define  GDS_TRI_LonMPL      GET_INT3(gds[23], gds[24], gds[25])
#define  GDS_TRI_BFlag       GET_UINT1(gds[27])

/* Spectral */
#define  GDS_PentaJ          GET_UINT2(gds[ 6], gds[ 7])
#define  GDS_PentaK          GET_UINT2(gds[ 8], gds[ 9])
#define  GDS_PentaM          GET_UINT2(gds[10], gds[11])
#define  GDS_RepType         GET_UINT1(gds[12])
#define  GDS_RepMode         GET_UINT1(gds[13])

/* Regular grid */
#define  GDS_NumLon          GET_UINT2(gds[ 6], gds[ 7])
#define  GDS_NumLat          GET_UINT2(gds[ 8], gds[ 9])
#define  GDS_FirstLat        GET_INT3(gds[10], gds[11], gds[12])
#define  GDS_FirstLon        GET_INT3(gds[13], gds[14], gds[15])
#define  GDS_ResFlag         GET_UINT1(gds[16])
#define  GDS_LastLat         GET_INT3(gds[17], gds[18], gds[19])
#define  GDS_LastLon         GET_INT3(gds[20], gds[21], gds[22])
#define  GDS_LonIncr         GET_INT2(gds[23], gds[24])
#define  GDS_LatIncr         GET_INT2(gds[25], gds[26])
#define  GDS_NumPar          GET_INT2(gds[25], gds[26])
#define  GDS_ScanFlag        GET_UINT1(gds[27])
#define  GDS_LatSP           GET_INT3(gds[32], gds[33], gds[34])
#define  GDS_LonSP           GET_INT3(gds[35], gds[36], gds[37])
#define  GDS_RotAngle        GET_Real(&(gds[38]))

/* Section 3: Bit Map Section (BMS) */

#define  BMS_Len	     ((bms) == NULL ? 0 : (int) (bms[0]<<16)+(bms[1]<<8)+bms[2])
#define  BMS_UnusedBits      (bms[3])
#define  BMS_Numeric         
#define  BMS_Bitmap	     ((bms) == NULL ? NULL : (bms)+6)
#define  BMS_BitmapSize      (((((bms[0]<<16)+(bms[1]<<8)+bms[2]) - 6)<<3) - bms[3])

/* Section 4: Binary Data Section (BDS) */

#define  BDS_Len	    ((int) ((bds[0]<<16)+(bds[1]<<8)+bds[2]))
#define  BDS_Flag	    (bds[3])
#define  BDS_BinScale       GET_INT2(bds[ 4], bds[ 5])
#define  BDS_RefValue       decfp2((int)bds[ 6], GET_UINT3(bds[ 7], bds[ 8], bds[ 9]))
#define  BDS_NumBits        ((int) bds[10])
#define  BDS_RealCoef       decfp2((int)bds[11], GET_UINT3(bds[12], bds[13], bds[14]))
#define  BDS_Power          GET_INT2(bds[13], bds[14])

/* Section 5: End Section (ES) */

#endif  /* _GRIBDECODE_H */
#ifndef _FILE_H
#define _FILE_H

#include <stdio.h>
#include <sys/types.h>


#define  FILE_UNDEFID   -1

#define  FILE_TYPE_STD   1
#define  FILE_TYPE_MMAP  2

const
char  *fileLibraryVersion(void);

void   fileDebug(int debug);

int    fileSetBufferType(int fileID, int type);
void   fileSetBufferSize(int fileID, long buffersize);

int    fileOpen(const char *filename, const char *mode);
void   fileClose(int fileID);

char  *fileInqName(int fileID);
int    fileInqMode(int fileID);

int    fileFlush(int fileID);
void   fileClearerr(int fileID);
int    fileEOF(int fileID);
void   fileRewind(int fileID);

off_t  fileGetPos(int fileID);
int    fileSetPos(int fileID, off_t offset, int whence);

int    fileGetc(int fileID);

size_t fileRead(int fileID, void *ptr, size_t size);
size_t fileWrite(int fileID, const void *ptr, size_t size);

#endif  /* _FILE_H */
#ifndef _PBIO_H
#define _PBIO_H

void  pbOpen   (int *unit, const char *name, const char *mode, int *iret);
void  pbClose  (int  unit, int *iret);
void  pbRead   (int unit, void *buffer, int nbytes, int *iret);
void  pbWrite  (int unit, void *buffer, int nbytes, int *iret);
void  pbSeek   (int unit, int offset, int whence, int *iret);
void  pbTell   (int unit, int *iret);
void  pbFlush  (int unit);
void  pbSize   (int unit, int *plen);
void  pbGrib   (int unit, void *buffer, int inlen, int *outlen, int *iret);
void  pbSetRaw (int unit, int raw);

#endif  /* _PBIO_H */
#ifndef _GRIBFORTRAN_H
#define _GRIBFORTRAN_H

#if ! defined   (_PBIO_H)
#  include "pbio.h"
#endif
#if ! defined   (_GRIB_H)
#  include "grib.h"
#endif

void gribEx(int *isec0, int *isec1, int *isec2, void *fsec2, int *isec3,
	    void *fsec3, int *isec4, void *fsec4, int klenp, int *kgrib,
	    int kleng, int *kword, char *hoper, int *kret);

void gribprs0  (int *isec0);
void gribprs1  (int *isec0, int *isec1);
void gribprs2  (int *isec0, int *isec2, void   *fsec2);
void gribprs2dp(int *isec0, int *isec2, double *fsec2);
void gribprs2sp(int *isec0, int *isec2, float  *fsec2);
void gribprs3  (int *isec0, int *isec3, void   *fsec3);
void gribprs3dp(int *isec0, int *isec3, double *fsec3);
void gribprs3sp(int *isec0, int *isec3, float  *fsec3);
void gribprs4  (int *isec0, int *isec4, void   *fsec4);
void gribprs4dp(int *isec0, int *isec4, double *fsec4);
void gribprs4sp(int *isec0, int *isec4, float  *fsec4);
void gribprs4w (int *isec4);

void gribsdbg(int debug);
void gribsrnd(int round);
void gribsref(void *fref);
void gribsrefdp(double fref);
void gribsrefsp(float  fref);
void gribsvck(int vcheck);

#endif  /* _GRIBFORTRAN_H */

#define  NINT(x)  ((x) < 0 ? (int)((x)-.5) : (int)((x)+.5))

static int confp3_init = 0;
static double rpow16m70tab[128];
static double rlog16;


static void
init_confp3()
{
  int iexp;

  rlog16 = 1.0 / log(16.0);

  for ( iexp = 0; iexp < 128; iexp++ )
    rpow16m70tab[iexp] = 1.0 / pow(16.0, (double)(iexp - 70));
  
  confp3_init = 1;
}

void
confp3(double pval, int *kexp, int *kmant, int kbits, int kround)
{
  /*

    Purpose:
    --------

    Convert floating point number from machine
    representation to GRIB representation.

    Input Parameters:
    -----------------

       pval    - Floating point number to be converted.
       kbits   - Number of bits in computer word.
       kround  - Conversion type.
                 0 , Closest number in GRIB format less than
                     original number.
                 1 , Closest number in GRIB format to the
                     original number (equal to, greater than or
                     less than original number).

    Output Parameters:
    ------------------

       kexp    - 8 Bit signed exponent.
       kmant   - 24 Bit mantissa.

    Method:
    -------

    Floating point number represented as 8 bit signed
    exponent and 24 bit mantissa in integer values.

    Externals.
    ----------

    decfp2    - Decode from IBM floating point format.

    Reference:
    ----------

    WMO Manual on Codes re GRIB representation.

    Comments:
    ---------

    Routine aborts if an invalid conversion type parameter
    is used or if a 24 bit mantissa is not produced.

    Author:
    -------
     
    John Hennessy   ECMWF   18.06.91

    Modifications:
    --------------

    Uwe Schulzweida   MPIfM   01/04/2001

    Convert to C from EMOS library version 130

    Uwe Schulzweida   MPIfM   02/08/2002

     - speed up by factor 1.6 on NEC SX6
        - replace 1.0 / pow(16.0, (double)(iexp - 70)) by rpow16m70tab[iexp]

    Uwe Schulzweida   MPIfM   25/01/2003

     - more speed up
        - replace 1.0 / pow(16.0, (double)(iexp - 70)) by rpow16m70tab[iexp]
	- replace 1.0 / log(16.0) by rlog16
  */

  static char func[] = "confp3";
  double zval, rpowref;
  double zref, zeps;
  int iexp, isign;
  int iround;
  extern int GRB_Debug;

  /* ----------------------------------------------------------------- */
  /*   Section 1 . Initialise                                          */
  /* ----------------------------------------------------------------- */

  /*  Check conversion type parameter. */

  iround = kround;
  if ( iround != 0 && iround != 1 )
    {
      Error(func, "Invalid conversion type = %d", iround);

      /*  If not aborting, arbitrarily set rounding to 'up'. */
     iround = 1;
    }

  /* ----------------------------------------------------------------- */
  /*   Section 2 . Convert value of zero.                              */
  /* ----------------------------------------------------------------- */

  if ( ! fabs(pval) )
    {
      *kexp  = 0;
      *kmant = 0;
      iexp   = 0;
      isign  = 0;
      goto LABEL900;
    }

  /* ----------------------------------------------------------------- */
  /*   Section 3 . Convert other values.                               */
  /* ----------------------------------------------------------------- */

  if ( ! confp3_init ) init_confp3();

  zeps = 1.0e-12;
  if ( kbits == 32 ) zeps = 1.0e-8;
  zref = pval;

  /*  Sign of value. */

  isign = 0;
  if ( zref < 0.0 )
    {
      isign = 128;
      zref  = - zref;
    }

  /*  Exponent. */

  iexp = (int) (log(zref) * rlog16 + 64.0 + 1.0 + zeps);

  if ( iexp < 0   ) iexp = 0;
  if ( iexp > 127 ) iexp = 127;

  /*
  rpowref = zref / pow(16.0, (double)(iexp - 70));
  */

  rpowref = zref * rpow16m70tab[iexp];

  /*  Mantissa. */

  if ( iround == 0 )
    {
      /*  Closest number in GRIB format less than original number. */
      /*  Truncate for positive numbers. */
      /*  Round up for negative numbers. */

      if ( isign == 0 )
	*kmant = (int) rpowref;
      else
	*kmant = NINT(rpowref + 0.5);
    }
  else
    {
      /*  Closest number in GRIB format to the original number   */
      /*  (equal to, greater than or less than original number). */

      *kmant = NINT(rpowref);
    }

  /*  Check that mantissa value does not exceed 24 bits. */
  /*  If it does, adjust the exponent upwards and recalculate */
  /*  the mantissa. */
  /*  16777215 = 2**24 - 1 */

  if ( *kmant > 16777215 )
    {

    LABEL350:

      ++iexp;

      /*  Check for exponent overflow during adjustment  */

      if ( iexp > 127 )
	{
          Message(func, "Exponent overflow");
          Message(func, "Original number = %30.20f", pval);
          Message(func, "Sign = %3d, Exponent = %3d, Mantissa = %12d",
		  isign, iexp, *kmant);

	  Error(func, "Exponent overflow");

	  /*  If not aborting, arbitrarily set value to zero  */

          Message(func, "Value arbitrarily set to zero.");
          *kexp  = 0;
          *kmant = 0;
          iexp  = 0;
          isign = 0;
          goto LABEL900;
	}

      rpowref = zref * rpow16m70tab[iexp];

      if ( iround == 0 )
	{
	  /*  Closest number in GRIB format less than original number. */
	  /*  Truncate for positive numbers. */
	  /*  Round up for negative numbers. */

	  if ( isign == 0 )
	    *kmant = (int) rpowref;
	  else
	    *kmant = NINT(rpowref + 0.5);
	}
      else
	{
	  /*  Closest number in GRIB format to the original number */
	  /*  (equal to, greater or less than original number). */

	  *kmant = NINT(rpowref);
	}

      /*  Repeat calculation (with modified exponent) if still have */
      /*  mantissa overflow. */

      if ( *kmant > 16777215 ) goto LABEL350;
    }

  /*  Add sign bit to exponent. */

  *kexp = iexp + isign;

  /* ----------------------------------------------------------------- */
  /*   Section 9. Return                                               */
  /* ----------------------------------------------------------------- */

LABEL900:

  if ( GRB_Debug )
    {
      Message(func, "Conversion type parameter = %4d", kround);
      Message(func, "Original number = %30.20f", pval);

      zval = decfp2(*kexp, *kmant);

      Message(func, "Converted to      %30.20f", zval);
      Message(func, "Sign = %3d, Exponent = %3d, Mantissa = %12d",
	      isign, iexp, *kmant);
    }

  return;
} /* confp3 */

static double pow16m64tab[128];
static int decfp2_init = 0;

static void init_decfp2(void)
{
  int iexp;

  for ( iexp = 0; iexp < 128; iexp++ )
    pow16m64tab[iexp] = pow(16.0, (double)(iexp - 64));
  
  decfp2_init = 1;
}

double decfp2(int kexp, int kmant)
{
  /*

    Purpose:
    --------

    Convert GRIB representation of a floating point
    number to machine representation.

    Input Parameters:
    -----------------

    kexp    - 8 Bit signed exponent.
    kmant   - 24 Bit mantissa.

    Output Parameters:
    ------------------

    Return value   - Floating point number represented
                     by kexp and kmant.

    Method:
    -------

    Floating point number represented as 8 bit exponent
    and 24 bit mantissa in integer values converted to
    machine floating point format.

    Externals:
    ----------

    None.

    Reference:
    ----------

    WMO Manual on Codes re GRIB representation.

    Comments:
    ---------

    Rewritten from DECFP, to conform to programming standards.
    Sign bit on 0 value now ignored, if present.
    If using 32 bit reals, check power of 16 is not so small as to
    cause overflows (underflows!); this causes warning to be given
    on Fujitsus.

    Author:
    -------

    John Hennessy   ECMWF   18.06.91

    Modifications:
    --------------

    Uwe Schulzweida   MPIfM   01/04/2001

     - Convert to C from EMOS library version 130
     
    Uwe Schulzweida   MPIfM   02/08/2002

     - speed up by factor 2 on NEC SX6
        - replace pow(2.0, -24.0) by constant POW_2_M24
        - replace pow(16.0, (double)(iexp - 64)) by pow16m64tab[iexp]
  */

  static char func[] = "decfp2";
  double pval;
  int iexp, isign;
  extern int GRB_Debug;

  /* ----------------------------------------------------------------- */
  /*   Section 1 . Convert value of 0.0. Ignore sign bit.              */
  /* ----------------------------------------------------------------- */

  if ( GRB_Debug ) Message(func, "KEXP = %d  KMANT = %d", kexp, kmant);
  /*
  if ( (kexp == 128 || kexp == 0) && kmant == 0 )
  */
  if ( (kexp == 128) || (kexp == 0) || (kexp == 255) )
    {
      pval = 0.0;
      goto LABEL900;
    }

  /* ----------------------------------------------------------------- */
  /*   Section 2 . Convert other values.                               */
  /* ----------------------------------------------------------------- */

  /*  Sign of value. */

  iexp = kexp;
  isign = 1;

  if ( iexp >= 128 )
    {
      iexp -= 128;
      isign = -1;
    }

  /*  Decode value. */

  /*
  pval = isign * pow(2.0, -24.0) * kmant * pow(16.0, (double)(iexp - 64));
  */
#define POW_2_M24  0.000000059604644775390625  /*  pow(2.0, -24.0) */

  if ( ! decfp2_init ) init_decfp2();

  pval = isign * POW_2_M24 * kmant * pow16m64tab[iexp];

  /* ----------------------------------------------------------------- */
  /*   Section 9. Return to calling routine.                           */
  /* ----------------------------------------------------------------- */

LABEL900:

  if ( GRB_Debug ) Message(func, "Returned value = %f", pval);

  return (pval);
} /* decfp2 */




int gribRefDate(int *isec1)
{
  int date, ryear, rmonth, rday;

  ryear = ISEC1_Year;
  if ( ryear != 255 )
    ryear = ((ISEC1_Century-1)*100 + ISEC1_Year);
  else
    ryear = 1;

  rmonth  = ISEC1_Month;
  rday    = ISEC1_Day;

  date = ryear*10000 + rmonth*100 + rday;

  return (date) ;
}

int gribRefTime(int *isec1)
{
  int time, rhour, rminute;

  rhour   = ISEC1_Hour;
  rminute = ISEC1_Minute;

  time = rhour*100 + rminute;

  return (time) ;
}

int gribTimeIsFC(int *isec1)
{
  int isFC = FALSE;

  if ( ISEC1_TimePeriod1 > 0 )
    {
      if ( ISEC1_TimeRange == 0 || ISEC1_TimeRange == 10 ) isFC = TRUE;
    }

  return (isFC);
}

void gribDateTime(int *isec1, int *date, int *time)
{
  static char func[] = "gribDateTime";
  static int lprint = TRUE;
  int ryear, rmonth, rday, rhour, rminute;
  double second = 0.0;
  double rvalue;
  double add;
  sUnit unit = {0, 1};

  ryear = ISEC1_Year;
  if ( ryear != 255 )
    ryear = ((ISEC1_Century-1)*100 + ISEC1_Year);
  else
    ryear = 1;

  rmonth  = ISEC1_Month;
  rday    = ISEC1_Day;

  rhour   = ISEC1_Hour;
  rminute = ISEC1_Minute;
  /*
  printf("ref %d/%d/%d %d:%d\n", ryear, rmonth, rday, rhour, rminute);
  */
  if ( ISEC1_TimePeriod1 > 0 )
    {
      if ( ISEC1_TimeRange == 0 || ISEC1_TimeRange == 10 )
	{
	  InvCalendar(ryear, rmonth, rday, rhour, rminute, second, unit, &rvalue);

	  add = 0;
	  switch ( ISEC1_TimeUnit )
	    {
	    case ISEC1_TABLE4_MINUTE:  add =    60 * ISEC1_TimePeriod1; break;
	    case ISEC1_TABLE4_HOUR:    add =  3600 * ISEC1_TimePeriod1; break;
	    case ISEC1_TABLE4_DAY:     add = 86400 * ISEC1_TimePeriod1; break;
	    default:
	      if ( lprint )
		{
		  gprintf(func, "Time unit %d unsupported", ISEC1_TimeUnit);
		  lprint = FALSE;
		}
	    }
	  rvalue += add;

	  Calendar(rvalue, unit, &ryear, &rmonth, &rday, &rhour, &rminute, &second);
	}
    }
  /*
  printf("new %d/%d/%d %d:%d\n", ryear, rmonth, rday, rhour, rminute);
  */
  *date = ryear*10000 + rmonth*100 + rday;
  *time = rhour*100 + rminute;

  return;
}

void gprintf(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

   fprintf(grprsm, "%-18s : ", caller);
  vfprintf(grprsm, fmt, args);
   fprintf(grprsm, "\n");

  va_end(args);
}


void
gribExDP(int *isec0, int *isec1, int *isec2, double *fsec2, int *isec3,
	 double *fsec3, int *isec4, double *fsec4, int klenp, int *kgrib,
	 int kleng, int *kword, char *hoper, int *kret)
{
  static char func[] = "gribExDP";
  int yfunc = *hoper;

  if ( yfunc == 'D' || yfunc == 'J' || yfunc == 'R' )
    gribDecode(isec0, isec1, isec2, fsec2, isec3,
	       fsec3, isec4, fsec4, klenp, kgrib,
	       kleng, kword, yfunc, kret);
  else if ( *hoper == 'C' )
    gribEncode(isec0, isec1, isec2, fsec2, isec3,
	       fsec3, isec4, fsec4, klenp, kgrib,
	       kleng, kword, yfunc, kret);
  else
    {
      Error(func, "oper %c unsupported\n", yfunc);
      *kret=-9;
    }
}

void
gribExSP(int *isec0, int *isec1, int *isec2, float *fsec2sp, int *isec3,
	 float *fsec3sp, int *isec4, float *fsec4sp, int klenp, int *kgrib,
	 int kleng, int *kword, char *hoper, int *kret)
{
  static char func[] = "gribExSP";
  int inum, j;
  double fsec2dp[1024];
  double fsec3dp[2];
  double *fsec4dp = NULL;
  int yfunc = *hoper;

  if ( yfunc == 'C' )
    {
      inum = 10 + isec2[11];
      for ( j = 0; j < inum; j++ ) fsec2dp[j] = fsec2sp[j];

      fsec3dp[0] = fsec3sp[0];
      fsec3dp[1] = fsec3sp[1];

      inum = isec4[0];
      fsec4dp = (double *) malloc(inum*sizeof(double));
      if ( fsec4dp == NULL ) SysError(func, "No Memory!");

      for ( j = 0; j < inum; j++ ) fsec4dp[j] = fsec4sp[j];

      gribExDP(isec0, isec1, isec2, fsec2dp, isec3,
	       fsec3dp, isec4, fsec4dp, klenp, kgrib,
	       kleng, kword, hoper, kret);

      free(fsec4dp);
    }
  else if ( yfunc == 'D' || yfunc == 'J' || yfunc == 'R' )
    {
      fsec4dp = (double *) malloc(klenp*sizeof(double));
      if ( fsec4dp == NULL ) SysError(func, "No Memory!");

      gribExDP(isec0, isec1, isec2, fsec2dp, isec3,
	       fsec3dp, isec4, fsec4dp, klenp, kgrib,
	       kleng, kword, hoper, kret);

      inum = 10 + isec2[11];
      for ( j = 0; j < inum; j++ ) fsec2sp[j] = fsec2dp[j];

      fsec3sp[0] = fsec3dp[0];
      fsec3sp[1] = fsec3dp[1];

      inum = isec4[0];
      for ( j = 0; j < inum; j++ ) fsec4sp[j] = fsec4dp[j];

      free(fsec4dp);
    }
  else
    {
      Error(func, "oper %c unsupported\n", yfunc);
      *kret=-9;
    }
}

int GRB_Debug = 0;    /* If set to 1, debugging */

void
gribSetDebug(int debug)
{
  static char func[] = "gribDebug";

  GRB_Debug = debug;

  if ( GRB_Debug )
    Message(func, "debug level %d", debug);
}

void
gribSetRound(int round)
{
}

void
gribSetRefDP(double refval)
{
}

void
gribSetRefSP(float refval)
{
  gribSetRefDP((double) refval);
}

void
gribSetValueCheck(int vcheck)
{
}



void gribPrintSec0(int *isec0)
{
  /*

    Print the information in the Indicator
    Section (Section 0) of decoded GRIB data.

    Input Parameters:

       isec0 - Array of decoded integers from Section 0


    Converted from EMOS routine GRPRS0.

       Uwe Schulzweida   MPIfM   01/04/2001

  */

  grsdef();

  fprintf(grprsm, " \n");
  fprintf(grprsm, " Section 0 - Indicator Section.       \n");
  fprintf(grprsm, " -------------------------------------\n");
  fprintf(grprsm, " Length of GRIB message (octets).     %9d\n", ISEC0_GRIB_Len);
  fprintf(grprsm, " GRIB Edition Number.                 %9d\n", ISEC0_GRIB_Version);
}

void gribPrintSec1(int *isec0, int *isec1)
{
  /*

    Print the information in the Product Definition
    Section (Section 1) of decoded GRIB data.

    Input Parameters:

       isec0 - Array of decoded integers from Section 0

       isec1 - Array of decoded integers from Section 1

    Comments:

       When decoding data from Experimental Edition or Edition 0,
       routine GRIBEX adds the additional fields available in
       Edition 1.


    Converted from EMOS routine GRPRS1.

       Uwe Schulzweida   MPIfM   01/04/2001

  */

  int iprev, icurr, icount, ioffset;
  int ibit, ierr, iout, iyear;
  int jloop, jiloop;
  float value;

  char hversion[9];
  /*
  char hfirst[121], hsecond[121], hthird[121], hfourth[121];
  */

  grsdef();

  /*
    -----------------------------------------------------------------
    Section 0 . Print required information.
    -----------------------------------------------------------------
  */

  fprintf(grprsm, " \n");
  fprintf(grprsm, " Section 1 - Product Definition Section.\n");
  fprintf(grprsm, " ---------------------------------------\n");

  fprintf(grprsm, " Code Table 2 Version Number.         %9d\n", isec1[0]);
  fprintf(grprsm, " Originating centre identifier.       %9d\n", isec1[1]);
  fprintf(grprsm, " Model identification.                %9d\n", isec1[2]);
  fprintf(grprsm, " Grid definition.                     %9d\n", isec1[3]);

  ibit = 8;
  prtbin(isec1[4], ibit, &iout, &ierr);
  fprintf(grprsm, " Flag (Code Table 1)                   %8.8d\n", iout);
  fprintf(grprsm, " Parameter identifier (Code Table 2). %9d\n", isec1[5]);

  /*
      IERR = CHKTAB2(ISEC1,HFIRST,HSECOND,HTHIRD,HFOURTH)
      IF( IERR .EQ. 0 ) THEN
       DO JLOOP = 121, 1, -1
          IF( HSECOND(JLOOP:JLOOP).NE.' ' ) THEN
            IOFFSET = JLOOP
            GOTO 110
          ENDIF
        ENDDO
        GOTO 120
 110    CONTINUE
        WRITE(*,'(2H ",A,1H")') HSECOND(1:IOFFSET)
 120    CONTINUE
      ENDIF
  */

  if ( isec1[5] != 127 )
    {
      fprintf(grprsm, " Type of level (Code Table 3).        %9d\n", isec1[6]);
      fprintf(grprsm, " Value 1 of level (Code Table 3).     %9d\n", isec1[7]);
      fprintf(grprsm, " Value 2 of level (Code Table 3).     %9d\n", isec1[8]);
    }
  else
    {
      fprintf(grprsm, " Satellite identifier.                %9d\n", isec1[6]);
      fprintf(grprsm, " Spectral band.                       %9d\n", isec1[7]);
    }

  iyear = isec1[9];
  if ( iyear != 255 )
    {
      iyear  = ((isec1[20]-1)*100 + isec1[9]);
      fprintf(grprsm, " Year of reference time of data.      %9d  (%4d)\n", isec1[9], iyear);
    }
  else
    {
      fprintf(grprsm, " Year of reference time of data MISSING  (=255)\n");
    }

  fprintf(grprsm, " Month of reference time of data.     %9d\n", isec1[10]);
  fprintf(grprsm, " Day of reference time of data.       %9d\n", isec1[11]);
  fprintf(grprsm, " Hour of reference time of data.      %9d\n", isec1[12]);
  fprintf(grprsm, " Minute of reference time of data.    %9d\n", isec1[13]);
  fprintf(grprsm, " Time unit (Code Table 4).            %9d\n", isec1[14]);
  fprintf(grprsm, " Time range one.                      %9d\n", isec1[15]);
  fprintf(grprsm, " Time range two.                      %9d\n", isec1[16]);
  fprintf(grprsm, " Time range indicator (Code Table 5)  %9d\n", isec1[17]);
  fprintf(grprsm, " Number averaged.                     %9d\n", isec1[18]);
  fprintf(grprsm, " Number missing from average.         %9d\n", isec1[19]);
  /*
     All ECMWF data in GRIB Editions before Edition 1 is decoded
     as 20th century data. Other centres are decoded as missing.
  */
  if ( isec0[1] < 1 && isec1[1] != 98 )
    fprintf(grprsm, " Century of reference time of data.   Not given\n");
  else
    fprintf(grprsm, " Century of reference time of data.   %9d\n", isec1[20]);

  /*   Print sub-centre  */
  fprintf(grprsm, " Sub-centre identifier.               %9d\n", isec1[21]);

  /*   Decimal scale factor  */
  fprintf(grprsm, " Units decimal scaling factor.        %9d\n", isec1[22]);

  /*
    -----------------------------------------------------------------
    Section 1 . Print local DWD information.
    -----------------------------------------------------------------
  */
  if ( (isec1[ 1] == 78 && isec1[36] == 253) ||
       (isec1[ 1] == 78 && isec1[36] == 254) )
    {
      fprintf(grprsm, " DWD local usage identifier.          %9d\n", isec1[36]);
      if ( isec1[36] == 253 )
	fprintf(grprsm, " (Database labelling and ensemble forecast)\n");
      if ( isec1[36] == 254 )
	fprintf(grprsm, " (Database labelling)\n");

      fprintf(grprsm, " Year of database entry                     %3d  (%4d)\n", isec1[43], 1900+isec1[43]);
      fprintf(grprsm, " Month of database entry                    %3d\n", isec1[44]);
      fprintf(grprsm, " Day of database entry                      %3d\n", isec1[45]);
      fprintf(grprsm, " Hour of database entry                     %3d\n", isec1[46]);
      fprintf(grprsm, " Minute of database entry                   %3d\n", isec1[47]);
      fprintf(grprsm, " DWD experiment number                %9d\n",isec1[48]);
      fprintf(grprsm, " DWD run type                         %9d\n",isec1[49]);
      if ( isec1[36] == 253 ) 
	{
	  fprintf(grprsm, " User id                              %9d\n",isec1[50]);
	  fprintf(grprsm, " Experiment identifier                %9d\n",isec1[51]);
	  fprintf(grprsm, " Ensemble identification type         %9d\n",isec1[52]);
	  fprintf(grprsm, " Number of ensemble members           %9d\n",isec1[53]);
	  fprintf(grprsm, " Actual number of ensemble member     %9d\n",isec1[54]);
	  fprintf(grprsm, " Model version                            %2d.%2.2d\n",isec1[55],isec1[56]);
	}
    }

  /*
    -----------------------------------------------------------------
    Section 2 . Print local ECMWF information.
    -----------------------------------------------------------------
  */
  /*
    Regular MARS labelling, or reformatted Washington EPS products.
  */
  if ( (isec1[ 1] == 98 && isec1[23] ==  1) ||
       (isec1[21] == 98 && isec1[23] ==  1) ||
       (isec1[ 1] ==  7 && isec1[21] == 98) )
    {
      /*   Parameters common to all definitions.  */

      fprintf(grprsm, " ECMWF local usage identifier.        %9d\n", isec1[36]);
      if ( isec1[36] == 1 )
	fprintf(grprsm, " (Mars labelling or ensemble forecast)\n");
      if ( isec1[36] == 2 )
        fprintf(grprsm, " (Cluster means and standard deviations)\n");
      if ( isec1[36] == 3 )
        fprintf(grprsm, " (Satellite image data)\n");
      if ( isec1[36] == 4 )
        fprintf(grprsm, " (Ocean model data)\n");
      if ( isec1[36] == 5 )
        fprintf(grprsm, " (Forecast probability data)\n");
      if ( isec1[36] == 6 )
        fprintf(grprsm, " (Surface temperature data)\n");
      if ( isec1[36] == 7 )
        fprintf(grprsm, " (Sensitivity data)\n");
      if ( isec1[36] == 8 )
        fprintf(grprsm, " (ECMWF re-analysis data)\n");
      if ( isec1[36] == 9 )
        fprintf(grprsm, " (Singular vectors and ensemble perturbations)\n");
      if ( isec1[36] == 10 )
        fprintf(grprsm, " (EPS tubes)\n");
      if ( isec1[36] == 11 )
        fprintf(grprsm, " (Supplementary data used by analysis)\n");
      if ( isec1[36] == 13 )
        fprintf(grprsm, " (Wave 2D spectra direction and frequency)\n");

      fprintf(grprsm, " Class.                               %9d\n", isec1[37]);
      fprintf(grprsm, " Type.                                %9d\n", isec1[38]);
      fprintf(grprsm, " Stream.                              %9d\n", isec1[39]);

      sprintf(hversion, "%8d", isec1[40]);
      fprintf(grprsm, " Version number or Experiment identifier.  %4s\n", &hversion[4]);
      /*
	ECMWF Local definition 1.
	(MARS labelling or ensemble forecast data)
      */
      if ( isec1[36] == 1 )
	{
	  fprintf(grprsm, " Forecast number.                     %9d\n", isec1[41]);
	  if ( isec1[39] != 1090 )
	    fprintf(grprsm, " Total number of forecasts.           %9d\n", isec1[42]);

	  return;
	}
      /*
	ECMWF Local definition 2.
	(Cluster means and standard deviations)
      */
      if ( isec1[36] == 2 )
	{
	  fprintf(grprsm, " Cluster number.                      %9d\n", isec1[41]);
	  fprintf(grprsm, " Total number of clusters.            %9d\n", isec1[42]);
	  fprintf(grprsm, " Clustering method.                   %9d\n", isec1[43]);
	  fprintf(grprsm, " Start time step when clustering.     %9d\n", isec1[44]);
	  fprintf(grprsm, " End time step when clustering.       %9d\n", isec1[45]);
	  fprintf(grprsm, " Northern latitude of domain.         %9d\n", isec1[46]);
	  fprintf(grprsm, " Western longitude of domain.         %9d\n", isec1[47]);
	  fprintf(grprsm, " Southern latitude of domain.         %9d\n", isec1[48]);
	  fprintf(grprsm, " Eastern longitude of domain.         %9d\n", isec1[49]);
	  fprintf(grprsm, " Operational forecast in cluster      %9d\n", isec1[50]);
	  fprintf(grprsm, " Control forecast in cluster          %9d\n", isec1[51]);
	  fprintf(grprsm, " Number of forecasts in cluster.      %9d\n", isec1[52]);

	  for (jloop = 0; jloop < isec1[52]; jloop++)
	    fprintf(grprsm, " Forecast number                      %9d\n", isec1[jloop+53]);

	  return;
	}
      /*
	ECMWF Local definition 3.
	(Satellite image data)
      */
      if ( isec1[36] == 3 )
	{
	  fprintf(grprsm, " Satellite spectral band.             %9d\n", isec1[41]);
	  fprintf(grprsm, " Function code.                       %9d\n", isec1[42]);
	  return;
	}
      /*
	ECMWF Local definition 4.
	(Ocean model data)
      */
      if ( isec1[36] == 4 )
	{
	  fprintf(grprsm, " Satellite spectral band.             %9d\n", isec1[41]);
	  if ( isec1[39] != 1090 )
	    fprintf(grprsm, " Function code.                       %9d\n", isec1[42]);
	  fprintf(grprsm, " Coordinate structure definition.\n");
	  fprintf(grprsm, " Fundamental spatial reference system.%9d\n", isec1[43]);
	  fprintf(grprsm, " Fundamental time reference.          %9d\n", isec1[44]);
	  fprintf(grprsm, " Space unit flag.                     %9d\n", isec1[45]);
	  fprintf(grprsm, " Vertical coordinate definition.      %9d\n", isec1[46]);
	  fprintf(grprsm, " Horizontal coordinate definition.    %9d\n", isec1[47]);
	  fprintf(grprsm, " Time unit flag.                      %9d\n", isec1[48]);
	  fprintf(grprsm, " Time coordinate definition.          %9d\n", isec1[49]);
	  fprintf(grprsm, " Position definition.     \n");
	  fprintf(grprsm, " Mixed coordinate field flag.         %9d\n", isec1[50]);
	  fprintf(grprsm, " Coordinate 1 flag.                   %9d\n", isec1[51]);
	  fprintf(grprsm, " Averaging flag.                      %9d\n", isec1[52]);
	  fprintf(grprsm, " Position of level 1.                 %9d\n", isec1[53]);
	  fprintf(grprsm, " Position of level 2.                 %9d\n", isec1[54]);
	  fprintf(grprsm, " Coordinate 2 flag.                   %9d\n", isec1[55]);
	  fprintf(grprsm, " Averaging flag.                      %9d\n", isec1[56]);
	  fprintf(grprsm, " Position of level 1.                 %9d\n", isec1[57]);
	  fprintf(grprsm, " Position of level 2.                 %9d\n", isec1[58]);
	  fprintf(grprsm, " Grid Definition.\n");
	  fprintf(grprsm, " Coordinate 3 flag (x-axis)           %9d\n", isec1[59]);
	  fprintf(grprsm, " Coordinate 4 flag (y-axis)           %9d\n", isec1[60]);
	  fprintf(grprsm, " Coordinate 4 of first grid point.    %9d\n", isec1[61]);
	  fprintf(grprsm, " Coordinate 3 of first grid point.    %9d\n", isec1[62]);
	  fprintf(grprsm, " Coordinate 4 of last grid point.     %9d\n", isec1[63]);
	  fprintf(grprsm, " Coordinate 3 of last grid point.     %9d\n", isec1[64]);
	  fprintf(grprsm, " i - increment.                       %9d\n", isec1[65]);
	  fprintf(grprsm, " j - increment.                       %9d\n", isec1[66]);
	  fprintf(grprsm, " Flag for irregular grid coordinates. %9d\n", isec1[67]);
	  fprintf(grprsm, " Flag for normal or staggered grids.  %9d\n", isec1[68]);
	  fprintf(grprsm, " Further information.\n");
	  fprintf(grprsm, " Further information flag.            %9d\n", isec1[69]);
	  fprintf(grprsm, " Auxiliary information.\n");
	  fprintf(grprsm, " No. entries in horizontal coordinate %9d\n", isec1[70]);
	  fprintf(grprsm, " No. entries in mixed coordinate defn.%9d\n", isec1[71]);
	  fprintf(grprsm, " No. entries in grid coordinate list. %9d\n", isec1[72]);
	  fprintf(grprsm, " No. entries in auxiliary array.      %9d\n", isec1[73]);
	  /*
	    Horizontal coordinate supplement.
	  */
	  fprintf(grprsm, " Horizontal coordinate supplement.\n");
	  if ( isec1[70] == 0 )
	    {
	      fprintf(grprsm, "(None).\n");
	    }
	  else
	    {
	      fprintf(grprsm, "Number of items = %d\n", isec1[70]);
	      for (jloop = 0; jloop < isec1[70]; jloop++)
		fprintf(grprsm, "         %12d\n", isec1[74+jloop]);
	    }
	  /*
	    Mixed coordinate definition.
	  */
	  fprintf(grprsm, " Mixed coordinate definition.\n");
	  if ( isec1[71] == 0 )
	    {
	      fprintf(grprsm, "(None).\n");
	    }
	  else
	    {
	      fprintf(grprsm, "Number of items = %d\n", isec1[71]);
	      ioffset = 74 + isec1[70];
	      for (jloop = 0; jloop < isec1[71]; jloop++)
		fprintf(grprsm, "         %12d\n", isec1[ioffset+jloop]);
	    }
	  /*
	    Grid coordinate list.
	  */
	  fprintf(grprsm, " Grid coordinate list. \n");
	  if ( isec1[72] == 0 )
	    {
	      fprintf(grprsm, "(None).\n");
	    }
	  else
	    {
	      fprintf(grprsm, "Number of items = %d\n", isec1[72]);
	      ioffset = 74 + isec1[70] + isec1[71];
	      for (jloop = 0; jloop < isec1[72]; jloop++)
		fprintf(grprsm, "         %12d\n", isec1[ioffset+jloop]);
	    }
	  /*
	    Auxiliary array.
	  */
	  fprintf(grprsm, " Auxiliary array.      \n");
	  if ( isec1[73] == 0 )
	    {
	      fprintf(grprsm, "(None).\n");
	    }
	  else
	    {
	      fprintf(grprsm, "Number of items = %d\n", isec1[73]);
	      ioffset = 74 + isec1[70] + isec1[71] + isec1[72];
	      for (jloop = 0; jloop < isec1[73]; jloop++)
		fprintf(grprsm, "         %12d\n", isec1[ioffset+jloop]);
	    }
	  /*
	    Post-auxiliary array.
	  */
	  fprintf(grprsm, " Post-auxiliary array. \n");
	  ioffset = 74 + isec1[70] + isec1[71] + isec1[72] + isec1[73];
	  if ( isec1[ioffset] == 0 )
	    {
	      fprintf(grprsm, "(None).\n");
	    }
	  else
	    {
	      fprintf(grprsm, "Number of items = %d\n", isec1[ioffset]);
	      for (jloop = 1; jloop < isec1[ioffset]; jloop++)
		fprintf(grprsm, "         %12d\n", isec1[ioffset+jloop]);
	    }

	  return;
	}
      /*
	ECMWF Local definition 5.
	(Forecast probability data)
      */
      if ( isec1[36] == 5 )
	{
	  fprintf(grprsm, " Forecast probability number          %9d\n", isec1[41]);
	  fprintf(grprsm, " Total number of forecast probabilities %7d\n", isec1[42]);
	  fprintf(grprsm, " Threshold units decimal scale factor %9d\n", isec1[43]);
	  fprintf(grprsm, " Threshold indicator(1=lower,2=upper,3=both) %2d\n", isec1[44]);
	  if ( isec1[44]  !=  2 )
	    fprintf(grprsm, " Lower threshold value                %9d\n", isec1[45]);
	  if ( isec1[44]  !=  1 )
	    fprintf(grprsm, " Upper threshold value                %9d\n", isec1[46]);
	  return;
	}
      /*
	ECMWF Local definition 6.
	(Surface temperature data)
      */
      if ( isec1[36] == 6 )
	{
	  iyear = isec1[43];
	  if ( iyear > 100 )
	    {
	      if ( iyear < 19000000 ) iyear = iyear + 19000000;
	      fprintf(grprsm, " Date of SST field used               %9d\n", iyear);
	    }
	  else
	    fprintf(grprsm, "Date of SST field used               Not given\n");
	}
      if ( isec1[44] == 0 )
	fprintf(grprsm, " Type of SST field (= climatology)    %9d\n", isec1[44]);
      if ( isec1[44] == 1 )
	fprintf(grprsm, " Type of SST field (= 1/1 degree)     %9d\n", isec1[44]);
      if ( isec1[44] == 2 )
	fprintf(grprsm, " Type of SST field (= 2/2 degree)     %9d\n", isec1[44]);

      fprintf(grprsm, " Number of ICE fields used:           %9d\n", isec1[45]);

      for (jloop = 1; jloop <= isec1[45]; jloop++)
	{
	  iyear = isec1[44+(jloop*2)];
	  if ( iyear > 100 )
	    {
              if ( iyear < 19000000 ) iyear = iyear + 19000000;
	      fprintf(grprsm, " Date of ICE field%3d                 %9d\n", jloop, iyear);
	      fprintf(grprsm, " Satellite number (ICE field%3d)      %9d\n", jloop,
		     isec1[45+(jloop*2)]);
	    }
	  else
	    fprintf(grprsm, "Date of SST field used               Not given\n");
	}
      /*
	ECMWF Local definition 7.
	(Sensitivity data)
      */
      if ( isec1[36] == 7 )
	{
	  if ( isec1[38]  ==  51 )
	    fprintf(grprsm, " Forecast number                      %9d\n", isec1[41]);
	  if ( isec1[38]  !=  51 )
	    fprintf(grprsm, " Iteration number                     %9d\n", isec1[41]);
	  if ( isec1[38]  !=  52 )
	    fprintf(grprsm, " Total number of diagnostics          %9d\n", isec1[42]);
	  if ( isec1[38]  ==  52 )
	    fprintf(grprsm, " No.interations in diag. minimisation %9d\n", isec1[42]);
	  fprintf(grprsm, " Domain(0=Global,1=Europe,2=N.Hem.,3=S.Hem.) %2d\n", isec1[43]);
	  fprintf(grprsm, " Diagnostic number                    %9d\n", isec1[44]);
	}
      /*
	ECMWF Local definition 8.
	(ECMWF re-analysis data)
      */
      if ( isec1[36] == 8 )
	{
	  if ( (isec1[39] == 1043) ||
	       (isec1[39] == 1070) ||
	       (isec1[39] == 1071) )
	    {
	      fprintf(grprsm, " Interval between reference times     %9d\n", isec1[41]);
	      for (jloop = 43; jloop <= 54; jloop++)
		{
		  jiloop = jloop + 8;
		  fprintf(grprsm, " ERA section 1 octet %2d.              %9d\n",
			 jiloop, isec1[jloop-1]);
		}
	    }
	  else
	    {
	      for (jloop = 42; jloop <= 54; jloop++)
		{
		  jiloop = jloop + 8;
		  fprintf(grprsm, " ERA section 1 octet %2d.              %9d\n",
			 jiloop, isec1[jloop-1]);
		}
	    }
	  return;
	}

      if ( isec1[38] > 4  && isec1[38] < 9 )
	{
	  fprintf(grprsm, " Simulation number.                   %9d\n", isec1[41]);
	  fprintf(grprsm, " Total number of simulations.         %9d\n", isec1[42]);
	}
      /*
	ECMWF Local definition 9.
	(Singular vectors and ensemble perturbations)
      */
      if ( isec1[36] == 9 )
	{
	  if ( isec1[38] == 60 )
	    fprintf(grprsm, " Perturbed ensemble forecast number   %9d\n", isec1[41]);
	  if ( isec1[38] == 61 )
	    fprintf(grprsm, " Initial state perturbation number    %9d\n", isec1[41]);
	  if ( isec1[38] == 62 )
	    fprintf(grprsm, " Singular vector number               %9d\n", isec1[41]);
	  if ( isec1[38] == 62 )
	    {
	      fprintf(grprsm, " Number of iterations                 %9d\n", isec1[42]);
	      fprintf(grprsm, " Number of singular vectors computed  %9d\n", isec1[43]);
	      fprintf(grprsm, " Norm used at initial time            %9d\n", isec1[44]);
	      fprintf(grprsm, " Norm used at final time              %9d\n", isec1[45]);
	      fprintf(grprsm, " Multiplication factor                %9d\n", isec1[46]);
    	      fprintf(grprsm, " Latitude of north-west corner        %9d\n", isec1[47]);
    	      fprintf(grprsm, " Longitude of north-west corner       %9d\n", isec1[48]);
	      fprintf(grprsm, " Latitude of south-east corner        %9d\n", isec1[49]);
	      fprintf(grprsm, " Longitude of south-east corner       %9d\n", isec1[50]);
	      fprintf(grprsm, " Accuracy                             %9d\n", isec1[51]);
	      fprintf(grprsm, " Number of singular vectors evolved   %9d\n", isec1[52]);
	      fprintf(grprsm, " Ritz number one                      %9d\n", isec1[53]);
	      fprintf(grprsm, " Ritz number two                      %9d\n", isec1[54]);
	    }
	}
      /*
	ECMWF Local definition 10.
	(EPS tubes)
      */
      if ( isec1[36] == 10 )
	{
	  fprintf(grprsm, " Tube number                          %9d\n", isec1[41]);
          fprintf(grprsm, " Total number of tubes                %9d\n", isec1[42]);
          fprintf(grprsm, " Central cluster definition           %9d\n", isec1[43]);
          fprintf(grprsm, " Parameter                            %9d\n", isec1[44]);
          fprintf(grprsm, " Type of level                        %9d\n", isec1[45]);
          fprintf(grprsm, " Northern latitude of domain of tubing%9d\n", isec1[46]);
          fprintf(grprsm, " Western longitude of domain of tubing%9d\n", isec1[47]);
          fprintf(grprsm, " Southern latitude of domain of tubing%9d\n", isec1[48]);
          fprintf(grprsm, " Eastern longitude of domain of tubing%9d\n", isec1[49]);
          fprintf(grprsm, " Tube number of operational forecast  %9d\n", isec1[50]);
          fprintf(grprsm, " Tube number of control forecast      %9d\n", isec1[51]);
          fprintf(grprsm, " Height/pressure of level             %9d\n", isec1[52]);
          fprintf(grprsm, " Reference step                       %9d\n", isec1[53]);
          fprintf(grprsm, " Radius of central cluster            %9d\n", isec1[54]);
          fprintf(grprsm, " Ensemble standard deviation          %9d\n", isec1[55]);
          fprintf(grprsm, " Dist.of tube extreme to ensemble mean%9d\n", isec1[56]);
          fprintf(grprsm, " Number of forecasts in the tube      %9d\n", isec1[57]);

          fprintf(grprsm, " List of ensemble forecast numbers:\n");
          for (jloop = 1; jloop <=  isec1[57]; jloop++)
	    fprintf(grprsm, "    %9d\n", isec1[57+jloop]);
	}
      /*
	ECMWF Local definition 11.
	(Supplementary data used by the analysis)
      */
      if ( isec1[36] == 11 )
	{
	  fprintf(grprsm, " Details of analysis which used the supplementary data:\n");
	  fprintf(grprsm, "   Class                              %9d\n", isec1[41]);
	  fprintf(grprsm, "   Type                               %9d\n", isec1[42]);
	  fprintf(grprsm, "   Stream                             %9d\n", isec1[43]);
	  sprintf(hversion, "%8d", isec1[44]);
	  fprintf(grprsm, "   Version number/experiment identifier:   %4s\n", &hversion[4]);
	  iyear = isec1[45];
	  if ( iyear > 50 )
	    iyear = iyear + 1900;
	  else
	    iyear = iyear + 2000;

	  fprintf(grprsm, "   Year                               %9d\n", iyear);
	  fprintf(grprsm, "   Month                              %9d\n", isec1[46]);
	  fprintf(grprsm, "   Day                                %9d\n", isec1[47]);
	  fprintf(grprsm, "   Hour                               %9d\n", isec1[48]);
	  fprintf(grprsm, "   Minute                             %9d\n", isec1[49]);
	  fprintf(grprsm, "   Century                            %9d\n", isec1[50]);
	  fprintf(grprsm, "   Originating centre                 %9d\n", isec1[51]);
	  fprintf(grprsm, "   Sub-centre                         %9d\n", isec1[52]);
	}
      /*
	ECMWF Local definition 12.
      */
      if ( isec1[36] == 12 )
	{
	  fprintf(grprsm, " (Mean, average, etc)\n");
          fprintf(grprsm, " Start date of the period              %8d\n", isec1[41]);
          fprintf(grprsm, " Start time of the period                  %4.4d\n", isec1[42]);
          fprintf(grprsm, " Finish date of the period             %8d\n", isec1[43]);
          fprintf(grprsm, " Finish time of the period                 %4.4d\n", isec1[44]);
          fprintf(grprsm, " Verifying date of the period          %8d\n", isec1[45]);
          fprintf(grprsm, " Verifying time of the period              %4.4d\n", isec1[46]);
          fprintf(grprsm, " Code showing method                   %8d\n", isec1[47]);
          fprintf(grprsm, " Number of different time intervals used  %5d\n", isec1[48]);
          fprintf(grprsm, " List of different time intervals used:\n");
          iprev  = isec1[49];
          icurr  = 0;
          icount = 0;
          for (jloop = 1; jloop <= isec1[48]; jloop++)
	    {
	      icurr = isec1[48+jloop];
	      if ( icurr != iprev )
		{
		  if ( icount == 1 )
		    fprintf(grprsm, "  - interval %5.4d used       once\n", iprev);
		  if ( icount == 2 )
		    fprintf(grprsm, "  - interval %5.4d used       twice\n", iprev);
		  if ( icount > 2 )
		    fprintf(grprsm, "  - interval %5.4d used %5d times\n",  iprev, icount);
		  iprev  = icurr;
		  icount = 1;
		}
	      else
		icount = icount + 1;
	    }
	  if ( icount == 1 )
	    fprintf(grprsm, "  - interval %5.4d used       once\n", iprev);
	  if ( icount == 2 )
	    fprintf(grprsm, "  - interval %5.4d used       twice\n", iprev);
	  if ( icount > 2 )
	    fprintf(grprsm, "  - interval %5.4d used %5d times\n",  iprev, icount);
	}
      /*
	ECMWF Local definition 13.
	(Wave 2D spectra direction and frequency)
      */
      if ( isec1[36] == 13 )
	{
          fprintf(grprsm, " Direction number                     %9d\n", isec1[43]);
	  fprintf(grprsm, " Frequency number                     %9d\n", isec1[44]);
	  fprintf(grprsm, " Total number of directions           %9d\n", isec1[45]);
	  fprintf(grprsm, " Total number of frequencies          %9d\n", isec1[46]);
	  fprintf(grprsm, " Scale factor applied to directions   %9d\n", isec1[47]);
	  fprintf(grprsm, " Scale factor applied to frequencies  %9d\n", isec1[48]);
	  fprintf(grprsm, " List of directions:\n");
          for (jloop = 1; jloop <= isec1[45]; jloop++)
            {
	      value = (float)(isec1[48+jloop])/(float)(isec1[47]);
	      if ( isec1[43] == jloop )
		fprintf(grprsm, " %2.2d:%15.7f   <-- this field value\n",  jloop, value);
	      else
		fprintf(grprsm, "%2.2d:%15.7f\n",  jloop, value);
            }
	  fprintf(grprsm, " List of frequencies:\n");
          for (jloop = 1; jloop <= isec1[46]; jloop++)
	    {
	      value = (float)(isec1[48+isec1[45]+jloop])/(float)(isec1[48]);
	      if ( isec1[44] == jloop )
		fprintf(grprsm, " %2.2d:%15.7f   <-- this field value\n",  jloop, value);
	      else
		fprintf(grprsm, "%2.2d:%15.7f\n",  jloop, value);

	      if ( isec1[49+isec1[45]+isec1[46]] != 0 )
		{
		  fprintf(grprsm, " System number (65535 = missing)      %9d\n",
			 isec1[49+isec1[45]+isec1[46]]);
		  fprintf(grprsm, " Method number (65535 = missing)      %9d\n",
			 isec1[50+isec1[45]+isec1[46]]);
		}
	    }
	  /*
	    ECMWF Local definition 14.
	    (Brightness temperature)
	  */
	  if ( isec1[36] == 14 )
	    {
	      fprintf(grprsm, " Channel number                       %9d\n", isec1[43]);
	      fprintf(grprsm, " Scale factor applied to frequencies  %9d\n", isec1[44]);
	      fprintf(grprsm, " Total number of frequencies          %9d\n", isec1[45]);
	      fprintf(grprsm, " List of frequencies:\n");
              for (jloop = 1; jloop <= isec1[45]; jloop++)
		{
		  value = (float)(isec1[45+jloop])/(float)(isec1[44]);
		  if ( isec1[43] == jloop )
		    fprintf(grprsm, " %3d:%15.9f   <-- this channel\n", jloop, value);
		  else
		    fprintf(grprsm, " %3d:%15.9f\n", jloop, value);
		}
	    }
	  /*
	    ECMWF Local definition 15.
	    (Ocean ensemble seasonal forecast)
	  */
	  if ( isec1[36] == 15 )
	    {
	      fprintf(grprsm, " Ensemble member number               %9d\n", isec1[41]);
	      fprintf(grprsm, " System number                        %9d\n", isec1[42]);
	      fprintf(grprsm, " Method number                        %9d\n", isec1[43]);
	    }
	  /*
	    ECMWF Local definition 16.
	    (Seasonal forecast monthly mean atmosphere data)
	  */
        if ( isec1[36] == 16 )
	  {
	    fprintf(grprsm, " Ensemble member number               %9d\n", isec1[41]);
	    fprintf(grprsm, " System number                        %9d\n", isec1[43]);
	    fprintf(grprsm, " Method number                        %9d\n", isec1[44]);
	    fprintf(grprsm, " Verifying month                      %9d\n", isec1[45]);
	    fprintf(grprsm, " Averaging period                     %9d\n", isec1[46]);
	  }
	/*
	  ECMWF Local definition 17.
	  (Sst or sea-ice used by analysis)
	*/
        if ( isec1[36] == 17 )
	  {
	    iyear = isec1[43];
	    if ( iyear > 100 )
	      {
		if ( iyear < 19000000 ) iyear = iyear + 19000000;
		fprintf(grprsm, " Date of sst/ice field used           %9d\n", iyear);
	      }
	    else
              fprintf(grprsm, " Date of sst/ice field used           Not given\n");
      
	    if ( isec1[44] == 0 )
	      fprintf(grprsm, " Type of sst/ice field (= climatology)%9d\n", isec1[44]);
	    if ( isec1[44] == 1 )
	      fprintf(grprsm, " Type of sst/ice field (= 1/1 degree) %9d\n", isec1[44]);
	    if ( isec1[44] == 2 )
	      fprintf(grprsm, " Type of sst/ice field (= 2/2 degree) %9d\n", isec1[44]);

	    fprintf(grprsm, " Number of ICE fields used:           %9d\n", isec1[45]);

	    for (jloop = 1; jloop < isec1[45]; jloop++)
	      {
		iyear = isec1[44+(jloop*2)];
		if ( iyear > 100 )
		  {
		    if ( iyear < 19000000 ) iyear = iyear + 19000000;
		    fprintf(grprsm, " Date of ICE field%3d                 %9d\n", jloop,
			   iyear);
		    fprintf(grprsm, " Satellite number (ICE field%3d)      %9d\n", jloop,
			   isec1[45+(jloop*2)]);
		  }
		else
		  fprintf(grprsm, "Date of sst/ice field used           Not given\n");
	      } 
	  }
	}
    }
  /*
    -----------------------------------------------------------------
    Section 3 . Print Washington ensemble product information.
    -----------------------------------------------------------------
  */
  /*
    Washington EPS products (but not reformatted Washington EPS
    products.
  */
  if ( (isec1[1] == 7 && isec1[23] == 1) && (! isec1[21] == 98) )
    {
      /*   CALL KWPRS1 (iSEC0,iSEC1)*/
    }
}

void printQuasi(int *isec2)
{
  /*

    Print the qusai-regular information in the Grid Description
    Section (Section 2) of decoded GRIB data.

    Input Parameters:

       isec2 - Array of decoded integers from Section 2.

    Comments:

       Only data representation types catered for are Gaussian
       grid, latitude/longitude grid, Spherical Harmonics,
       Polar stereographic and Space view perspective.

    Converted from EMOS routine PTQUASI.

       Uwe Schulzweida   MPIfM   01/04/2001

  */

  char yout[64];
  int nextlat, nrepeat, latcnt;
  int j;
  int ntos;

  /*
    -----------------------------------------------------------------
    Section 1. Print quasi-grid data.
    -----------------------------------------------------------------
  */
  /*
    See if scanning is north->south or south->north
  */
  fprintf(grprsm, "  Number of points along a parallel varies.\n");

  ntos = ( fmod((double) isec2[10], 128.) < 64 );

  if ( ntos )
    fprintf(grprsm, "  Number of points.   Parallel. (North to South)\n");
  else
    fprintf(grprsm, "  Number of points.   Parallel. (South to North)\n");

  /*  Display number of points for each latitude */
  latcnt  = isec2[2];
  nextlat = 0;
  memset(yout, ' ', (size_t) 11);

  for ( j = 0; j < latcnt; j++ )
    {
      nextlat = nextlat + 1;
      sprintf(yout, "%4d", nextlat);

      /*       Finished?  */
      if ( nextlat > latcnt ) break;
      if ( nextlat == latcnt )
	{
	  fprintf(grprsm, " %5d                %-12s\n", isec2[nextlat+21], yout);
	  break;
	}
      /*
	Look for neighbouring latitudes with same number of points
      */
      nrepeat = 0;

    LABEL110:
      /*
	If neighbouring latitudes have same number of points
	increase the repeat count.
      */
      if ( isec2[nextlat+21+1] == isec2[nextlat+21] )
	{
          nrepeat = nrepeat + 1;
          nextlat = nextlat + 1;
	  if ( nextlat < latcnt ) goto LABEL110;
	}
      /*
	Display neighbouring latitudes with same number of points as
	'nn to mm'.
      */
      if ( nrepeat >= 1 )
	{
	  strncpy(yout+4, " to", 3);
	  sprintf(yout+7, "%5d", nextlat);
        }
      fprintf(grprsm, " %5d                %-12s\n", isec2[nextlat+21], yout);
      memset(yout, ' ', (size_t) 11);
    }
}

void gribPrintSec2DP(int *isec0, int *isec2, double *fsec2)
{
  /*

    Print the information in the Grid Description
    Section (Section 2) of decoded GRIB data.

    Input Parameters:

       isec0  - Array of decoded integers from Section 0

       isec2  - Array of decoded integers from Section 2

       fsec2  - Array of decoded floats from Section 2

    Comments:

       Only data representation types catered for are Gaussian
       grid, latitude/longitude grid, Spherical Harmonics,
       Polar stereographic and Space view perspective.


    Converted from EMOS routine GRPRS2.

       Uwe Schulzweida   MPIfM   01/04/2001

  */

  int i, ibit, iedit, ierr, iout, iresol;

  grsdef();
  /*
    -----------------------------------------------------------------
    Section 1 . Print GRIB Edition number.
    -----------------------------------------------------------------
  */
  iedit = isec0[1];
  fprintf(grprsm, " \n");
  fprintf(grprsm, " Section 2 - Grid Description Section.\n");
  fprintf(grprsm, " -------------------------------------\n");
  /*
    -----------------------------------------------------------------
    Section 2 . Print spherical harmonic data.
    -----------------------------------------------------------------
  */
  if ( isec2[0] == 50 || isec2[0] == 60 || 
       isec2[0] == 70 || isec2[0] == 80 )
    {
      fprintf(grprsm, " Data represent type = spectral     (Table 6) %9d\n", isec2[0]);
      fprintf(grprsm, " J - Pentagonal resolution parameter.         %9d\n", isec2[1]);
      fprintf(grprsm, " K - Pentagonal resolution parameter.         %9d\n", isec2[2]);
      fprintf(grprsm, " M - Pentagonal resolution parameter.         %9d\n", isec2[3]);
      fprintf(grprsm, " Representation type (Table 9)                %9d\n", isec2[4]);
      fprintf(grprsm, " Representation mode (Table 10).              %9d\n", isec2[5]);
      for (i = 7; i <= 11; i++)
        fprintf(grprsm, " Not used.                                    %9d\n", isec2[i-1]);
      fprintf(grprsm, " Number of vertical coordinate parameters.    %9d\n", isec2[11]);
      goto LABEL800;
    }
  /*
    -----------------------------------------------------------------
    Section 3 . Print Gaussian grid data.
    -----------------------------------------------------------------
  */
  if ( isec2[0] ==  4 || isec2[0] == 14 || 
       isec2[0] == 24 || isec2[0] == 34 )
    {
      fprintf(grprsm, " (Southern latitudes and Western longitudes are negative.)\n");
      fprintf(grprsm, " Data represent type = gaussian     (Table 6) %9d\n", isec2[0]);
      /*
	Quasi-regular grids introduced in Edition 1.
      */
      if ( isec2[16] == 0 || iedit < 1 )
	fprintf(grprsm, " Number of points along a parallel.           %9d\n", isec2[1]);
      else
      	printQuasi(isec2);

      fprintf(grprsm, " Number of points along a meridian.           %9d\n", isec2[2]);
      fprintf(grprsm, " Latitude of first grid point.                %9d\n", isec2[3]);
      fprintf(grprsm, " Longitude of first grid point.               %9d\n", isec2[4]);

      ibit = 8;
      iresol = isec2[5] + isec2[17] + isec2[18];
      prtbin(iresol, ibit, &iout, &ierr);

      fprintf(grprsm, " Resolution and components flag.               %8.8d\n", iout);
      fprintf(grprsm, " Latitude of last grid point.                 %9d\n", isec2[6]);
      fprintf(grprsm, " Longitude of last grid point.                %9d\n", isec2[7]);
      /*
	Print increment if given.
      */
      if ( isec2[5] == 128 )
	fprintf(grprsm, " i direction (East-West) increment.           %9d\n", isec2[8]);
      else
	fprintf(grprsm, " i direction (East-West) increment            Not given\n");

      fprintf(grprsm, " Number of parallels between pole and equator.%9d\n", isec2[9]);

      ibit = 8;
      prtbin(isec2[10], ibit, &iout, &ierr);

      fprintf(grprsm, " Scanning mode flags (Code Table 8)            %8.8d\n", iout);
      fprintf(grprsm, " Number of vertical coordinate parameters.    %9d\n", isec2[11]);
      goto LABEL800;
    }
  /*
    -----------------------------------------------------------------
    Section 4 . Print Latitude / longitude grid data.
    -----------------------------------------------------------------
  */
  if ( isec2[0] ==  0 || isec2[0] == 10 || 
       isec2[0] == 20 || isec2[0] == 30 )
    {
      fprintf(grprsm, " (Southern latitudes and Western longitudes are negative.)\n");
      fprintf(grprsm, " Data represent type = lat/long     (Table 6) %9d\n", isec2[0]);
      /*
	Quasi-regular lat/long grids also possible.
      */
      if ( isec2[16] == 0 )
	fprintf(grprsm, " Number of points along a parallel.           %9d\n", isec2[1]);
      else
        printQuasi(isec2);

      fprintf(grprsm, " Number of points along a meridian.           %9d\n", isec2[2]);
      fprintf(grprsm, " Latitude of first grid point.                %9d\n", isec2[3]);
      fprintf(grprsm, " Longitude of first grid point.               %9d\n", isec2[4]);

      ibit = 8;
      iresol = isec2[5] + isec2[17] + isec2[18];
      prtbin(iresol, ibit, &iout, &ierr);

      fprintf(grprsm, " Resolution and components flag.               %8.8d\n", iout);
      fprintf(grprsm, " Latitude of last grid point.                 %9d\n", isec2[6]);
      fprintf(grprsm, " Longitude of last grid point.                %9d\n", isec2[7]);
      /*
	Print increment if given.
      */
      if ( isec2[8] < 0 )
	fprintf(grprsm, " i direction (East-West) increment            Not given\n");
      else
	fprintf(grprsm, " i direction (East-West) increment.           %9d\n", isec2[8]);

      if ( isec2[9] < 0 )
	fprintf(grprsm, " j direction (North-South) increment          Not given\n");
      else
	fprintf(grprsm, " j direction (North-South) increment.         %9d\n", isec2[9]);
    
      ibit = 8;
      prtbin(isec2[10], ibit, &iout, &ierr);

      fprintf(grprsm, " Scanning mode flags (Code Table 8)            %8.8d\n", iout);
      fprintf(grprsm, " Number of vertical coordinate parameters.    %9d\n", isec2[11]);
      goto LABEL800;
    }
  /*
    -----------------------------------------------------------------
    Section 5 . Print polar stereographic data.
    -----------------------------------------------------------------
  */
  if ( isec2[0] == 5 )
    {
      fprintf(grprsm, " (Southern latitudes and Western longitudes are negative.)\n");
      fprintf(grprsm, " Data represent type = polar stereo (Table 6) %9d\n", isec2[0]);
      fprintf(grprsm, " Number of points along X axis.               %9d\n", isec2[1]);
      fprintf(grprsm, " Number of points along Y axis.               %9d\n", isec2[2]);
      fprintf(grprsm, " Latitude of first grid point.                %9d\n", isec2[3]);
      fprintf(grprsm, " Longitude of first grid point.               %9d\n", isec2[4]);
      ibit = 8;
      iresol = isec2[17] + isec2[18];
      prtbin(iresol, ibit, &iout, &ierr);
      fprintf(grprsm, " Resolution and components flag.               %8.8d\n", iout);
      fprintf(grprsm, " Orientation of the grid.                     %9d\n", isec2[6]);
      fprintf(grprsm, " X direction increment.                       %9d\n", isec2[8]);
      fprintf(grprsm, " Y direction increment.                       %9d\n", isec2[9]);
      ibit = 8;
      prtbin(isec2[10], ibit, &iout, &ierr);
      fprintf(grprsm, " Scanning mode flags (Code Table 8)            %8.8d\n", iout);
      fprintf(grprsm, " Number of vertical coordinate parameters.    %9d\n", isec2[11]);
      fprintf(grprsm, " Projection centre flag.                      %9d\n", isec2[12]);
      goto LABEL800;
    }
  /*
    -----------------------------------------------------------------
    Section 6 . Print Lambert conformal data.
    -----------------------------------------------------------------
  */
  if ( isec2[0] == 3 )
    {
      fprintf(grprsm, " (Southern latitudes and Western longitudes are negative.)\n");
      fprintf(grprsm, " Data represent type = Lambert      (Table 6) %9d\n", isec2[0]);
      fprintf(grprsm, " Number of points along X axis.               %9d\n", isec2[1]);
      fprintf(grprsm, " Number of points along Y axis.               %9d\n", isec2[2]);
      fprintf(grprsm, " Latitude of first grid point.                %9d\n", isec2[3]);
      fprintf(grprsm, " Longitude of first grid point.               %9d\n", isec2[4]);
      ibit = 8;
      iresol = isec2[17] + isec2[18] + isec2[5];
      prtbin(iresol, ibit, &iout, &ierr);
      fprintf(grprsm, " Resolution and components flag.               %8.8d\n", iout);
      fprintf(grprsm, " Orientation of the grid.                     %9d\n", isec2[6]);
      fprintf(grprsm, " X direction increment.                       %9d\n", isec2[8]);
      fprintf(grprsm, " Y direction increment.                       %9d\n", isec2[9]);
      ibit = 8;
      prtbin(isec2[10], ibit, &iout, &ierr);
      fprintf(grprsm, " Scanning mode flags (Code Table 8)            %8.8d\n", iout);
      fprintf(grprsm, " Number of vertical coordinate parameters.    %9d\n", isec2[11]);
      fprintf(grprsm, " Projection centre flag.                      %9d\n", isec2[12]);
      fprintf(grprsm, " Latitude intersection 1 - Latin 1 -.         %9d\n", isec2[13]);
      fprintf(grprsm, " Latitude intersection 2 - Latin 2 -.         %9d\n", isec2[14]);
      fprintf(grprsm, " Latitude of Southern Pole.                   %9d\n", isec2[19]);
      fprintf(grprsm, " Longitude of Southern Pole.                  %9d\n", isec2[20]);
      goto LABEL800;
    }
  /*
    -----------------------------------------------------------------
    Section 7 . Print space view perspective or orthographic data.
    -----------------------------------------------------------------
  */
  if ( isec2[0] == 90 )
    {
      fprintf(grprsm, " (Southern latitudes and Western longitudes are negative.)\n");
      fprintf(grprsm, " Data represent type = space/ortho  (Table 6) %9d\n", isec2[0]);
      fprintf(grprsm, " Number of points along X axis.               %9d\n", isec2[1]);
      fprintf(grprsm, " Number of points along Y axis.               %9d\n", isec2[2]);
      fprintf(grprsm, " Latitude of sub-satellite point.             %9d\n", isec2[3]);
      fprintf(grprsm, " Longitude of sub-satellite point.            %9d\n", isec2[4]);
      iresol = isec2[17] + isec2[18];
      fprintf(grprsm, " Diameter of the earth in x direction.        %9d\n", isec2[6]);
      fprintf(grprsm, " Y coordinate of sub-satellite point.         %9d\n", isec2[9]);
      ibit = 8;
      prtbin(isec2[10], ibit, &iout, &ierr);
      fprintf(grprsm, " Scanning mode flags (Code Table 8)            %8.8d\n", iout);
      fprintf(grprsm, " Number of vertical coordinate parameters.    %9d\n", isec2[11]);
      fprintf(grprsm, " Orientation of the grid.                     %9d\n", isec2[6]);
      fprintf(grprsm, " Altitude of the camera.                      %9d\n", isec2[13]);
      fprintf(grprsm, " Y coordinate of origin of sector image.      %9d\n", isec2[14]);
      fprintf(grprsm, " X coordinate of origin of sector image.      %9d\n", isec2[15]);
      goto LABEL800;
    }
  /*
    -----------------------------------------------------------------
    Section 7.5 . Print ocean data
    -----------------------------------------------------------------
  */
  /*
  if ( isec2[0] == 192 && isec1[1] == 98 )
    {
      fprintf(grprsm, " Data represent type = ECMWF ocean  (Table 6) %9d\n", isec2[0]);
      if ( isec2[1] ==  32767 )
	fprintf(grprsm, " Number of points along the first axis.       Not used\n");
      else
	fprintf(grprsm, " Number of points along the first axis.       %9d\n", isec2[1]);

      if ( isec2[2] ==  32767 )
	fprintf(grprsm, " Number of points along the second axis.      Not used\n");
      else
	fprintf(grprsm, " Number of points along the second axis.      %9d\n", isec2[2]);

      ibit = 8;
      prtbin(isec2[10], ibit, &iout, &ierr);
      fprintf(grprsm, " Scanning mode flags (Code Table 8)            %8.8d\n", iout);
      goto LABEL800;
    }
    */
  /*
    -----------------------------------------------------------------
    Section 7.6 . Print triangular data
    -----------------------------------------------------------------
  */
  if ( isec2[0] == 192 /* && isec1[1] == 78 */ )
    {
      fprintf(grprsm, " Data represent type = triangular   (Table 6) %9d\n", isec2[0]);
      fprintf(grprsm, " Number of factor 2 in factorisation of Ni.   %9d\n", isec2[1]);
      fprintf(grprsm, " Number of factor 3 in factorisation of Ni.   %9d\n", isec2[2]);
      fprintf(grprsm, " Number of diamonds (Nd).                     %9d\n", isec2[3]);
      fprintf(grprsm, " Number of triangular subdivisions of the\n");
      fprintf(grprsm, "           icosahedron (Ni).                  %9d\n", isec2[4]);
      fprintf(grprsm, " Flag for orientation of diamonds (Table A).  %9d\n", isec2[5]);
      fprintf(grprsm, " Latitude of pole point.                      %9d\n", isec2[6]);
      fprintf(grprsm, " Longitude of pole point.                     %9d\n", isec2[7]);
      fprintf(grprsm, " Longitude of the first diamond.              %9d\n", isec2[8]);
      fprintf(grprsm, " Flag for storage sequence (Table B).         %9d\n", isec2[9]);
      fprintf(grprsm, " Number of vertical coordinate parameters.    %9d\n", isec2[11]);
      goto LABEL800;
    }
  /*
    -----------------------------------------------------------------
    Drop through to here => representation type not catered for.
    -----------------------------------------------------------------
  */
  fprintf(grprsm, "GRPRS2 :Data representation type not catered for -%d\n", isec2[0]);

  goto LABEL900;
  /*
    -----------------------------------------------------------------
    Section 8 . Print vertical coordinate parameters,
                rotated grid information,
                stretched grid information, if any.
    -----------------------------------------------------------------
  */
 LABEL800:;
  /*
    Vertical coordinate parameters ...
  */
  if ( isec2[11] != 0 )
    {
      fprintf(grprsm, " \n");
      fprintf(grprsm, " Vertical Coordinate Parameters.\n");
      fprintf(grprsm, " -------------------------------\n");
      for ( i = 10; i < isec2[11]+10; i++ )
	fprintf(grprsm, "    %20.12f\n", fsec2[i]);
    }
  /*
    Rotated and stretched grids introduced in Edition 1.
  */
  if ( iedit < 1 ) goto LABEL900;
  /*
    Rotated grid information ...
  */
  if ( isec2[0] == 10 || isec2[0] == 30 || 
       isec2[0] == 14 || isec2[0] == 34 || 
       isec2[0] == 60 || isec2[0] == 80 || 
       isec2[0] == 30 )
    {
      fprintf(grprsm, " \n");
      fprintf(grprsm, " Latitude of southern pole of rotation.       %9d\n", isec2[12]);
      fprintf(grprsm, " Longitude of southern pole of rotation.      %9d\n", isec2[13]);
      fprintf(grprsm, " Angle of rotation.                     %20.10f\n", fsec2[0]);
    }
  /*
    Stretched grid information ...
  */
  if ( isec2[0] == 20 || isec2[0] == 30 || 
       isec2[0] == 24 || isec2[0] == 34 || 
       isec2[0] == 70 || isec2[0] == 80 )
    {
      fprintf(grprsm, " \n");
      fprintf(grprsm, " Latitude of pole of stretching.              %9d\n", isec2[14]);
      fprintf(grprsm, " Longitude of pole of stretching.             %9d\n", isec2[15]);
      fprintf(grprsm, " Stretching factor.                     %20.10f\n", fsec2[1]);
    }

 LABEL900:;

  return;
}

void gribPrintSec2SP(int *isec0, int *isec2, float  *fsec2sp)
{
  static char func[] = "grprs2sp";
  int inum;
  int j;
  double *fsec2;

  inum = 10 + isec2[11];

  fsec2 = (double *) malloc(inum*sizeof(double));
  if ( fsec2 == NULL ) SysError(func, "No Memory!");

  for ( j = 0; j < inum; j++ )
     fsec2[j] = fsec2sp[j];
  
  gribPrintSec2DP(isec0, isec2, fsec2);

  free(fsec2);
}

void gribPrintSec3DP(int *isec0, int *isec3, double *fsec3)
{
  /*

    Print the information in the Bit-Map Section
    (Section 3) of decoded GRIB data.

    Input Parameters:

       isec0  - Array of decoded integers from Section 0

       isec3  - Array of decoded integers from Section 3

       fsec3  - Array of decoded floats from Section 3


    Converted from EMOS routine GRPRS3.

       Uwe Schulzweida   MPIfM   01/04/2001

  */

  grsdef();

  fprintf(grprsm, " \n");
  fprintf(grprsm, " Section 3 - Bit-map Section.\n");
  fprintf(grprsm, " -------------------------------------\n");

  if ( isec3[0] != 0 )
    fprintf(grprsm, " Predetermined bit-map number.                %9d\n", isec3[0]);
  else
    fprintf(grprsm, " No predetermined bit-map.\n");

  fprintf(grprsm, " Missing data value for integer data.    %14d\n", isec3[1]);

  fprintf(grprsm, " Missing data value for real data. %20.6g\n", fsec3[1]);
}

void gribPrintSec3SP(int *isec0, int *isec3, float  *fsec3sp)
{
  double fsec3[2];

  fsec3[0] = fsec3sp[0];
  fsec3[1] = fsec3sp[1];
  
  gribPrintSec3DP(isec0, isec3, fsec3);
}

void gribPrintSec4DP(int *isec0, int *isec4, double *fsec4)
{
  /*

    Print the information in the Binary Data Section
    (Section 4) of decoded GRIB data.

    Input Parameters:

       isec0  - Array of decoded integers from Section 0

       isec4  - Array of decoded integers from Section 4

       fsec4  - Array of decoded floats from Section 4


    Converted from EMOS routine GRPRS4.

       Uwe Schulzweida   MPIfM   01/04/2001

  */

  int inum;
  int j;

  grsdef();

  /*
    -----------------------------------------------------------------
    Section 1 . Print integer information from isec4.
    -----------------------------------------------------------------
  */
  fprintf(grprsm, " \n");
  fprintf(grprsm, " Section 4 - Binary Data  Section.\n");
  fprintf(grprsm, " -------------------------------------\n");

  fprintf(grprsm, " Number of data values coded/decoded.         %9d\n", isec4[0]);
  fprintf(grprsm, " Number of bits per data value.               %9d\n", isec4[1]);
  fprintf(grprsm, " Type of data       (0=grid pt, 128=spectral).%9d\n", isec4[2]);
  fprintf(grprsm, " Type of packing    (0=simple, 64=complex).   %9d\n", isec4[3]);
  fprintf(grprsm, " Type of data       (0=float, 32=integer).    %9d\n", isec4[4]);
  fprintf(grprsm, " Additional flags   (0=none, 16=present).     %9d\n", isec4[5]);
  fprintf(grprsm, " Reserved.                                    %9d\n", isec4[6]);
  fprintf(grprsm, " Number of values   (0=single, 64=matrix).    %9d\n", isec4[7]);
  fprintf(grprsm, " Secondary bit-maps (0=none, 32=present).     %9d\n", isec4[8]);
  fprintf(grprsm, " Values width       (0=constant, 16=variable).%9d\n", isec4[9]);
  /*
    If complex packing ..
  */
  if ( isec4[3] == 64 )
    {
      if ( isec4[2] == 128 )
	{
	  fprintf(grprsm, " Byte offset of start of packed data (N).     %9d\n", isec4[15]);
	  fprintf(grprsm, " Power (P * 1000).                            %9d\n", isec4[16]);
	  fprintf(grprsm, " Pentagonal resolution parameter J for subset.%9d\n", isec4[17]);
	  fprintf(grprsm, " Pentagonal resolution parameter K for subset.%9d\n", isec4[18]);
	  fprintf(grprsm, " Pentagonal resolution parameter M for subset.%9d\n", isec4[19]);
	}
      else
	{
	  fprintf(grprsm, " Bits number of 2nd order values    (none=>0).%9d\n", isec4[10]);
	  fprintf(grprsm, " General extend. 2-order packing (0=no,8=yes).%9d\n", isec4[11]);
	  fprintf(grprsm, " Boustrophedonic ordering        (0=no,4=yes).%9d\n", isec4[12]);
	  fprintf(grprsm, " Spatial differencing order          (0=none).%9d\n", isec4[13]+isec4[14]);
        }
    }
  /*
    Number of non-missing values
  */
  if ( isec4[20] != 0 )
    fprintf(grprsm, " Number of non-missing values                 %9d\n", isec4[20]);
  /*
    Information on matrix of values , if present.
  */
  if ( isec4[7] == 64 )
    {
      fprintf(grprsm, " First dimension (rows) of each matrix.       %9d\n", isec4[49]);
      fprintf(grprsm, " Second dimension (columns) of each matrix.   %9d\n", isec4[50]);
      fprintf(grprsm, " First dimension coordinate values definition.%9d\n", isec4[51]);
      fprintf(grprsm, " (Code Table 12)\n");
      fprintf(grprsm, " NC1 - Number of coefficients for 1st dimension.%7d\n", isec4[52]);
      fprintf(grprsm, " Second dimension coordinate values definition.%8d\n", isec4[53]);
      fprintf(grprsm, " (Code Table 12)\n");
      fprintf(grprsm, " NC2 - Number of coefficients for 2nd dimension.%7d\n", isec4[54]);
      fprintf(grprsm, " 1st dimension physical signifance (Table 13). %8d\n", isec4[55]);
      fprintf(grprsm, " 2nd dimension physical signifance (Table 13).%8d\n", isec4[56]);
    }
  /*
    -----------------------------------------------------------------
    Section 2. Print values from fsec4.
    -----------------------------------------------------------------
  */

  inum = isec4[0];
  if ( inum <  0 ) inum = - inum;
  if ( inum > 20 ) inum = 20;
  /*
    Print first inum values.
  */
  fprintf(grprsm, " \n");
  fprintf(grprsm, " First %4d data values.\n", inum);

  if ( isec4[4] == 0 )
    {
      /*
	Print real values ...
      */
      for ( j = 0; j < inum; j++ )
	{
	  if ( fabs(fsec4[j]) > 0 )
	    {
	      if ( fabs(fsec4[j]) >= 0.1 && fabs(fsec4[j]) <= 1.e8 )
		fprintf(grprsm, " %#16.8G    \n", fsec4[j]);
	      else
		fprintf(grprsm, " %#20.8E\n", fsec4[j]);
	    }
	  else
	    fprintf(grprsm, " %#16.0f    \n", fabs(fsec4[j]));
	}
    }
  else
    {
      /*
	Print integer values ...
      */
      fprintf(grprsm, " Print of integer values not supported\n");
      /*
        CALL SETPAR(IBIT,IDUM,IDUM)
        DO 212 J=1,INUM
           INSPT = 0
           CALL INXBIT(IVALUE,1,INSPT,FSEC4(J),1,IBIT,IBIT,'C',IRET)
           WRITE (*,9033) IVALUE
 9033 FORMAT(' ',I15)
  212   CONTINUE
      ENDIF
      */
    }
}

void gribPrintSec4SP(int *isec0, int *isec4, float  *fsec4sp)
{
  int inum;
  int j;
  double fsec4[20];

  inum = isec4[0];
  if ( inum <  0 ) inum = -inum;
  if ( inum > 20 ) inum = 20;

  for ( j = 0; j < inum; j++ ) fsec4[j] = fsec4sp[j];
  
  gribPrintSec4DP(isec0, isec4, fsec4);
}

void gribPrintSec4Wave(int *isec4)
{
  /*

    Print the wave coordinate information in the Binary Data
    Section (Section 4) of decoded GRIB data.

    Input Parameters:

       isec4 - Array of decoded integers from Section 4

    Comments:

       Wave coordinate information held in isec4 are 32-bit floats,
       hence the PTEMP and NTEMP used for printing are 4-byte variables.


    Converted from EMOS routine GRPRS4W.

       Uwe Schulzweida   MPIfM   01/04/2001

  */
  int    jloop;
  int    ntemp[100];
  float *ptemp;

  grsdef();

  /*
    -----------------------------------------------------------------
    Section 1 . Print integer information from isec4.
    -----------------------------------------------------------------
  */
  fprintf(grprsm, " Coefficients defining first dimension coordinates:\n");
  for ( jloop = 0; jloop < isec4[52]; jloop++ )
    {
      ntemp[jloop] = isec4[59 + jloop];
      ptemp = (float *) &ntemp[jloop];
      fprintf(grprsm, "%20.10f\n", *ptemp);
    }
  fprintf(grprsm, " Coefficients defining second dimension coordinates:\n");
  for ( jloop = 0; jloop < isec4[54]; jloop++ )
    {
      ntemp[jloop] = isec4[59 + isec4[52] + jloop];
      ptemp = (float *) &ntemp[jloop];
      fprintf(grprsm, "%20.10f\n", *ptemp);
    }
}


static int encode_init = 0;
static double  IntPower2[158];
static double rIntPower2[127];

static void init_encode(void)
{
  int jloop;

  IntPower2[0] = 1.0;

#if defined (__uxp__)
#pragma loop scalar  /* bug on vpp */
#endif
  for ( jloop = 1; jloop < 158; jloop++ )
    {
      IntPower2[jloop] = IntPower2[jloop-1] * 2.0;
    }

  rIntPower2[0] = 1.0;

#if defined (__uxp__)
#pragma loop scalar  /* bug on vpp */
#endif
  for ( jloop = 1; jloop < 127; jloop++ )
    {
      rIntPower2[jloop] = rIntPower2[jloop-1] * 0.5;
    }
  
  encode_init = 1;
}


int  BitsPerInt = sizeof(int) * 8;

#define PutnZero(n) \
{ \
  int i; \
  for ( i = z; i < z+n; i++ ) lGrib[i] = 0; \
  z += n; \
}

#define Put1Byte(Value)  (lGrib[z++] = (Value))
#define Put2Byte(Value) ((lGrib[z++] = (Value) >>  8), \
                         (lGrib[z++] = (Value)))
#define Put3Byte(Value) ((lGrib[z++] = (Value) >> 16), \
                         (lGrib[z++] = (Value) >>  8), \
                         (lGrib[z++] = (Value)))

#define Put1Real(Value) \
{ \
  int Exponent, Mantissa; \
  int one = 1; \
  confp3(Value, &Exponent, &Mantissa, BitsPerInt, one); \
  Put1Byte(Exponent); \
  Put3Byte(Mantissa); \
}

/* GRIB block 0 - indicator block */

void encodeIS(GRIBPACK *lGrib, int *gribLen)
{
  int z = *gribLen;

  lGrib[0] = 'G';
  lGrib[1] = 'R';
  lGrib[2] = 'I';
  lGrib[3] = 'B';

  /* 
   * lGrib[4]-lGrib[6] contains full length of grib record. 
   * included before finished CODEGB
   */

  z = 7;   
  Put1Byte(1); 
  z = 8;

  *gribLen = z;
}

/* GRIB block 5 - end block */

void encodeES(GRIBPACK *lGrib, int *gribLen)
{
  int z = *gribLen;

  lGrib[z++] = '7';
  lGrib[z++] = '7';
  lGrib[z++] = '7';
  lGrib[z++] = '7';

  lGrib[4] = z >> 16;
  lGrib[5] = z >>  8;
  lGrib[6] = z;

  while ( z & 7 ) lGrib[z++] = 0;

  *gribLen = z;
}

/* GRIB block 1 - product definition block. */

#define DWD_extension_253_len 37
#define DWD_extension_254_len 26

int getLocalExtLen(int *isec1)
{
  int extlen = 0;

  if ( ISEC1_LocalFLag )
    {
      if ( ISEC1_CenterID == 78 ) 
	{
	  if ( isec1[36] == 254 ) 
	    {
	      extlen = DWD_extension_254_len;
	    }
	  else if ( isec1[36] == 253 )
	    { 
	      extlen = DWD_extension_253_len;
	    }
	}
    }

  return (extlen);
}

int getPdsLen(int *isec1)
{
  int pdslen = 28;

  pdslen += getLocalExtLen(isec1);

  return (pdslen);
}

void encodePDS_DWD_local_Extension_254(GRIBPACK *lGrib, int *zs, int *isec1)
{
  int i, localextlen, isvn;
  int z = *zs;

  localextlen = getLocalExtLen(isec1);
  for ( i = 0; i < localextlen-2; i++ )
    {
      Put1Byte(isec1[24+i]);
    }

  isvn = isec1[49] << 15 | isec1[48]; /* DWD experiment identifier    */
  Put2Byte(isvn);             /* DWD run type (0=main, 2=ass, 3=test) */

  *zs = z;
}

void encodePDS_DWD_local_Extension_253(GRIBPACK *lGrib, int *zs, int *isec1)
{
  int i, localextlen, isvn;
  int z = *zs;

  localextlen = DWD_extension_254_len;
  for ( i = 0; i < localextlen-2; i++ )
    {
      Put1Byte(isec1[24+i]);
    }

  isvn = isec1[49] << 15 | isec1[48]; /* DWD experiment identifier    */
  Put2Byte(isvn);             /* DWD run type (0=main, 2=ass, 3=test) */
  Put1Byte(isec1[50]);        /* 55 User id, specified by table       */
  Put2Byte(isec1[51]);        /* 56 Experiment identifier             */
  Put2Byte(isec1[52]);        /* 58 Ensemble identification by table  */
  Put2Byte(isec1[53]);        /* 60 Number of ensemble members        */
  Put2Byte(isec1[54]);        /* 62 Actual number of ensemble member  */
  Put1Byte(isec1[55]);        /* 64 Model major version number        */ 
  Put1Byte(isec1[56]);        /* 65 Model minor version number        */ 

  *zs = z;
}

/* GRIB BLOCK 1 - PRODUCT DESCRIPTION SECTION */

void encodePDS(GRIBPACK *lpds, int pdsLen, int *isec1)
{
  GRIBPACK *lGrib = lpds;
  int z = 0;

  Put3Byte(pdsLen);               /*  0 Length of Block 1        */
  Put1Byte(ISEC1_CodeTable);      /*  3 Local table number       */
  Put1Byte(ISEC1_CenterID);       /*  4 Identification of centre */
  Put1Byte(ISEC1_ModelID);        /*  5 Identification of model  */
  Put1Byte(ISEC1_GridDefinition); /*  6 Grid definition          */
  Put1Byte(ISEC1_Sec2Or3Flag);    /*  7 Block 2 included         */
  Put1Byte(ISEC1_Parameter);      /*  8 Parameter Code           */
  Put1Byte(ISEC1_LevelType);      /*  9 Type of level            */
  if ( (ISEC1_LevelType !=  20) &&
       (ISEC1_LevelType !=  99) &&
       (ISEC1_LevelType != 100) &&
       (ISEC1_LevelType != 103) &&
       (ISEC1_LevelType != 105) &&
       (ISEC1_LevelType != 107) &&
       (ISEC1_LevelType != 109) &&
       (ISEC1_LevelType != 111) &&
       (ISEC1_LevelType != 113) &&
       (ISEC1_LevelType != 115) &&
       (ISEC1_LevelType != 117) &&
       (ISEC1_LevelType != 125) &&
       (ISEC1_LevelType != 127) &&
       (ISEC1_LevelType != 160) &&
       (ISEC1_LevelType != 210) )
    {
      Put1Byte(ISEC1_Level1);
      Put1Byte(ISEC1_Level2);
    }
  else
    {
      Put2Byte(ISEC1_Level1);     /* 10 Level                    */    
    }
  Put1Byte(ISEC1_Year);           /* 12 Year of Century          */
  Put1Byte(ISEC1_Month);          /* 13 Month                    */
  Put1Byte(ISEC1_Day);            /* 14 Day                      */
  Put1Byte(ISEC1_Hour);           /* 15 Hour                     */
  Put1Byte(ISEC1_Minute);         /* 16 Minute                   */
  if ( ISEC1_AvgNum > 0 )
    {
      Put1Byte(2);                    /* 17 Time unit                */
      Put1Byte(0);                    /* 18 Time 1                   */
      Put1Byte(0);                    /* 19 Time 2                   */
      Put1Byte(3);                    /* 20 Timerange flag           */
      Put2Byte(ISEC1_AvgNum);     /* 21 Average                  */
    }
  else
    {
      Put1Byte(ISEC1_TimeUnit);   /* 17 Time unit                */
      if ( ISEC1_TimeRange == 10 )
	{
	  Put1Byte(ISEC1_TimePeriod1);
	  Put1Byte(ISEC1_TimePeriod2);
	}
      else if ( ISEC1_TimeRange == 113 || ISEC1_TimeRange ==   0 )
	{
	  Put1Byte(ISEC1_TimePeriod1);
          Put1Byte(0);
	}
      else if ( ISEC1_TimeRange ==   4 || ISEC1_TimeRange ==   2 )
	{
          Put1Byte(0);
	  Put1Byte(ISEC1_TimePeriod2);
	}
      else
	{
          Put1Byte(0);
          Put1Byte(0); 
	}
      Put1Byte(ISEC1_TimeRange);      /* 20 Timerange flag           */
      Put2Byte(ISEC1_AvgNum);         /* 21 Average                  */
    }

  Put1Byte(ISEC1_AvgMiss);            /* 23 Missing from averages    */
  Put1Byte(ISEC1_Century);            /* 24 Century                  */
  Put1Byte(ISEC1_SubCenterID);        /* 25 Subcenter                */
  Put2Byte(ISEC1_DecScaleFactor);     /* 26 Decimal scale factor     */

  if ( ISEC1_LocalFLag )
    {
      if ( ISEC1_CenterID == 78 ) 
	{
	  if ( isec1[36] == 254 ) 
	    {
	      encodePDS_DWD_local_Extension_254(lGrib, &z, isec1);
	    }
	  else if ( isec1[36] == 253 )
	    { 
	      encodePDS_DWD_local_Extension_253(lGrib, &z, isec1);
	    }
	}
      else
	{
	  int i, localextlen;
	  localextlen = getLocalExtLen(isec1);
	  for ( i = 0; i < localextlen; i++ )
	    {
	      Put1Byte(isec1[24+i]);
	    }
	}
    }
}

/* GRIB BLOCK 2 - GRID DESCRIPTION SECTION */

void encodeGDS(GRIBPACK *lGrib, int *gribLen, int *isec1, int *isec2, double *fsec2)
{
  int z = *gribLen;
  int Exponent, Mantissa;
  int one = 1;
  int i;
  int ival;
  int BlockLength = 32;

  BlockLength += ISEC2_NumVCP * 4;
  if ( ISEC2_Reduced ) BlockLength += 2 * ISEC2_NumLat;

  Put3Byte(BlockLength);        /*  0- 2 Length of Block 2 Byte 0 */

  Put1Byte(ISEC2_NumVCP);       /*  3    NV */
  if ( ISEC2_NumVCP || ISEC2_Reduced )
    Put1Byte(          33);     /*  4    PV */
  else
    Put1Byte(         255);     /*  4   255 */

  Put1Byte(ISEC2_GridType);     /*  5  LatLon=0 Gauss=4 Spectral=50 */

  if ( ISEC2_GridType == GTYPE_SPECTRAL )
    {
      Put2Byte(ISEC2_PentaJ);   /*  6- 7 Pentagonal resolution J  */
      Put2Byte(ISEC2_PentaK);   /*  8- 9 Pentagonal resolution K  */
      Put2Byte(ISEC2_PentaM);   /* 10-11 Pentagonal resolution M  */
      Put1Byte(ISEC2_RepType);  /* 12    Representation type      */
      Put1Byte(ISEC2_RepMode);  /* 13    Representation mode      */
      PutnZero(18);             /* 14-31 reserved                 */
    }
  else if ( ISEC2_GridType == GTYPE_TRIANGULAR )
    {
      Put2Byte(ISEC2_TRI_NI2);
      Put2Byte(ISEC2_TRI_NI3);
      Put3Byte(ISEC2_TRI_ND);
      Put3Byte(ISEC2_TRI_NI);
      Put1Byte(ISEC2_TRI_AFlag);
      ival = ISEC2_TRI_LatPP;
      if ( ival < 0 ) ival = 8388608 - ival;
      Put3Byte(ival);
      ival = ISEC2_TRI_LonPP;
      if ( ival < 0 ) ival = 8388608 - ival;
      Put3Byte(ival);
      ival = ISEC2_TRI_LonMPL;
      if ( ival < 0 ) ival = 8388608 - ival;
      Put3Byte(ival);
      Put1Byte(ISEC2_TRI_BFlag);
      PutnZero(5);
    }
  else
    {
      Put2Byte(ISEC2_NumLon);        /*  6- 7 Longitudes               */

      Put2Byte(ISEC2_NumLat);          /*  8- 9 Latitudes                */
      ival = ISEC2_FirstLat;
      if ( ival < 0 ) ival = 8388608 - ival;
      Put3Byte(ival);                  /* 10-12 Latitude  of Origin      */
      ival = ISEC2_FirstLon;
      if ( ival < 0 ) ival = 8388608 - ival;
      Put3Byte(ival);                  /* 13-15 Longitude of Origin      */
      Put1Byte(ISEC2_ResFlag);         /* 16    Resolution flag          */

      ival = ISEC2_LastLat;
      if ( ival < 0 ) ival = 8388608 - ival;
      Put3Byte(ival);                  /* 17-19 Latitude  of Extreme     */
      ival = ISEC2_LastLon;
      if ( ival < 0 ) ival = 8388608 - ival;
      Put3Byte(ival);                  /* 20-22 Longitude of Extreme     */
      ival = ISEC2_LonIncr;
      Put2Byte(ival);                  /* 23-24 i - direction increment  */
      if ( ISEC2_GridType == GTYPE_GAUSSIAN )
	Put2Byte(ISEC2_NumPar);        /* 25-26 Latitudes Pole->Equator  */
      else
	{
	  ival = ISEC2_LatIncr;
	  Put2Byte(ival);              /* 25-26 j - direction increment  */
	}
      Put1Byte(ISEC2_ScanFlag);        /* 27    Scanning mode            */
      PutnZero(4);                     /* 28-31 reserved                 */
    }

#if defined (SX) || defined (ES)
#pragma vdir novector     /* vectorization gives wrong results on NEC */
#endif
  for ( i = 0; i < ISEC2_NumVCP; ++i )
    {
      confp3(fsec2[10+i], &Exponent, &Mantissa, BitsPerInt, one);
      Put1Byte(Exponent);
      Put3Byte(Mantissa);
    }

  if ( ISEC2_Reduced )
    for ( i = 0; i < ISEC2_NumLat; i++ ) Put2Byte(ISEC2_RowLon(i));

  *gribLen = z;
}

/* GRIB BLOCK 3 - BIT MAP SECTION */

void encodeBMS(GRIBPACK *lGrib, int *gribLen, double *fsec3, int *isec4, double *data, int *datasize)
{
  static char func[] = "encodeBMS";
  GRIBPACK *bitmap;
  int bitmapSize;
  int i;
  int bmsLen, bmsUnusedBits;
  int fsec4size;
  int z = *gribLen;
  unsigned int *imask;
  /*  unsigned int c, imask; */

  bitmapSize = ISEC4_NumValues;
  bitmap = &lGrib[z+6];
  fsec4size = 0;
  /*
  for ( i = 0; i < bitmapSize/8; i++ ) bitmap[i] = 0;

  for ( i = 0; i < bitmapSize; i++ )
    {
      if ( data[i] != FSEC3_MissVal )
	{
	  data[fsec4size++] = data[i];
	  bitmap[i/8] |= 1<<(7-(i&7));
	}
    }
  */

  imask = (unsigned int *) malloc(bitmapSize*sizeof(int));
  memset(imask, 0, bitmapSize*sizeof(int));

#if defined (CRAY)
#pragma _CRI ivdep
#endif
#if defined (SX) || defined (ES)
#pragma vdir nodep
#endif
#ifdef __uxpch__
#pragma loop novrec
#endif
  for ( i = 0; i < bitmapSize; i++ )
    {
      if ( data[i] != FSEC3_MissVal )
	{
	  data[fsec4size++] = data[i];
	  imask[i] = 1;
	}
    }

#if defined (CRAY)
#pragma _CRI ivdep
#endif
#if defined (SX) || defined (ES)
#pragma vdir nodep
#endif
#ifdef __uxpch__
#pragma loop novrec
#endif
  for ( i = 0; i < bitmapSize/8; i++ )
    {
      bitmap[i] = (imask[i*8+0] << 7) | (imask[i*8+1] << 6) |
	          (imask[i*8+2] << 5) | (imask[i*8+3] << 4) |
	          (imask[i*8+4] << 3) | (imask[i*8+5] << 2) |
	          (imask[i*8+6] << 1) | (imask[i*8+7]);
    }

  free(imask);

  bmsLen = bitmapSize/8 + 6;
  bmsUnusedBits = bitmapSize%8;

  Put3Byte(bmsLen);   /*  0- 2 Length of Block 3 Byte 0 */
  Put1Byte(bmsUnusedBits);
  Put2Byte(0);

  *gribLen += bmsLen;

  *datasize = fsec4size;
}

/* GRIB BLOCK 4 - BINARY DATA SECTION */

int encodeBDS(GRIBPACK *lGrib, int *gribLen, int *isec2, int *isec4, int datasize, double data[])
{
  /* Uwe Schulzweida, 11/04/2003 : Check that number of bits per value is not exceeded */
  /* Uwe Schulzweida,  6/05/2003 : Copy result to fpval to prevent integer overflow */

  static char func[] = "encodeBDS";
  static int lwarn_cplx = TRUE;
  int z = *gribLen;
  int  i, jloop;
  int  BlockLength, PackStart, Flag;
  int  binscale = 0, pval;
  int  byte_per_value;
  int  nbpv, ibits = BitsPerInt;
  int  max_nbpv_pow2;
  double fpval;
  double factor = 1, fmin, fmax, zref;
  double range, rangec;
  double jpepsln = 1.0e-12;     /* -----> tolerance used to check equality     */
                                /*        of floating point numbers - needed   */
		                /*        on some platforms (eg vpp700, linux) */

  if ( ! encode_init ) init_encode();

  if ( isec4[3] == 64 && lwarn_cplx )
    {
      lwarn_cplx = FALSE;
      Message(func, "complex packing of spectral data unsupported, using simple packing!");
    }

  byte_per_value = ISEC4_NumBits >> 3;

  if ( ISEC2_GridType == GTYPE_SPECTRAL )
    {
      PackStart   = 1;
      Flag        = 128 + 8;
      BlockLength = 16 + byte_per_value*(datasize - 1);
    }
  else
    {
      PackStart   = 0;
      Flag        = 8;
      BlockLength = 12 + byte_per_value*datasize;
    }

  fmin = fmax = data[PackStart];

#if   defined (CRAY)
#pragma _CRI ivdep
#elif defined (SX) || defined (ES)
#pragma vdir nodep
#elif defined (__uxp__)
#pragma loop novrec
#endif
  for ( i = PackStart+1; i < datasize; ++i )
    {
      if ( fmin > data[i] ) fmin = data[i];
      if ( fmax < data[i] ) fmax = data[i];
    }

  zref = fmin;

  /*
    Adjust number of bits per value if full integer length to
    avoid hitting most significant bit (sign bit).
  */
  nbpv = ISEC4_NumBits;
  /* if( nbpv == ibits ) nbpv = nbpv - 1; */
  /*
    Calculate the binary scaling factor to spread the range of
    values over the number of bits per value.
    Limit scaling to 2**-126 to 2**127 (using IEEE 32-bit floats
    as a guideline).           
  */
  range = fabs(fmax - fmin);
  if ( fabs(fmin) < FLT_MIN ) fmin = 0;
  /*
    Have to allow tolerance in comparisons on some platforms
    (eg vpp700 and linux), such as 0.9999999999999999 = 1.0,
    to avoid clipping ranges which are a power of 2.
  */
  if ( range <= jpepsln )
    {
      binscale = 0;
    }
  else if ( (fmin != 0.0) && (fabs(range/fmin) <= jpepsln) )
    {
      binscale = 0;
    }
  else if ( fabs(range-1.0) <= jpepsln )
    {
      binscale = 1 - nbpv;
    }
  else if ( range > 1.0 )
    {
      rangec = range + jpepsln;
      for ( jloop = 1; jloop < 128; jloop++ )
	{
	  if ( IntPower2[jloop] > rangec ) break;
	}
      if ( jloop == 128 )
	{
	  gprintf(func, "Problem calculating binary scale value for encode");
	  gprintf(func, "> range %g rangec %g fmin %g fmax %g", range, rangec, fmin, fmax);
	  return (707);
	}
      else
	{
	  binscale = jloop - nbpv;
	}
    }
  else
    {
      rangec = range - jpepsln;
      for ( jloop = 1; jloop < 127; jloop++ )
	{
	  if ( rIntPower2[jloop] < rangec ) break;
	}
      if ( jloop == 127 )
	{
	  gprintf(func, "Problem calculating binary scale value for encode");
	  gprintf(func, "< range %g rangec %g fmin %g fmax %g", range, rangec, fmin, fmax);
	  return (707);
	}
      else
	{
	  binscale = 1 - jloop - nbpv;
	}
    }

  if ( binscale != 0 )
    {
      if ( binscale < 0 ) factor =  IntPower2[-binscale];
      else                factor = rIntPower2[ binscale];
    }

  max_nbpv_pow2 = (int) IntPower2[nbpv] - 1;

  ref2ibm(&zref, BitsPerInt);

  Put3Byte(BlockLength);      /*  0-2 Length of Block 4        */
  Put1Byte(Flag);             /*  3   Flag & Unused bits       */
  if (binscale < 0) binscale = 32768 - binscale;
  Put2Byte(binscale);         /*  4-5 Scale factor             */
  Put1Real(zref);             /*  6-9 Reference value          */
  Put1Byte(nbpv);             /*   10 Packing size             */

  if ( PackStart ) Put1Real(data[0]);

  if      ( ISEC4_NumBits ==  8 )
    {
#if   defined (CRAY)
#pragma _CRI ivdep
#elif defined (SX) || defined (ES)
#pragma vdir nodep
#elif defined (__uxp__)
#pragma loop novrec
#endif
      for ( i = PackStart ; i < datasize; i++ )
	{
          fpval = ((data[i] - zref) * factor + 0.5);
	  if ( fpval > max_nbpv_pow2 ) fpval = max_nbpv_pow2;
	  if ( fpval < 0 ) fpval = 0;
	  pval = (int) fpval;
	  lGrib[z  ] = pval;
          z++;
	}
    }
  else if ( ISEC4_NumBits == 16 )
    {
#if   defined (CRAY)
#pragma _CRI ivdep
#elif defined (SX) || defined (ES)
#pragma vdir nodep
#elif defined (__uxp__)
#pragma loop novrec
#endif
      for ( i = PackStart ; i < datasize; i++ )
	{
	  fpval = ((data[i] - zref) * factor + 0.5);
	  if ( fpval > max_nbpv_pow2 ) fpval = max_nbpv_pow2;
	  if ( fpval < 0 ) fpval = 0;
	  pval = (int) fpval;
	  lGrib[z  ] = pval >>  8;
          lGrib[z+1] = pval;
	  z += 2;
	}
    }
  else if ( ISEC4_NumBits == 24 )
    {
#if   defined (CRAY)
#pragma _CRI ivdep
#elif defined (SX) || defined (ES)
#pragma vdir nodep
#elif defined (__uxp__)
#pragma loop novrec
#endif
      for ( i = PackStart ; i < datasize; i++ )
	{
          fpval = ((data[i] - zref) * factor + 0.5);
	  if ( fpval > max_nbpv_pow2 ) fpval = max_nbpv_pow2;
	  if ( fpval < 0 ) fpval = 0;
	  pval = (int) fpval;
          lGrib[z  ] =  pval >> 16;
          lGrib[z+1] =  pval >>  8;
          lGrib[z+2] =  pval;
          z += 3;
	}
    }
  else
    {
      Error(func, "Unimplemented packing factor %d", ISEC4_NumBits);
    }

  Put1Byte(0);              /*  Fillbyte                     */

  *gribLen = z;

  return (0);
}

void gribEncode(int *isec0, int *isec1, int *isec2, double *fsec2, int *isec3,
		double *fsec3, int *isec4, double *fsec4, int klenp, int *kgrib,
		int kleng, int *kword, int efunc, int *kret)
{
  static char func[] = "gribEncode";
  int gribLen = 0; /* Counter of GRIB length for output */ 
  long isLen, pdsLen;
  GRIBPACK *lpds;
  char *CGrib;
  int iret;
  int fsec4size = 0;
  int bmsIncluded;
  size_t len;
  GRIBPACK *lGrib;

  CGrib = (char *) kgrib;

  len = (klenp * 4 + 2000); /* only for maximum 24 bit packing */

  lGrib = (GRIBPACK *) malloc(len*sizeof(GRIBPACK));
  if ( lGrib == NULL ) SysError(func, "No Memory!");

  isLen = 8;
  encodeIS(lGrib, &gribLen);
  lpds = &lGrib[isLen];
  pdsLen = getPdsLen(isec1);

  encodePDS(lpds, pdsLen, isec1);
  gribLen += pdsLen;

  encodeGDS(lGrib, &gribLen, isec1, isec2, fsec2);
  /*
    ----------------------------------------------------------------
    BMS Bit-Map Section Section (Section 3)
    ----------------------------------------------------------------
  */ 
  bmsIncluded = ISEC1_Sec2Or3Flag & 64;

  if ( bmsIncluded )
    {
      encodeBMS(lGrib, &gribLen, fsec3, isec4, fsec4, &fsec4size);
    }
  else
    {
      fsec4size = ISEC4_NumValues;
    }

  iret = encodeBDS(lGrib, &gribLen, isec2, isec4, fsec4size, fsec4);
  encodeES(lGrib, &gribLen);

  if ( (size_t) gribLen > len )
    Error(func, "lGrib buffer to small! len = %d  gribLen = %d", len, gribLen);

  if ( (size_t) gribLen > kleng*sizeof(int) )
    Error(func, "kgrib buffer to small! kleng = %d  gribLen = %d", kleng, gribLen);

  (void) PACK_GRIB(lGrib, (char *)CGrib, gribLen, -1L);

  free(lGrib);

  ISEC0_GRIB_Len     = gribLen;
  ISEC0_GRIB_Version = 1;

  *kword = gribLen / sizeof(int);
  if ( (size_t) gribLen != *kword * sizeof(int) ) *kword += 1;

  *kret = 0;
}



int gribVersion(unsigned char *is, size_t buffersize)
{
  static char func[] = "gribVersion";

  if ( buffersize < 8 )
    Error(func, "buffer to small (current size %d)!\n", (int) buffersize);

  return (IS_GRIB_Version);
}

double GET_Real(unsigned char *grib)
{
  int iexp, imant;

  iexp  = GET_UINT1(grib[0]);
  imant = GET_UINT3(grib[1], grib[2], grib[3]);

  return (decfp2(iexp, imant));
}

int decodeIS(unsigned char *is, int *isec0, int *iret)
{
  static char func[] = "decodeIS";
  int isLen = 0;
  int grib1offset;
  int lgrib = FALSE, lbudg = FALSE, ltide = FALSE;

  /*
    Octets 1 - 4 : The letters G R I B.
    Four 8 bit fields.
  */
  /*
    Check letters -> GRIB, BUDG or TIDE.
  */
  /*
    Check that 'GRIB' is found where expected.
  */
  if ( is[0] == 'G' && is[1] == 'R' && is[2] == 'I' && is[3] == 'B' ) lgrib = TRUE;
  /*
    ECMWF pseudo-grib data uses 'BUDG' and 'TIDE'.
  */
  if ( is[0] == 'B' && is[1] == 'U' && is[2] == 'D' && is[3] == 'G' ) lbudg = TRUE;
  if ( is[0] == 'T' && is[1] == 'I' && is[2] == 'D' && is[3] == 'E' ) ltide = TRUE;
  /*
    Data is not GRIB or pseudo-grib.
  */
  if ( lgrib == FALSE && lbudg == FALSE && ltide == FALSE )
    {
      *iret = 305;
      gprintf(func, "Input data is not GRIB or pseudo-grib.");
      gprintf(func, "Return code = %d", *iret);
    }
  if ( lbudg == TRUE || ltide == TRUE )
    {
      *iret = 305;
      gprintf(func, "Pseudo-grib data unsupported.");
      gprintf(func, "Return code = %d", *iret);
    }

  /*
    Octets 5 - 7 : Length of message.
    One 24 bit field.
  */
  ISEC0_GRIB_Len = IS_GRIB_Len;
  /*
    Octet 8 : GRIB Edition Number.
    One 8 bit field.
  */
  ISEC0_GRIB_Version = IS_GRIB_Version;

  if ( ISEC0_GRIB_Version > 1 )
    Error(func, "GRIB version %d unsupported!", ISEC0_GRIB_Version);

  grib1offset = ISEC0_GRIB_Version * 4;

  isLen = 4 + grib1offset;

  return (isLen);
}

void decodePDS_DWD_local_Extension_254(unsigned char *pds, int *isec1)
{
  int i, isvn;

  isec1[36] = GET_UINT1(pds[40]); /* extension identifier */
  for ( i = 0; i < 11; i++) 
    { 
      isec1[37+i] =  GET_UINT1(pds[41+i]);
    } 

  isvn = GET_UINT2(pds[52],pds[53]);
  
  isec1[48] =  isvn % 0x8000;              /* DWD experiment identifier            */
  isec1[49] =  isvn >> 15;                 /* DWD run type (0=main, 2=ass, 3=test) */

}



void decodePDS_DWD_local_Extension_253(unsigned char *pds, int *isec1)
{
  int i, isvn;

  isec1[36] = GET_UINT1(pds[40]); /* extension identifier */
  for ( i = 0; i < 11; i++) 
    { 
      isec1[37+i] =  GET_UINT1(pds[41+i]);
    } 

  isvn = GET_UINT2(pds[52],pds[53]);
  
  isec1[48] =  isvn % 0x8000;              /* DWD experiment identifier            */
  isec1[49] =  isvn >> 15;                 /* DWD run type (0=main, 2=ass, 3=test) */
  isec1[50] =  GET_UINT1(pds[54]);         /* User id, specified by table          */
  isec1[51] =  GET_UINT2(pds[55],pds[56]); /* Experiment identifier                */
  isec1[52] =  GET_UINT2(pds[57],pds[58]); /* Ensemble identification by table     */
  isec1[53] =  GET_UINT2(pds[59],pds[60]); /* Number of ensemble members           */
  isec1[54] =  GET_UINT2(pds[61],pds[62]); /* Actual number of ensemble member     */
  isec1[55] =  GET_UINT1(pds[63]);         /* Model major version number           */
  isec1[56] =  GET_UINT1(pds[64]);         /* Model minor version number           */

}

int decodePDS(unsigned char *pds, int *isec0, int *isec1)
{
  int pdsLen;

  pdsLen = PDS_Len;

  ISEC1_CodeTable      = PDS_CodeTable;
  ISEC1_CenterID       = PDS_CenterID;
  ISEC1_ModelID        = PDS_ModelID;
  ISEC1_GridDefinition = PDS_GridDefinition;
  ISEC1_Sec2Or3Flag    = PDS_Sec2Or3Flag;
  ISEC1_Parameter      = PDS_Parameter;
  ISEC1_LevelType      = PDS_LevelType;

  if ( (ISEC1_LevelType !=  20) && 
       (ISEC1_LevelType !=  99) && 
       (ISEC1_LevelType != 100) && 
       (ISEC1_LevelType != 103) && 
       (ISEC1_LevelType != 105) && 
       (ISEC1_LevelType != 107) && 
       (ISEC1_LevelType != 109) && 
       (ISEC1_LevelType != 111) && 
       (ISEC1_LevelType != 113) && 
       (ISEC1_LevelType != 115) && 
       (ISEC1_LevelType != 117) && 
       (ISEC1_LevelType != 125) && 
       (ISEC1_LevelType != 127) && 
       (ISEC1_LevelType != 160) && 
       (ISEC1_LevelType != 210) )
    {
      ISEC1_Level1 = PDS_Level1;
      ISEC1_Level2 = PDS_Level2;
    }
  else
    {
      ISEC1_Level1 = PDS_Level;
      ISEC1_Level2 = 0;
    }

  ISEC1_Year           = PDS_Year;
  ISEC1_Month          = PDS_Month;
  ISEC1_Day            = PDS_Day;
  ISEC1_Hour           = PDS_Hour;
  ISEC1_Minute         = PDS_Minute;
  ISEC1_TimeUnit       = PDS_TimeUnit;
  ISEC1_TimePeriod1    = PDS_TimePeriod1;
  ISEC1_TimePeriod2    = PDS_TimePeriod2;
  ISEC1_TimeRange      = PDS_TimeRange;
  ISEC1_AvgNum         = PDS_AvgNum;
  ISEC1_AvgMiss        = PDS_AvgMiss;

  if ( ISEC0_GRIB_Version == 1 )
    {
      ISEC1_Century        = PDS_Century;
      ISEC1_SubCenterID    = PDS_Subcenter;
      ISEC1_DecScaleFactor = PDS_DecimalScale;
    }
  else
    {
      ISEC1_Century        = 1;
      ISEC1_SubCenterID    = 0;
      ISEC1_DecScaleFactor = 0;
    }

  ISEC1_LocalFLag = 0;
  if ( pdsLen > 28 )
    {
      int localextlen;
      localextlen = pdsLen-28;

      ISEC1_LocalFLag = 1;

      if ( ISEC1_CenterID == 78 ) 
	{
	  if ( pds[40] == 254 ) 
	    {
	      decodePDS_DWD_local_Extension_254(pds, isec1);
	    }
	  else if ( pds[40] == 253 )
	    { 
	      decodePDS_DWD_local_Extension_253(pds, isec1);
	    }
	}
      else   
	{
	  int i;
	  for ( i = 0; i < localextlen; i++ )
	    {
	      isec1[24+i] = pds[28+i];
	    }
	}
    }
  return (pdsLen);
}

int decodeGDS(unsigned char  *gds, int *isec0, int *isec2, double *fsec2, int dfunc)
{
  static char func[] = "decodeGDS";
  int imisng;
  int  ReducedGrid = FALSE, VertCoorTab = FALSE;
  int  locnv = 0, locnl;
  int  jlenl;
  int  i;
  int iexp, imant;
  int ipvpl, ipl;
  int gdsLen = 0;
#if defined (SX) || defined (ES)
  unsigned char *igrib;
  GRIBPACK *lgrib = NULL;
  size_t lGribLen = 0;
#endif

  imisng = 0;

  gdsLen = GDS_Len;

  ipvpl = GDS_PVPL;
  if ( ipvpl == 0 ) ipvpl = 255;

  if ( ipvpl != 255 )
    { /* Either vct or regrid */
      if ( GDS_NV != 0 )
	{ /* we have vct */
	  VertCoorTab = TRUE;
	  ipl =  4*GDS_NV + ipvpl - 1;
	  if ( ipl < gdsLen )
	    {
	      ReducedGrid = TRUE;
	    }
	}
      else
	{
	  VertCoorTab = FALSE;
	  ReducedGrid = TRUE;
	}
      /*	  ReducedGrid = (gdsLen - 32 - 4*GDS_NV); */
    }
 
  if ( ISEC0_GRIB_Version == 0 )
    {
      if ((gdsLen - 32) > 0) VertCoorTab = TRUE;
      else                   VertCoorTab = FALSE;
    }
  
  if ( ReducedGrid )
    {
      ISEC2_Reduced = TRUE;
      locnl = GDS_PVPL - 1 + (VertCoorTab * 4 * GDS_NV);
      jlenl = (gdsLen - locnl)  >> 1;
      for ( i = 0; i < jlenl; i++ )
	{
	  ISEC2_RowLon(i) = GET_UINT2(gds[locnl+2*i], gds[locnl+2*i+1]);
	}
    }

  ISEC2_GridType = GDS_GridType;

  /*
     Gaussian grid definition.
  */
  if ( ISEC2_GridType == GTYPE_LATLON    ||
       ISEC2_GridType == GTYPE_GAUSSIAN  ||
       ISEC2_GridType == GTYPE_LATLON_ROT )
    {
      if ( ! ReducedGrid ) ISEC2_NumLon = GDS_NumLon;
      ISEC2_NumLat    = GDS_NumLat;
      ISEC2_FirstLat  = GDS_FirstLat;
      ISEC2_FirstLon  = GDS_FirstLon;
      ISEC2_ResFlag   = GDS_ResFlag;
      ISEC2_LastLat   = GDS_LastLat;
      ISEC2_LastLon   = GDS_LastLon;
      ISEC2_LonIncr   = GDS_LonIncr;
      ISEC2_NumPar    = GDS_NumPar;
      ISEC2_ScanFlag  = GDS_ScanFlag;
      ISEC2_LatSP     = GDS_LatSP;
      ISEC2_LonSP     = GDS_LonSP;
      FSEC2_RotAngle  = GDS_RotAngle;

      /*
	if ( Lons != Longitudes || Lats != Latitudes )
	Error(func, "Latitude/Longitude Conflict");
      */
    }
  else if ( ISEC2_GridType == GTYPE_GAUSSIAN     ||
	    ISEC2_GridType == GTYPE_GAUSSIAN_ROT ||
	    ISEC2_GridType == GTYPE_GAUSSIAN_STR ||
	    ISEC2_GridType == GTYPE_GAUSSIAN_ROTSTR )
    {
      /*
      iret = decodeGDS_GG(gds, gdspos, isec0, isec2, imisng);
      */
    }
  else if ( ISEC2_GridType == GTYPE_LATLON     ||
	    ISEC2_GridType == GTYPE_LATLON_ROT ||
	    ISEC2_GridType == GTYPE_LATLON_STR ||
	    ISEC2_GridType == GTYPE_LATLON_ROTSTR )
    {
      /*
      iret = decodeGDS_LL(gds, gdspos, isec0, isec2, imisng);
      */
    }
  else if ( ISEC2_GridType == GTYPE_SPECTRAL )
    {
      ISEC2_PentaJ  = GDS_PentaJ; /* Truncation */
      ISEC2_PentaK  = GDS_PentaK;
      ISEC2_PentaM  = GDS_PentaM;
      ISEC2_RepType = GDS_RepType;
      ISEC2_RepMode = GDS_RepMode;
      isec2[ 6] = 0;
      isec2[ 7] = 0;
      isec2[ 8] = 0;
      isec2[ 9] = 0;
      isec2[10] = 0;
      /*
      iret = decodeGDS_SH(gds, gdspos, isec0, isec2, imisng);
      */
    }
  else if ( ISEC2_GridType == GTYPE_TRIANGULAR )
    {
      ISEC2_TRI_NI2    = GDS_TRI_NI2;
      ISEC2_TRI_NI3    = GDS_TRI_NI3;
      ISEC2_TRI_ND     = GDS_TRI_ND;
      ISEC2_TRI_NI     = GDS_TRI_NI;
      ISEC2_TRI_AFlag  = GDS_TRI_AFlag;
      ISEC2_TRI_LatPP  = GDS_TRI_LatPP;
      ISEC2_TRI_LonPP  = GDS_TRI_LonPP;
      ISEC2_TRI_LonMPL = GDS_TRI_LonMPL;
      ISEC2_TRI_BFlag  = GDS_TRI_BFlag;
      /*
      iret = decodeGDS_TR(gds, gdspos, isec0, isec2, imisng);
      */
    }
  else
    Message(func, "Gridtype %d unsupported", ISEC2_GridType);

      /*    vertical coordinate parameters for hybrid levels.     */
      /*    get number of vertical coordinate parameters, if any. */

  ISEC2_NumVCP = 0;

  isec2[17] = 0;
  isec2[18] = 0;

  if ( VertCoorTab == TRUE )
    {
      if ( ISEC0_GRIB_Version  == 0 )
	{
	  locnv = 32;
	  ISEC2_NumVCP = (gdsLen - 32) >> 2;
	}
      else
	{
	  locnv = GDS_PVPL - 1;
	  ISEC2_NumVCP = GDS_NV;
	}
#if defined (SX) || defined (ES)
      lGribLen = 4*ISEC2_NumVCP;	      
      lgrib    = (GRIBPACK *) malloc(lGribLen*sizeof(GRIBPACK));

      igrib = &gds[locnv];
      if ( ISEC2_NumVCP > 0 ) (void) UNPACK_GRIB(igrib, lgrib, lGribLen, -1L);
      for ( i = 0; i < ISEC2_NumVCP; i++ )
	{
	  iexp   = (lgrib[4*i  ]);
	  imant  =((lgrib[4*i+1]) << 16) +
	          ((lgrib[4*i+2]) <<  8) +
	           (lgrib[4*i+3]);
	  fsec2[10+i] = pow(2.0, -24.0) * imant * pow(16.0, (double)(iexp - 64));
	}

      free(lgrib);
#else
      for ( i = 0; i < ISEC2_NumVCP; i++ )
	{
	  iexp   = (gds[locnv+4*i  ]);
	  imant  =((gds[locnv+4*i+1]) << 16) +
	          ((gds[locnv+4*i+2]) <<  8) +
	           (gds[locnv+4*i+3]);
	  fsec2[10+i] = decfp2(iexp,imant);
	}
#endif
    }

  return (gdsLen);
}

int decodeBDS(unsigned char *bds, int *isec2, int *isec4, double *fsec4, int fsec4len, int dfunc, int *iret)
{
  static char func[] = "decodeBDS";
  unsigned char *igrib;
  GRIBPACK *lgrib = NULL;
  int cplx, jup, kup, mup;
  int locnd;
  int jlend, jlenc;
  int i;
  int iflag, irep, lnil, jscale, imiss;
  int ioff;
  double fmin = 0., zscale = 0.;
  int iexp, imant;
  int pcStart, pcScale;
  double *fpdata = fsec4;

  int bdsLen;

  *iret = 0;
  igrib = bds;

  memset(isec4, 0, 42*sizeof(int));

  /* get length of binary data block. */

  bdsLen = BDS_Len;

  /* 4 bit flag / 4 bit count of unused bits at end of block octet. */

  iflag = BDS_Flag;

  /* 0------- grid point           */
  /* 1------- spherical harmonics  */

  irep  = iflag >> 7;

  if ( irep == 0 ) isec4[2] = 0;
  else             isec4[2] = 128;

  /* -0------  simple packing */
  /* -1------ complex packing */

  cplx  = (iflag >> 6) & 1;

  if ( cplx > 0 ) isec4[3] = 64;
  else            isec4[3] =  0;

  /* ----++++ number of unused bits at end of section) */

  lnil  = iflag & 15;
  
  /* scale factor (2 bytes) */;

  jscale = BDS_BinScale;

  /* check for missing data indicators. */

  iexp  = bds[ 6];
  imant = GET_UINT3(bds[ 7], bds[ 8], bds[ 9]);

  imiss = (jscale == 65535 && iexp == 255 && imant == 16777215);

  /* convert reference value and scale factor. */

  if ( imiss == 0 )
    {
      fmin = BDS_RefValue;
      zscale = pow(2.0, (double) jscale);
    }

  /* get number of bits in each data value. */

  ISEC4_NumBits = BDS_NumBits;

  /* octet number of start of packed data */
  /* calculated from start of block 4 - 1 */

  locnd = 11;

  /* if data is in spherical harmonic form, distinguish   */
  /* between simple/complex packing (cplx = 0/1)          */

  if ( irep == 1 && cplx == 0 )
    {
      /*    no unpacked binary data present */

      jup = kup = mup = 0;

      /*    octet number of start of packed data */
      /*    calculated from start of block 4 - 1 */

      locnd = 15;

      /*    get real (0,0) coefficient in grib format and     */
      /*    convert to floating point.                        */

      if ( dfunc != 'J' )
	{
	  if ( imiss ) *fpdata++ = 0.0;
	  else         *fpdata++ = BDS_RealCoef;
	}
    }

  ioff = irep;

  if ( irep == 1 && cplx == 1 )
    {
      /*    scaling factor */
      isec4[16] = BDS_Power;

      /*    pentagonal resolution parameters of the */
      /*    unpacked section of data field          */

      jup = bds[15];
      kup = bds[16];
      mup = bds[17];

      isec4[17] = jup;
      isec4[18] = kup;
      isec4[19] = mup;

      /*    unpacked binary data */

      locnd = 18;

      if ( dfunc != 'J' )
	for (i = 0; i < ((jup+1)*(jup+2)); i++)
	  {
	    iexp   = (bds[locnd+4*i  ]);
	    imant  =((bds[locnd+4*i+1]) << 16) +
	            ((bds[locnd+4*i+2]) <<  8) +
                     (bds[locnd+4*i+3]);

	    if ( imiss ) *fpdata++ = 0.0;
	    else         *fpdata++ = decfp2(iexp,imant);
	  }
      locnd = 18 + 4*(jup+1)*(jup+2);
      ioff = (jup+1)*(jup+2);
    }

  /* decode data values to floating point and store in fsec4. */
  /* first calculate the number of data values.                */
  /* Take into account that spherical harmonics can be packed  */
  /* simple (cplx = 0) or complex (cplx = 1)                   */

  jlend = bdsLen - locnd;
  if ( ISEC4_NumBits == 0 )
    {
      if ( jlend > 1 )
	Error(func, "Number of bits per data value = 0");

      jlend = ISEC2_NumLon*ISEC2_NumLat;
    }
  else
    jlend = (jlend * 8 - lnil) / ISEC4_NumBits;


  ISEC4_NumValues        = jlend + ioff;
  ISEC4_NumNonMissValues = 0;

  if ( dfunc == 'J' ) return (bdsLen);

  /* check length of output array. */
  
  if ( ! (dfunc == 'J') )
    if ( jlend+ioff > fsec4len )
      {
	*iret = 710;
	gprintf(func, " Output array too small. Length = %d", fsec4len);
	gprintf(func, " Number of values = %d", jlend+ioff);
	gprintf(func, " Return code =  %d", *iret);
	return (0);
      }

  if ( imiss ) memset((char *)fpdata, 0, jlend*sizeof(double));
  else
    {
      igrib += locnd;
      jlenc = jlend * ISEC4_NumBits / 8;
      if ( ISEC4_NumBits ==  8 || ISEC4_NumBits == 16 ||
	   ISEC4_NumBits == 24 || ISEC4_NumBits == 32 )
	{
	  if ( jlenc > 0 ) 
	    {
	      lgrib = (GRIBPACK *) malloc(jlenc*sizeof(GRIBPACK));
	      if ( lgrib == NULL ) SysError(func, "No Memory!");

	      (void) UNPACK_GRIB(igrib, lgrib, jlenc, -1L);
	      /*
		for ( i = 0; i < jlenc; ++i ) lgrib[i] = (long) igrib[i];
	      */
	    }
	}

      if ( ISEC4_NumBits ==  0 )
	{
	  for ( i = 0; i < jlend; i++ )
	    fpdata[i] = fmin;
	}
      else if ( ISEC4_NumBits ==  8 )
         for ( i = 0; i < jlend; i++ )
	   {
             fpdata[i] = fmin + zscale * lgrib[i];
	   }
      else if ( ISEC4_NumBits == 16 )
         for ( i = 0; i < jlend; i++ )
	   {
             fpdata[i] = fmin + zscale *
                         ((lgrib[2*i  ] <<  8) +  lgrib[2*i+1]);
	   }
      else if ( ISEC4_NumBits == 24 )
         for ( i = 0; i < jlend; i++ )
	   {
             fpdata[i] = fmin + zscale *
                         ((lgrib[3*i  ] << 16) + (lgrib[3*i+1] <<  8) +
                           lgrib[3*i+2]);
	   }
      else if ( ISEC4_NumBits == 32 )
         for ( i = 0; i < jlend; i++ )
	   {
             fpdata[i] = fmin + zscale *
                         ((lgrib[4*i  ] << 24) + (lgrib[4*i+1] << 16) +
                          (lgrib[4*i+2] <<  8) +  lgrib[4*i+3]);
	   }
      else if ( ISEC4_NumBits <= 25 )
	{
	  unsigned char *bits = igrib;
	  unsigned int jmask;
	  int tbits = 0;
	  int n_bits = ISEC4_NumBits;
	  int t_bits = 0;
	  jmask = (1 << n_bits) - 1;
	  /* code from wgrib routine BDS_unpack */
	  for ( i = 0; i < jlend; i++ )
	    {
	      while ( t_bits < n_bits )
		{
		  tbits = (tbits * 256) + *bits++;
		  t_bits += 8;
                }
	      t_bits -= n_bits;
	      fpdata[i] = (tbits >> t_bits) & jmask;
            }
	  /* at least this vectorizes :) */
	  for ( i = 0; i < jlend; i++ )
	    fpdata[i] = fmin + zscale*fpdata[i];
	}
      else
	{
	  fprintf(stderr," Unimplemented packing factor %d\n", ISEC4_NumBits);
	  exit(1);
	}
    }

  if ( lgrib ) free(lgrib);

  if ( irep == 1 && cplx == 1 )
    {
      pcStart = isec4[19];
      pcScale = isec4[16];
      scatterComplex(fsec4, pcStart, ISEC2_PentaJ, ISEC4_NumValues);
      scaleComplex(fsec4, pcStart, pcScale, ISEC2_PentaJ);
    }

  return (bdsLen);
}

void gribDecode(int *isec0, int *isec1, int *isec2, double *fsec2, int *isec3,
		double *fsec3, int *isec4, double *fsec4, int fsec4len, int *kgrib,
		int kleng, int *kword, int dfunc, int *iret)
{
  static char func[] = "gribDecode";
  UCHAR *is = NULL, *pds = NULL, *gds = NULL, *bms = NULL, *bds = NULL;
  int isLen = 0, pdsLen = 0, gdsLen = 0, bmsLen = 0, bdsLen = 0, esLen = 0;
  int gribLen = 0;
  int gdsIncluded = FALSE;
  int bmsIncluded = FALSE;
  int bitmapSize = 0;
  int ldebug = FALSE;
  int llarge = FALSE, l_iorj = FALSE;
  int lsect2 = FALSE, lsect3 = FALSE;
  int fsec4off;

  *iret = 0;

  grsdef();

  ISEC2_Reduced = FALSE;

  /*
    ----------------------------------------------------------------
    IS Indicator Section (Section 0)
    ----------------------------------------------------------------
  */
  is = (unsigned char *) &kgrib[0];

  isLen = decodeIS(is, isec0, iret);

  /*
    If count is negative, have to rescale by factor of -120.
    This is a fixup to get round the restriction on product lengths
    due to the count being only 24 bits. It is only possible because
    the (default) rounding for GRIB products is 120 bytes.
  */
  if (  ISEC0_GRIB_Len < 0 )
    {
      if ( ldebug )
	gprintf(func, "Special case, negative length multiplied by -120");
      llarge = TRUE;
      ISEC0_GRIB_Len *= (-120);
    }
  /*
    When decoding or calculating length, previous editions
    of the GRIB code must be taken into account.

    In the table below, covering sections 0 and 1 of the GRIB
    code, octet numbering is from the beginning of the GRIB
    message;
    * indicates that the value is not available in the code edition;
    R indicates reserved, should be set to 0;
    Experimental edition is considered as edition -1.

    GRIB code edition -1 has fixed length of 20 octets for
    section 1, the length not included in the message.
    GRIB code edition 0 has fixed length of 24 octets for
    section 1, the length being included in the message.
    GRIB code edition 1 can have different lengths for section
    1, the minimum being 28 octets, length being included in
    the message.

                                         Octet numbers for code
                                                  editions

                 Contents.                   -1      0      1
                 ---------                ----------------------
       Letters GRIB                          1-4    1-4    1-4
       Total length of GRIB message.          *      *     5-7
       GRIB code edition number               *      *      8
       Length of Section 1.                   *     5-7    9-11
       Reserved octet (R).                    *      8(R)   *
       Version no. of Code Table 2.           *      *     12
       Identification of centre.              5      9     13
       Generating process.                    6     10     14
       Grid definition .                      7     11     15
       Flag (Code Table 1).                   8     12     16
       Indicator of parameter.                9     13     17
       Indicator of type of level.           10     14     18
       Height, pressure etc of levels.      11-12  15-16  19-20
       Year of century.                      13     17     21
       Month.                                14     18     22
       Day.                                  15     19     23
       Hour.                                 16     20     24
       Minute.                               17     21     25
       Indicator of unit of time.            18     22     26
       P1 - Period of time.                  19     23     27
       P2 - Period of time                  20(R)   24     28
       or reserved octet (R).
       Time range indicator.                21(R)   25     29
       or reserved octet (R).
       Number included in average.       22-23(R)  26-27  30-31
       or reserved octet (R).
       Number missing from average.         24(R)  28(R)   32
       or reserved octet (R).
       Century of data.                       *      *     33
       Designates sub-centre if not 0.        *      *     34
       Decimal scale factor.                  *      *    35-36
       Reserved. Set to 0.                    *      *    37-48
       (Need not be present)
       For originating centre use only.       *      *    49-nn
       (Need not be present)

    Identify which GRIB code edition is being decoded.

    In GRIB edition 1, the edition number is in octet 8.
    In GRIB edition 0, octet 8 is reserved and set to 0.
    In GRIB edition -1, octet 8 is a flag field and can have a
    a valid value of 0, 1, 2 or 3.

    However, GRIB edition number 0 has a fixed
    length of 24, included in the message, for section 1, so
    if the value extracted from octets 5-7 is 24 and that from
    octet 8 is 0, it is safe to assume edition 0 of the code.

  */
  if ( ISEC0_GRIB_Len == 24 && ISEC0_GRIB_Version == 0 )
    {
      /*
	Set length of GRIB message to missing data value.
      */
      ISEC0_GRIB_Len = 0;
    }
  /*
    If Grib Edition 1 and only length is required, go to section 9.
  */
  if ( dfunc == 'L' ) goto LABEL900;

  /*
    ----------------------------------------------------------------
    PDS Product Definition Section (Section 1)
    ----------------------------------------------------------------
  */ 
  pds = is + isLen;

  pdsLen = decodePDS(pds, isec0, isec1);

  /*
    ----------------------------------------------------------------
    GDS Grid Description Section (Section 2)
    ----------------------------------------------------------------
  */
  gdsIncluded = ISEC1_Sec2Or3Flag & 128;

  if ( gdsIncluded )
    {
      gds = is + isLen + pdsLen;

      gdsLen = decodeGDS(gds, isec0, isec2, fsec2, dfunc);
    }

  /*
    ----------------------------------------------------------------
    BMS Bit-Map Section Section (Section 3)
    ----------------------------------------------------------------
  */ 
  bmsIncluded = ISEC1_Sec2Or3Flag & 64;

  isec3[0] = 0;
  if ( bmsIncluded )
    {
      bms = is + isLen + pdsLen + gdsLen;

      bmsLen = BMS_Len;
      bitmapSize = ((bmsLen - 6)<<3) - BMS_UnusedBits;
      /*
      fprintf(stderr," bitmapSize = %d\n", bitmapSize);
      */
    }

  /*
    ----------------------------------------------------------------
    BDS Binary Data Section (Section 4)
    ----------------------------------------------------------------
  */
  bds = is + isLen + pdsLen + gdsLen + bmsLen;

  bdsLen = decodeBDS(bds, isec2, isec4, fsec4, fsec4len, dfunc, iret);

  if ( *iret != 0 ) return;

  if ( bitmapSize > 0 )
    {
      ISEC4_NumNonMissValues = ISEC4_NumValues;
      ISEC4_NumValues        = bitmapSize;

      if ( ! (dfunc == 'J') )
	{
	  int i, j;
	  GRIBPACK bitmap;
	  GRIBPACK *imask;

	  fsec4off = bitmapSize - ISEC4_NumNonMissValues;

	  /*
	  unsigned char *bitmap;
	  bitmap = BMS_Bitmap;
	  j = ISEC4_NumNonMissValues;
	  for ( i = ISEC4_NumValues-1; i >= 0; i-- )
	    {
	      if ( (bitmap[i/8]>>(7-(i&7)))&1 )
		fsec4[i] = fsec4[--j];
	      else
		fsec4[i] = FSEC3_MissVal;
	    }
	  */

	  imask = (GRIBPACK *) malloc(bitmapSize*sizeof(GRIBPACK));

	  (void) UNPACK_GRIB(BMS_Bitmap, imask, bitmapSize/8, -1L);

#if defined (CRAY)
#pragma _CRI ivdep
#endif
#if defined (SX) || defined (ES)
#pragma vdir nodep
#endif
#ifdef __uxpch__
#pragma loop novrec
#endif
	  for ( i = bitmapSize/8-1; i >= 0; i-- )
	    {
	      bitmap = imask[i];
	      imask[i*8+0] = 1 & (bitmap >> 7);
	      imask[i*8+1] = 1 & (bitmap >> 6);
	      imask[i*8+2] = 1 & (bitmap >> 5);
	      imask[i*8+3] = 1 & (bitmap >> 4);
	      imask[i*8+4] = 1 & (bitmap >> 3);
	      imask[i*8+5] = 1 & (bitmap >> 2);
	      imask[i*8+6] = 1 & (bitmap >> 1);
	      imask[i*8+7] = 1 & (bitmap);
	    }
	  /*
#if defined (CRAY)
#pragma _CRI ivdep
#endif
#if defined (SX) || defined (ES)
#pragma vdir nodep
#endif
#ifdef __uxpch__
#pragma loop novrec
#endif
	  for ( i = ISEC4_NumNonMissValues-1; i >= 0 ; i-- )
	    fsec4[fsec4off+i] = fsec4[i];
	  */
	  j = ISEC4_NumNonMissValues;

#if defined (CRAY)
#pragma _CRI ivdep
#endif
#if defined (SX) || defined (ES)
#pragma vdir nodep
#endif
#ifdef __uxpch__
#pragma loop novrec
#endif
	  for ( i = ISEC4_NumValues-1; i >= 0; i-- )
	    fsec4[i] = imask[i] ? fsec4[--j] : FSEC3_MissVal;

	  /*
	  for ( i = 0; i < ISEC4_NumValues; i++ )
	    fsec4[i] = imask[i] ? fsec4[fsec4off++] : FSEC3_MissVal;
	  */

	  free(imask);
	}
    }

  if ( dfunc == 'R' && ISEC2_Reduced )
    {
      double *ztemp;
      double rmiss;
      int  nlon, nlat;
      rmiss = -999.9;
      ISEC2_Reduced = 0;
      ISEC2_NumLon = ISEC2_NumLat*2;
      nlon = ISEC2_NumLon;
      nlat = ISEC2_NumLat;
      ISEC4_NumValues = nlon*nlat;

      ztemp = (double *) malloc(ISEC4_NumValues*sizeof(double));
      if ( ztemp == NULL ) SysError(func, "No Memory!");
      (void) qu2reg2(fsec4, ISEC2_RowLonPtr, nlat, nlon, ztemp, FSEC3_MissVal, iret);
      free(ztemp);
    }


  if ( ISEC0_GRIB_Version == 1 ) isLen = 8;
  esLen = 4;

  gribLen = isLen + pdsLen + gdsLen + bmsLen + bdsLen + esLen;

  if ( ISEC0_GRIB_Len )
    if ( gribLen != ISEC0_GRIB_Len )
      {
	Warning(func, "grib1Len = %d gribLen = %d", ISEC0_GRIB_Len, gribLen);
      }

  ISEC0_GRIB_Len = gribLen;

  *kword = gribLen / sizeof(int);
  if ( (size_t) gribLen != *kword * sizeof(int) ) *kword += 1;

  /*
    ----------------------------------------------------------------
    Section 9 . Abort/return to calling routine.
    ----------------------------------------------------------------
  */
 LABEL900:;

  if ( ldebug )
    {
      gprintf(func, "Section 9.");
      gprintf(func, "Output values set -");

      gribPrintSec0(isec0);
      gribPrintSec1(isec0, isec1);
      /*
	Print section 2 if present.
      */
      if ( lsect2 ) gribPrintSec2DP(isec0, isec2, fsec2);

      if ( ! l_iorj )
	{
	  /*
	    Print section 3 if present.
	  */
	  if ( lsect3 ) gribPrintSec3DP(isec0, isec3, fsec3);

	  gribPrintSec4DP(isec0, isec4, fsec4);
	  /*
	    Special print for 2D spectra wave field real values in
	    section 4
	  */
	  if ( (isec1[ 0] ==  140) && 
	       (isec1[ 1] ==   98) && 
	       (isec1[23] ==    1) && 
	       ((isec1[39] == 1045) || (isec1[39] == 1081))  && 
	       ((isec1[ 5] ==  250) || (isec1[ 5] ==  251)) )
	    gribPrintSec4Wave(isec4);
	}
    }
}
#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif




int gribOpen(const char *filename, const char *mode)
{
  int fileID;

  fileID = fileOpen(filename, mode);

#if defined (__sun)
  if ( fileID != FILE_UNDEFID && tolower(*mode) == 'r' )
    {
      fileSetBufferType(fileID, FILE_TYPE_MMAP);
    }
#endif

  return (fileID);  
}

void gribClose(int fileID)
{
  fileClose(fileID);
}

off_t gribGetPos(int fileID)
{
  return (fileGetPos(fileID));
}

int gribCheckFiletype(int fileID)
{
  static char func[] = "gribCheckFiletype";
  int ierr;
  int found = 0;
  char buffer[4];

  if ( fileRead(fileID, buffer, 4) != 4 ) return(found);

  if ( strncmp(buffer, "GRIB", 4) == 0 )
    {
      found = 1;
      if ( GRB_Debug ) Message(func, "found GRIB file = %s", fileInqName(fileID));
    }
  else
    {
      long offset;

      ierr = gribFileSeek(fileID, &offset);
      fileRewind(fileID);
      if ( !ierr )
	{
	  found = 1;
	  if ( GRB_Debug ) Message(func, "found seek GRIB file = %s", fileInqName(fileID));
	}
    }

  return (found);
}

int gribFileSeek(int fileID, long *offset)
{
  /* position file pointer after GRIB */
  static char func[] = "gribFileSeek";
  int buffersize = 4096;
  unsigned char buffer[4096];
  int retry = 4096;
  int i;

  *offset = 0;

  buffer[0] = fileGetc(fileID);
  buffer[1] = fileGetc(fileID);
  buffer[2] = fileGetc(fileID);
  buffer[3] = fileGetc(fileID);
  /*
  fileRead(fileID, buffer, 4);
  */

  while ( retry-- )
    {
      for (i = 0; i < buffersize-4; i++)
	{
	  if (buffer[i  ] == 'G' && 
	      buffer[i+1] == 'R' &&
	      buffer[i+2] == 'I' &&
	      buffer[i+3] == 'B')
	    {
	      if ( GRB_Debug )
		Message(func, "record offset = %d", (int) *offset);
	      return (0);
	    }
	  else
	    {
	      if ( fileEOF(fileID) ) return (-1);
	      buffer[i+4] = fileGetc(fileID);
	      (*offset)++;
	    }
	}
      buffer[0] = buffer[i  ];
      buffer[1] = buffer[i+1];
      buffer[2] = buffer[i+2];
      buffer[3] = buffer[i+3];
    }

  if ( GRB_Debug )
    Message(func, "record offset = %d", (int) *offset);

  return (1);
}

int gribReadSize(int fileID)
{
  static char func[] = "gribReadSize";
  int gribversion, gribsize;
  off_t pos;

  pos = fileGetPos(fileID); 

  gribsize = (fileGetc(fileID) << 16) + (fileGetc(fileID) << 8) + fileGetc(fileID);

  gribversion = fileGetc(fileID);

  if ( gribsize == 24 )
    {
      if ( gribversion != 1 && gribversion != 2 ) gribversion = 0;
    }

  if ( GRB_Debug )
    Message(func, "gribversion = %d", gribversion);

  if ( gribversion == 0 )
    {
      int pdssize = 0, gdssize = 0, bmssize = 0, bdssize = 0;
      int issize = 4, essize = 4;
      int flag;

      pdssize = gribsize;
      fileSetPos(fileID, (off_t) 3, SEEK_CUR);
      if ( GRB_Debug )
	Message(func, "pdssize     = %d", pdssize);
      flag = fileGetc(fileID);
      if ( GRB_Debug )
	Message(func, "flag        = %d", flag);
  
      fileSetPos(fileID, (off_t) pdssize-8, SEEK_CUR);

      if ( flag & 128 )
	{
	  gdssize = (fileGetc(fileID) << 16) + (fileGetc(fileID) << 8) + fileGetc(fileID);
	  fileSetPos(fileID, (off_t) gdssize-3, SEEK_CUR);
	}
      if ( GRB_Debug )
	Message(func, "gdssize     = %d", gdssize);

      if ( flag & 64 )
	{
	  bmssize = (fileGetc(fileID) << 16) + (fileGetc(fileID) << 8) + fileGetc(fileID);
	  fileSetPos(fileID, (off_t) bmssize-3, SEEK_CUR);
	}
      if ( GRB_Debug )
	Message(func, "bmssize     = %d", bmssize);

      bdssize = (fileGetc(fileID) << 16) + (fileGetc(fileID) << 8) + fileGetc(fileID);

      if ( GRB_Debug )
	Message(func, "bdssize     = %d", bdssize);

      gribsize = issize + pdssize + gdssize + bmssize + bdssize + essize;
    }
  else if ( gribversion == 2 )
    {
      int i;
      /* we set gribsize the following way because it doesn't matter then
	 whether int is 4 or 8 bytes long - we don't have to care if the size
	 really fits: if it does not, the record can not be read at all */

      gribsize = 0;
      for ( i = 0; i < 8; i++ ) gribsize = (gribsize << 8) | fileGetc(fileID);
    }
  else if ( gribversion != 1 )
    {
      gribsize = 0;
      Error(func, "GRIB version %d unsupported!", gribversion);
    }

  if ( fileEOF(fileID) ) gribsize = 0;

  if ( GRB_Debug )
    Message(func, "gribsize    = %d", gribsize);

  fileSetPos(fileID, pos, SEEK_SET);

  return (gribsize);
}

int gribGetSize(int fileID)
{
  static char func[] = "gribGetSize";
  size_t recsize;
  long offset;
  int ierr;

  ierr = gribFileSeek(fileID, &offset); /* position file pointer after GRIB */
  if ( ierr > 0 ) Error(func, "GRIB record not found\n");
  if ( ierr == -1 )
    return (0);
  else if ( ierr == 1 )
    return (0);

  recsize = gribReadSize(fileID);

  if ( GRB_Debug )
    Message(func, "recsize = %d\n", recsize);

  fileSetPos(fileID, (off_t) -4, SEEK_CUR);

  return (recsize);
}

int gribRead(int fileID, unsigned char *buffer, size_t *buffersize)
{
  static char func[] = "gribRead";
  long offset;
  int ierr = 0;
  size_t nread, recsize, recsize0;

  ierr = gribFileSeek(fileID, &offset); /* position file pointer after GRIB */
  if ( ierr > 0 ) Error(func, "GRIB record not found\n");
  if ( ierr == -1 )
    {
      *buffersize = 0;
      return (-1);
    }
  else if ( ierr == 1 )
    {
      *buffersize = 0;
      return (-2);
    }

  recsize = gribReadSize(fileID);

  buffer[0] = 'G';
  buffer[1] = 'R';
  buffer[2] = 'I';
  buffer[3] = 'B';

  recsize0 = recsize;

  if ( recsize > *buffersize )
    {
      recsize = *buffersize;
      ierr = -3;
    }

  *buffersize = recsize0;

  nread = fileRead(fileID, &buffer[4], recsize-4);

  if ( nread != recsize-4 ) ierr = 1;

  return (ierr);
}

int gribWrite(int fileID, unsigned char *buffer, size_t buffersize)
{
  static char func[] = "gribWrite";
  int  nwrite = 0;

  if( (nwrite = fileWrite(fileID, buffer, buffersize)) != (int) buffersize )
    {
      perror(func);
      nwrite = -1;
    }

  return ((int) nwrite);
}


/* ============ */
/* scaleComplex */
/* ============ */

void scaleComplex(double *fpdata, int pcStart, int pcScale, int truncation)
{
  static char func[] = "scaleComplex";
  double power;
  double *scale = (double *) malloc((truncation+1)*sizeof(double));
  int  n, m;
  int  index;

  if ( scale == NULL ) SysError(func, "No Memory!");

  if ( pcScale < -10000 || pcScale > 10000 )
    {
      fprintf(stderr, " scaleComplex: Invalid power given %6d\n", pcScale);
      return;
   }

  /* Setup scaling factors = n(n+1)^^p for n = 1 to truncation */

  if ( pcScale == 0 ) return;

  power = (double) pcScale / 1000.;
  scale[0] = 1.0;

  for ( n = 1; n <= truncation; n++ )
    {
      if (pcScale != 1000)
         scale[n] = 1.0 / pow((double) (n*(n+1)), power);
      else
         scale[n] = 1.0 /     (double) (n*(n+1));
    }

  /* Scale the values */

  index = 0;

  for (m = 0; m < pcStart;     m++)
  for (n = m; n <= truncation; n++)
    {
      if ( n >= pcStart )
	{
          fpdata[index  ] *= scale[n];
          fpdata[index+1] *= scale[n];
	}
      index += 2;
    }

  for (m = pcStart; m <= truncation; m++)
  for (n = m;       n <= truncation; n++)
    {
      fpdata[index  ] *= scale[n];
      fpdata[index+1] *= scale[n];
      index += 2;
    }

  free(scale);
}

/* ============== */
/* ScatterComplex */
/* ============== */

void scatterComplex(double *fpdata, int pcStart, int truncation, int dimSP)
{
  static char func[] = "scatterComplex";
  double *fphelp = (double *) malloc(dimSP*sizeof(double));
  int  m, n;
  int  index, inext;

  if ( fphelp == NULL ) SysError(func, "No Memory!");

  index = inext = 0;

  for (m = 0; m <= pcStart;    m++)
  for (n = m; n <= truncation; n++)
    {
      if ( pcStart >= n )
	{
          fphelp[index  ] = fpdata[inext++];
          fphelp[index+1] = fpdata[inext++];
	}
      index += 2;
    }
  index = 0;
  for (m = 0; m <= truncation; m++)
  for (n = m; n <= truncation; n++)
    {
      if ( n > pcStart )
	{
	  fphelp[index  ] = fpdata[inext++];
	  fphelp[index+1] = fpdata[inext++];
	}
      index += 2;
    }
  for (m = 0; m < dimSP; m++) fpdata[m] = fphelp[m];

  free(fphelp);
}

void scm0(double *pdl, double *pdr, double *pfl, double *pfr, int klg)
{
  /* System generated locals */
  double r_1;

  /* Local variables */
  double zfac, zeps, zbeta;
  int jl;
  double zalpha;

  /* **** SCM0   - Apply SCM0 limiter to derivative estimates. */
  /* output: */
  /*   pdl   = the limited derivative at the left edge of the interval */
  /*   pdr   = the limited derivative at the right edge of the interval */
  /* inputs */
  /*   pdl   = the original derivative at the left edge */
  /*   pdr   = the original derivative at the right edge */
  /*   pfl   = function value at the left edge of the interval */
  /*   pfr   = function value at the right edge of the interval */
  /*   klg   = number of intervals where the derivatives are limited */

  /*  define constants */

  zeps = 1.0e-12;
  zfac = (1.0 - zeps) * 3.0;

  for ( jl = 0; jl < klg; ++jl )
    {
      if ( (r_1 = pfr[jl] - pfl[jl], fabs(r_1)) > zeps )
	{
	  zalpha = pdl[jl] / (pfr[jl] - pfl[jl]);
	  zbeta  = pdr[jl] / (pfr[jl] - pfl[jl]);
	  if ( zalpha <= 0.0 ) pdl[jl] = 0.0;
	  if ( zbeta  <= 0.0 ) pdr[jl] = 0.0;
	  if ( zalpha > zfac ) pdl[jl] = zfac * (pfr[jl] - pfl[jl]);
	  if ( zbeta  > zfac ) pdr[jl] = zfac * (pfr[jl] - pfl[jl]);
	}
      else
	{
	  pdl[jl] = 0.0;
	  pdr[jl] = 0.0;
	}
    }
} /* scm0 */

int rowina2(double *p, int ko, int ki, double *pw,
	    int kcode, double msval, int *kret)
{
  /* System generated locals */
  int pw_dim1, pw_offset, i_1;

  /* Local variables */
  double zwt1, zrdi, zpos;
  int jl, ip;
  double zdo, zwt;

  /* Parameter adjustments */
  --p;
  pw_dim1 = ko + 3;
  pw_offset = pw_dim1;
  pw -= pw_offset;

  /* **** ROWINA2 - Interpolation of row of values. */
  /*     Input Parameters. */
  /*     ----------------- */
  /*     P      - Row of values to be interpolated. */
  /*              Dimension must be at least KO. */
  /*     KO     - Number of values required. */
  /*     KI     - Number of values in P on input. */
  /*     PW     - Working array. */
  /*              Dimension must be at least (0:KO+2,3). */
  /*     KCODE  - Interpolation required. */
  /*              1 , linear. */
  /*              3 , cubic. */
  /*     PMSVAL - Value used for missing data indicator. */

  /*     Output Parameters. */
  /*     ------------------ */
  /*     P     - Now contains KO values. */
  /*     KRET  - Return code */
  /*             0, OK */
  /*             Non-zero, error */

  /*     Author. */
  /*     ------- */
  /*     J.D.Chambers    ECMWF     22.07.94 */

  /*     ********************************    */
  /*     Section 1.  Linear interpolation .. */
  /*     ********************************    */

  *kret = 0;

  if ( kcode == 1 )
    {
      /*    Move input values to work array */
      for ( jl = 1; jl <= ki; ++jl )
	pw[jl + pw_dim1] = p[jl];

      /*    Arrange wrap-around value in work array */
      pw[ki + 1 + pw_dim1] = p[1];

      /*    Set up constants to be used to figure out weighting for */
      /*    values in interpolation. */
      zrdi = (double) ki;
      zdo = 1.0 / (double) ko;

      /*    Loop through the output points */
      for ( jl = 1; jl <= ko; ++jl )
	{

	  /*    Calculate weight from the start of row */
	  zpos = (jl - 1) * zdo;
	  zwt = zpos * zrdi;

	  /*    Get the current array position(minus 1) from the weight - */
	  /*    note the implicit truncation. */
	  ip = (int) zwt;

	  /*    If the left value is missing, use the right value */
	  if (pw[ip + 1 + pw_dim1] == msval)
	    {
	      p[jl] = pw[ip + 2 + pw_dim1];
	    }
	  /*    If the right value is missing, use the left value */
	  else if (pw[ip + 2 + pw_dim1] == msval)
	    {
	      p[jl] = pw[ip + 1 + pw_dim1];
	    }
	  /*    If neither missing, interpolate ... */
	  else
	    {

	      /*       Adjust the weight to range (0.0 to 1.0) */
	      zwt -= ip;

	      /*       Interpolate using the weighted values on either side */
	      /*       of the output point position */
	      p[jl] = (1.0 - zwt) * pw[ip + 1 + pw_dim1] +
		zwt * pw[ip + 2 + pw_dim1];
	    }
	}

      /*     *******************************    */
      /*     Section 2.  Cubic interpolation .. */
      /*     *******************************    */

    }
  else if ( kcode == 3 )
    {
      i_1 = ki;
      for ( jl = 1; jl <= i_1; ++jl )
	{
          if ( p[jl] == msval )
	    {
	      fprintf(stderr," ROWINA2: ");
	      fprintf(stderr," Cubic interpolation not supported");
	      fprintf(stderr," for fields containing missing data.\n");
	      *kret = 1;
	      goto L900;
	    }
          pw[jl + pw_dim1] = p[jl];
	}
      pw[pw_dim1] = p[ki];
      pw[ki + 1 + pw_dim1] = p[1];
      pw[ki + 2 + pw_dim1] = p[2];
      i_1 = ki;
      for ( jl = 1; jl <= i_1; ++jl )
	{
          pw[jl + (pw_dim1 << 1)] =
	        - pw[jl - 1 + pw_dim1] / 3.0 -
	          pw[jl     + pw_dim1] * 0.5 +
	          pw[jl + 1 + pw_dim1] - pw[jl + 2 + pw_dim1] / 6.0;
          pw[jl + 1 + pw_dim1 * 3] =
                  pw[jl - 1 + pw_dim1] / 6.0 -
                  pw[jl     + pw_dim1] +
                  pw[jl + 1 + pw_dim1] * 0.5 +
                  pw[jl + 2 + pw_dim1] / 3.0;
	}

      scm0(&pw[(pw_dim1 << 1) + 1], &pw[pw_dim1 * 3 + 2],
	   &pw[pw_dim1 + 1], &pw[pw_dim1 + 2], ki);

      zrdi = (double) ki;
      zdo = 1.0 / (double) ko;
      for ( jl = 1; jl <= ko; ++jl )
	{
          zpos = (jl - 1) * zdo;
          zwt = zpos * zrdi;
          ip = (int) zwt + 1;
          zwt = zwt + 1.0 - ip;
          zwt1 = 1.0 - zwt;
          p[jl] = ((3.0 - zwt1 * 2.0) * pw[ip + pw_dim1] +
                  zwt * pw[ip + (pw_dim1 << 1)]) * zwt1 * zwt1 +
                  ((3.0 - zwt * 2.0) * pw[ip + 1 + pw_dim1] -
                  zwt1 * pw[ip + 1 + pw_dim1 * 3]) * zwt * zwt;
	}

    }
  else
    {
      /*    **************************************    */
      /*    Section 3.  Invalid interpolation code .. */
      /*    **************************************    */
      fprintf(stderr," ROWINA2:");
      fprintf(stderr," Invalid interpolation code = %2d\n",kcode);
      *kret = 2;
    }

L900:
    return 0;
} /* rowina2 */

int qu2reg2(double *pfield, int *kpoint, int klat, int klon,
	    double *ztemp, double msval, int *kret)
{
   static char func[] = "qu2reg2";

   /* System generated locals */
   int i_1, i_2;
   int kcode = 1;

   /* Local variables */
   int ilii, ilio, icode;
   double *zline = NULL;
   double zwork[1929];   /* was [643][3] */
   int iregno, iquano, j210, j220, j230, j240, j225;


   zline = (double *) malloc(klon*sizeof(double));
   if ( zline == NULL ) SysError(func, "No Memory!");

   /* Parameter adjustments */
   --pfield;
   --kpoint;

/* **** QU2REG - Convert quasi-regular grid data to regular. */
/*     Input Parameters. */
/*     ----------------- */
/*     PFIELD     - Array containing quasi-regular grid */
/*                  data. */
/*     KPOINT     - Array containing list of the number of */
/*                  points on each latitude (or longitude) of */
/*                  the quasi-regular grid. */
/*     KLAT       - Number of latitude lines */
/*     KLON       - Number of longitude lines */
/*     KCODE      - Interpolation required. */
/*                  1 , linear - data quasi-regular on */
/*                               latitude lines. */
/*                  3 , cubic -  data quasi-regular on */
/*                               latitude lines. */
/*                  11, linear - data quasi-regular on */
/*                               longitude lines. */
/*                  13, cubic -  data quasi-regular on */
/*                               longitude lines. */
/*     PMSVAL     - Value used for missing data indicator. */
/*     Output Parameters. */
/*     ------------------ */
/*     KRET       - return code */
/*                  0 = OK */
/*                  non-zero indicates fatal error */
/*     PFIELD     - Array containing regular grid data. */
/*     Author. */
/*     ------- */
/*     J.D.Chambers     ECMWF      22.07.94 */
/*     J.D.Chambers     ECMWF      13.09.94 */
/*     Add return code KRET and remove calls to ABORT. */


/* ------------------------------ */
/* Section 1. Set initial values. */
/* ------------------------------ */

   *kret = 0;

/* Check input parameters. */

   if (kcode != 1 && kcode != 3 && kcode != 11 && kcode != 13) {
      fprintf(stderr," QU2REG :");
      fprintf(stderr," Invalid interpolation type code = %2d\n",kcode);
      *kret = 1;
      goto L900;
   }

/* Set array indices to 0. */

   ilii = 0;
   ilio = 0;

/* Establish values of loop parameters. */

   if (kcode > 10) {

/*    Quasi-regular along longitude lines. */

      iquano = klon;
      iregno = klat;
      icode = kcode - 10;
   } else {

/*    Quasi-regular along latitude lines. */

      iquano = klat;
      iregno = klon;
      icode = kcode;
   }

/*     -------------------------------------------------------- */
/**    Section 2. Interpolate field from quasi to regular grid. */
/*     -------------------------------------------------------- */

   i_1 = iquano;
   for (j230 = 1; j230 <= i_1; ++j230) {

      if (iregno != kpoint[j230]) {

/*       Line contains less values than required,so */
/*       extract quasi-regular grid values for a line */

         i_2 = kpoint[j230];
         for (j210 = 1; j210 <= i_2; ++j210) {
            ++ilii;
            zline[j210 - 1] = pfield[ilii];
         }

/*       and interpolate this line. */

         rowina2(zline, iregno, kpoint[j230], zwork, icode, msval, kret);
         if (*kret != 0) goto L900;

/*       Add regular grid values for this line to the
         temporary array. */

         i_2 = iregno;
         for (j220 = 1; j220 <= i_2; ++j220) {
            ++ilio;
            ztemp[ilio - 1] = zline[j220 - 1];
         }

      } else {

/*       Line contains the required number of values, so add */
/*       this line to the temporary array. */

         i_2 = iregno;
         for (j225 = 1; j225 <= i_2; ++j225) {
            ++ilio;
            ++ilii;
            ztemp[ilio - 1] = pfield[ilii];
         }
      }
   }

/* Copy temporary array to user array. */

   i_1 = klon * klat;
   for (j240 = 1; j240 <= i_1; ++j240) {
      pfield[j240] = ztemp[j240 - 1];
   }

/* -------------------------------------------------------- */
/* Section 9. Return to calling routine. Format statements. */
/* -------------------------------------------------------- */

L900:

   free(zline);

   return 0;
} /* qu2reg2_ */


FILE *grprsm;
double fref;
double fmaxval;
int nfref;
int nfmaxval;
int nrnd;
int ndbg;
int nvck;
int nonoff;
int noabort;
int num2ok;
int next2o;
int nloc2o;
int nsubce;


void
grsdef(void)
{
  /*
C---->
C**** GRSDEF - Initial (default) setting of common area variables
C              for GRIBEX package.
C
C     Purpose.
C     --------
C
C     Sets initial values for common area variables for all
C     routines of GRIBEX package, if not already done.
C
C**   Interface.
C     ----------
C
C     CALL GRSDEF
C
C     Input Parameters.
C     -----------------
C
C     None.
C
C     Output Parameters.
C     ------------------
C
C     None.
C
C     Method.
C     -------
C
C     Self-explanatory.
C
C     Externals.
C     ----------
C
C     None.
C
C     Reference.
C     ----------
C
C     See subroutine GRIBEX.
C
C     Comments.
C     ---------
C
C     None
C
C     Author.
C     -------
C
C     J. Clochard, Meteo France, for ECMWF - March 1998.
C
C     Modifications.
C     --------------
C
C     J. Clochard, Meteo France, for ECMWF - June 1999.
C     Add variable NSUBCE.
C     Use a static variable to determine if initialisation has already
C     been done. NUSER removed .
C     Reverse defaults for NEXT2O and NLOC2O, for consistency with
C     version 13.023 of software .
C
  */
  static char func[] = "grsdef";
  /*
C     ----------------------------------------------------------------
C*    Section 0 . Definition of variables.
C     ----------------------------------------------------------------
  */
  char *hndbg, *hnvck;
  char *env_stream;
  static int lfirst = TRUE;

  if ( ! lfirst ) return;

  /*
    ----------------------------------------------------------------
    Section 1 . Set values, conditionally.
    ----------------------------------------------------------------
  */
  /*
    Common area variables have not been set. Set them.
    
    User supplied reference value.
  */
  fref   = 0.0;
  /*
    Reference value supplied by user flag. Set to off.
  */
  nfref  = 0;
  /*
    User supplied maximum value.
  */
  fmaxval   = 0.0;
  /*
    Maximum value supplied by user flag. Set to off.
  */
  nfmaxval  = 0;
  /*
    Set rounding to 120 bytes on.
  */
  nrnd   = 1;
  /*
    Set debug print off.
  */
  ndbg   = 0;
  
  hndbg = getenv("GRIBEX_DEBUG");
  if ( hndbg != NULL )
    {
      if ( !strncmp(hndbg, "ON", 2) )
        ndbg = 1;
      else if( *hndbg == '1')
        ndbg = 1;
      else if( *hndbg == '2')
        ndbg = 2;
      else
        ndbg = 0;
    }
  /*
    Set GRIB value checking on.
  */
  nvck   = 1;
  
  hnvck = getenv("GRIBEX_CHECK");
  if ( hnvck )
    {
      if ( !strncmp(hnvck, "OFF", 3) )
        nvck = 0;
      else
        nvck = 1;
    }
  /*
    See if output stream needs changing
  */
  grprsm = stdout;
  env_stream = getenv("GRPRS_STREAM");
  if ( env_stream )
    {
      if ( isdigit((int) env_stream[0]) )
	{
	  int unit;
	  unit = atoi(env_stream);
	  if ( unit < 1 || unit > 99 )
	    Warning(func, "Invalid number for GRPRS_STREAM: %d\n", unit);
	  else if ( unit == 2 )
	    grprsm = stderr;
	  else if ( unit == 6 )
	    grprsm = stdout;
	  else
	    {
	      char filename[] = "unit.00";
	      sprintf(filename, "%2.2d", unit);
	      grprsm = fopen(filename, "w");
	      if ( ! grprsm )
		SysError(func, "GRPRS_STREAM = %d", unit);
	    }
	}
      else
	{
	  if ( env_stream[0] )
	    {
	      grprsm = fopen(env_stream, "w");
	      if ( ! grprsm )
		SysError(func, "GRPRS_STREAM = %s", env_stream);
	    }
	}
    }
  /*
    Set P factor switch to default, user supplies the P factor.
  */
  nonoff = 0;
  /*
    Set abort flag to NO abort
  */
  noabort = 1;
  /*
    Mark common area values set by user.
  */
  lfirst = FALSE;
  /*
    Exhaustive use of all possible second-order packing methods
    for HOPER='K'. Set to off.
  */
  num2ok  = 0;
  /*
    Use of extended second-order packing methods for grid-point
    encoding (HOPER='C' and 'K'). Set to on.
  */
  next2o  = 1;
  /*
    Use of non-local second-order packing methods for grid-point
    encoding (HOPER='C' and 'K'). Set to on.
  */
  nloc2o  = 1;
  /*
    Use of (all valid) sub-centre values for ECMWF fields encoding .
    encoding. Set to off.
  */
  nsubce  = 0;
}

#undef  IsBigendian
#define IsBigendian()  ( u_byteorder.c[sizeof(long) - 1] )

/* pack 8-bit bytes from 64-bit words to a packed buffer */
/* same as : for ( int i = 0; i < bc; ++i ) cp[i] = (char) up[i]; */

long packInt64(INT64 *up, char *cp, long bc, long tc)
{
#if defined (CRAY)
  (void) _pack(up, cp, bc, tc);
#else
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  char *cp0;
  INT64 upi, *up0, *ip0, *ip1, *ip2, *ip3, *ip4, *ip5, *ip6, *ip7;
  long head, trail, inner, i, j;
  int ipack = sizeof(INT64);
  
  /* Bytes until first word boundary in destination buffer */

  head = ( (long) cp ) & (ipack-1);
  if ( head != 0 ) head = ipack - head;

  inner = bc - head;

  /* Trailing bytes which do not make a full word */

  trail = inner & (ipack-1);

  /* Number of bytes/words to be processed in fast loop */

  inner -= trail;
  inner /= ipack;

  ip0 = up + head;
  ip1 = ip0 + 1;
  ip2 = ip0 + 2;
  ip3 = ip0 + 3;
  ip4 = ip0 + 4;
  ip5 = ip0 + 5;
  ip6 = ip0 + 6;
  ip7 = ip0 + 7;

  up0 = (INT64 *) (cp + head);

  /* Here we should process any bytes until the first word boundary 
   * of our destination buffer 
   * That code is missing so far  because our output buffer is 
   * word aligned by FORTRAN 
   */

  j = 0;

  if ( IsBigendian() )
    {
#if defined (CRAY)
#pragma _CRI ivdep
#endif
#if defined (SX) || defined (ES)
#pragma vdir nodep
#endif
#ifdef __uxpch__
#pragma loop novrec
#endif
      for ( i = 0 ; i < inner ; i++ )
	{
	  upi =             (   ip0[j]         << 56 ) 
	                 |  ( ( ip1[j] & 255 ) << 48 )
	                 |  ( ( ip2[j] & 255 ) << 40 )
	                 |  ( ( ip3[j] & 255 ) << 32 )
	                 |  ( ( ip4[j] & 255 ) << 24 ) ;
	  up0[i] = upi   |  ( ( ip5[j] & 255 ) << 16 )
	                 |  ( ( ip6[j] & 255 ) <<  8 )
	                 |    ( ip7[j] & 255 ) ;
	  j += ipack;
	}
    }
  else
    {
      for ( i = 0 ; i < inner ; i++ )
	{
	  upi =             (   ip7[j]         << 56 ) 
	                 |  ( ( ip6[j] & 255 ) << 48 )
                         |  ( ( ip5[j] & 255 ) << 40 )
                         |  ( ( ip4[j] & 255 ) << 32 )
                         |  ( ( ip3[j] & 255 ) << 24 ) ;
	  up0[i] = upi   |  ( ( ip2[j] & 255 ) << 16 )
                         |  ( ( ip1[j] & 255 ) <<  8 )
                         |    ( ip0[j] & 255 ) ;
	  j += ipack;
	}
    }

  cp0 = (char *) ( up0 + inner );
  if ( trail > 0 )
    {
      up0[inner] = 0;
      for ( i = 0 ; i < trail ; i ++ )
	{
	  *cp0 = (char ) ip0[ipack*inner+i];
	  cp0++;
	}
    }

  if ( tc != -1 )
    {
      bc++;
      *cp0 = (char) tc;
    }
#endif
  return (bc);
}

/* unpack 8-bit bytes from a packed buffer with 64-bit words */
/* same as : for ( int i = 0; i < bc; ++i ) up[i] = (INT64) cp[i]; */

long
unpackInt64(unsigned char *cp, INT64 *up, long bc, long tc)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  unsigned char *cp0;
  unsigned INT64 *up0;
  INT64 *ip0, *ip1, *ip2, *ip3, *ip4, *ip5, *ip6, *ip7;
  long head, trail, inner, i, j;
  long offset;
  int ipack = sizeof(INT64);

  /* Bytes until first word boundary in source buffer */

  head = ( (long) cp ) & (ipack-1);
  if ( head != 0 ) head = ipack - head;
  if ( head > bc ) head = bc;

  inner = bc - head;

  /* Trailing bytes which do not make a full word */
 
  trail = inner & (ipack-1);
 
  /* Number of bytes/words to be processed in fast loop */

  inner -= trail;
  inner /= ipack;

  ip0 = up + head;
  ip1 = ip0 + 1;
  ip2 = ip0 + 2;
  ip3 = ip0 + 3;
  ip4 = ip0 + 4;
  ip5 = ip0 + 5;
  ip6 = ip0 + 6;
  ip7 = ip0 + 7;

  up0 = (unsigned INT64 *) (cp + head);

  /* Process any bytes until the first word boundary 
   * of our source buffer 
   */
  for ( i = 0 ; i < head ; i++ ) up[i] = (INT64) cp[i];

  j = 0;

  if ( IsBigendian() )
    {
#if defined (CRAY)
#pragma _CRI ivdep
#endif
#if defined (SX) || defined (ES)
#pragma vdir nodep
#endif
#ifdef __uxpch__
#pragma loop novrec
#endif
      for ( i = 0 ; i < inner ; i++ )
	{
	  ip0[j] = (up0[i] >> 56) & 255;
	  ip1[j] = (up0[i] >> 48) & 255;
	  ip2[j] = (up0[i] >> 40) & 255;
	  ip3[j] = (up0[i] >> 32) & 255;
	  ip4[j] = (up0[i] >> 24) & 255;
	  ip5[j] = (up0[i] >> 16) & 255;
	  ip6[j] = (up0[i] >>  8) & 255;
	  ip7[j] = (up0[i])       & 255;

	  j += ipack;
	}
    }
  else
    {
      for ( i = 0 ; i < inner ; i++ )
	{
	  ip7[j] = (up0[i] >> 56) & 255;
	  ip6[j] = (up0[i] >> 48) & 255;
	  ip5[j] = (up0[i] >> 40) & 255;
	  ip4[j] = (up0[i] >> 32) & 255;
	  ip3[j] = (up0[i] >> 24) & 255;
	  ip2[j] = (up0[i] >> 16) & 255;
	  ip1[j] = (up0[i] >>  8) & 255;
	  ip0[j] = (up0[i])       & 255;

	  j += ipack;
	}
    }

  if ( trail > 0 )
    {
      offset = head + ipack*inner;
      cp0 = cp + offset;
      for ( i = 0 ; i < trail ; i++ ) up[i+offset] = (INT64) cp0[i];
    }
  /*
  if ( tc != -1 ) {
    bc++;
    *cp0 = (char) tc;
  }
  */
  return (bc);
}

/* pack 8-bit bytes from 32-bit words to a packed buffer */
/* same as : for ( int i = 0; i < bc; ++i ) cp[i] = (char) up[i]; */

#if  defined  (INT32)
long packInt32(INT32 *up, char *cp, long bc, long tc)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  char *cp0;
  INT32 *up0, *ip0, *ip1, *ip2, *ip3;
  long head, trail, inner, i, j;
  int ipack = sizeof(INT32);
  
  /* Bytes until first word boundary in destination buffer */

  head = ( (long) cp ) & (ipack-1);
  if ( head != 0 ) head = ipack - head;

  inner = bc - head;

  /* Trailing bytes which do not make a full word */

  trail = inner & (ipack-1);

  /* Number of bytes/words to be processed in fast loop */

  inner -= trail;
  inner /= ipack;

  ip0 = up + head;
  ip1 = ip0 + 1;
  ip2 = ip0 + 2;
  ip3 = ip0 + 3;

  up0 = (INT32 *) (cp + head);

  /* Here we should process any bytes until the first word boundary 
   * of our destination buffer 
   * That code is missing so far  because our output buffer is 
   * word aligned by FORTRAN 
   */

  j = 0;

  if ( IsBigendian() )
    {
#if defined (CRAY)
#pragma _CRI ivdep
#endif
#if defined (SX) || defined (ES)
#pragma vdir nodep
#endif
#ifdef __uxpch__
#pragma loop novrec
#endif
      for ( i = 0 ; i < inner ; i++ )
	{
	  up0[i] =          (   ip0[j]         << 24 ) 
	                 |  ( ( ip1[j] & 255 ) << 16 )
	                 |  ( ( ip2[j] & 255 ) <<  8 )
	                 |    ( ip3[j] & 255 ) ;
	  j += ipack;
	}
    }
  else
    {
      for ( i = 0 ; i < inner ; i++ )
	{
	  up0[i] =          (   ip3[j]         << 24 ) 
	                 |  ( ( ip2[j] & 255 ) << 16 )
                         |  ( ( ip1[j] & 255 ) <<  8 )
                         |    ( ip0[j] & 255 ) ;
	  j += ipack;
	}
    }

  cp0 = (char *) ( up0 + inner );
  if ( trail > 0 )
    {
      up0[inner] = 0;
      for ( i = 0 ; i < trail ; i ++ )
	{
	  *cp0 = (char ) ip0[ipack*inner+i];
	  cp0++;
	}
    }

  if ( tc != -1 )
    {
      bc++;
      *cp0 = (char) tc;
    }

  return (bc);
}
#endif

/* unpack 8-bit bytes from a packed buffer with 32-bit words */
/* same as : for ( int i = 0; i < bc; ++i ) up[i] = (INT32) cp[i]; */

#if  defined  (INT32)
long unpackInt32(unsigned char *cp, INT32 *up, long bc, long tc)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  unsigned char *cp0;
  unsigned INT32 *up0;
  INT32 *ip0, *ip1, *ip2, *ip3;
  long head, trail, inner, i, j;
  long offset;
  int ipack = sizeof(INT32);

  /* Bytes until first word boundary in source buffer */

  head = ( (long) cp ) & (ipack-1);
  if ( head != 0 ) head = ipack - head;
  if ( head > bc ) head = bc;

  inner = bc - head;

  /* Trailing bytes which do not make a full word */
 
  trail = inner & (ipack-1);
 
  /* Number of bytes/words to be processed in fast loop */

  inner -= trail;
  inner /= ipack;

  ip0 = up + head;
  ip1 = ip0 + 1;
  ip2 = ip0 + 2;
  ip3 = ip0 + 3;

  up0 = (unsigned INT32 *) (cp + head);

  /* Process any bytes until the first word boundary 
   * of our source buffer 
   */
  for ( i = 0 ; i < head ; i++ ) up[i] = (INT32) cp[i];

  j = 0;

  if ( IsBigendian() )
    {
#if defined (CRAY)
#pragma _CRI ivdep
#endif
#if defined (SX) || defined (ES)
#pragma vdir nodep
#endif
#ifdef __uxpch__
#pragma loop novrec
#endif
      for ( i = 0 ; i < inner ; i++ )
	{
	  ip0[j] = (up0[i] >> 24) & 255;
	  ip1[j] = (up0[i] >> 16) & 255;
	  ip2[j] = (up0[i] >>  8) & 255;
	  ip3[j] = (up0[i])       & 255;

	  j += ipack;
	}
    }
  else
    {
      for ( i = 0 ; i < inner ; i++ )
	{
	  ip3[j] = (up0[i] >> 24) & 255;
	  ip2[j] = (up0[i] >> 16) & 255;
	  ip1[j] = (up0[i] >>  8) & 255;
	  ip0[j] = (up0[i])       & 255;

	  j += ipack;
	}
    }

  if ( trail > 0 )
    {
      offset = head + ipack*inner;
      cp0 = cp + offset;
      for ( i = 0 ; i < trail ; i++ ) up[i+offset] = (INT32) cp0[i];
    }
  /*
  if ( tc != -1 ) {
    bc++;
    *cp0 = (char) tc;
  }
  */
  return (bc);
}
#endif

void prtbin(int kin, int knbit, int *kout, int *kerr)
{
  /*

    Produces a decimal number with ones and zeroes
    corresponding to the ones and zeroes of the input
    binary number.
    eg input number 1011 binary, output number 1011 decimal.


    Input Parameters:
    
       kin   - Integer variable containing binary number.

       knbit - Number of bits in binary number.

    Output Parameters:

       kout  - Integer variable containing decimal value
               with ones and zeroes corresponding to those of
	       the input binary number.

       kerr  - 0, If no error.
               1, Number of bits in binary number exceeds
	          maximum allowed or is less than 1.


    Converted from EMOS routine PRTBIN.

       Uwe Schulzweida   MPIfM   01/04/2001

  */
  int idec;
  int ik;
  int itemp;
  int j;

  /*
    Check length of binary number to ensure decimal number
    generated will fit in the computer word - in this case will
    it fit in a Cray 48 bit integer?
  */
  if ( knbit < 1 || knbit > 14 )
    {
      *kerr = 1;
      printf(" prtbin : Error in binary number length - %3d bits.\n", knbit);
      return;
    }
  else
    *kerr = 0;
  /*
    -----------------------------------------------------------------
    Section 1. Generate required number.
    -----------------------------------------------------------------
  */
  *kout = 0;
  ik    = kin;
  idec  = 1;

  for ( j = 0; j < knbit; j++ )
    {
      itemp = ik - ( (ik/2)*2 );
      *kout = (*kout) + itemp * idec;
      ik    = ik / 2;
      idec  = idec * 10;
    }

  return;
}

void
ref2ibm(double *pref, int kbits)
{
  /*

    Purpose:
    --------

    Code and check reference value in IBM format

    Input Parameters:
    -----------------

    pref       - Reference value
    kbits      - Number of bits per computer word.

    Output Parameters:
    ------------------

    pref       - Reference value

    Method:
    -------

    Codes in IBM format, then decides to ensure that reference 
    value used for packing is not different from that stored
    because of packing differences.

    Externals.
    ----------

    confp3    - Encode into IBM floating point format.
    decfp2    - Decode from IBM floating point format.

    Reference:
    ----------

    None.

    Comments:
    --------

    None.

    Author:
    -------

    J.D.Chambers     ECMWF      17:05:94

    Modifications:
    --------------

    Uwe Schulzweida   MPIfM   01/04/2001

    Convert to C from EMOS library version 130

  */

  static char func[] = "ref2ibm";
  static int itrnd;
  static int kexp, kmant;
  static double ztemp, zdumm;
  extern int GRB_Debug;

  /* ----------------------------------------------------------------- */
  /*   Section 1. Convert to and from IBM format.                      */
  /* ----------------------------------------------------------------- */

  /*  Convert floating point reference value to IBM representation. */

  itrnd = 1;
  zdumm = ztemp = *pref;
  confp3(zdumm, &kexp, &kmant, kbits, itrnd);

  if ( kexp == 0 && kmant == 0 ) return;

  /*  Set reference value to that actually stored in the GRIB code. */

  *pref = decfp2(kexp, kmant);

  /*  If the nearest number which can be represented in */
  /*  GRIB format is greater than the reference value,  */
  /*  find the nearest number in GRIB format lower      */
  /*  than the reference value.                         */

  if ( ztemp < *pref )
    {
      /*  Convert floating point to GRIB representation */
      /*  using truncation to ensure that the converted */
      /*  number is smaller than the original one.      */

      itrnd = 0;
      zdumm = *pref = ztemp;
      confp3(zdumm, &kexp, &kmant, kbits, itrnd);

      /*  Set reference value to that stored in the GRIB code. */

      *pref = decfp2(kexp, kmant);

      if ( ztemp < *pref )
	{
	  if ( GRB_Debug )
	    {
	      Message(func, "Reference value error.");
	      Message(func, "Notify Met.Applications Section.");
	      Message(func, "ZTEMP = ", ztemp);
	      Message(func, "PREF = ", pref);
	    }
	  *pref = ztemp;
	}
    }

  return;
} /* ref2ibm */
#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif



#undef 	 ABS
#define	 ABS(a)		((a) < 0 ? -(a) : (a))

/*
 * The following two functions convert between Julian day number and
 * Gregorian/Julian dates (Julian dates are used prior to October 15,
 * 1582; Gregorian dates are used after that).  Julian day number 0 is
 * midday, January 1, 4713 BCE.  The Gregorian calendar was adopted
 * midday, October 15, 1582.
 *
 * Author: Robert Iles, March 1994
 *
 * C Porter: Steve Emmerson, October 1995
 *
 * Original: http://www.nag.co.uk:70/nagware/Examples/calendar.f90
 *
 * There is no warranty on this code.
 */


/*
 * Convert a Julian day number to a Gregorian/Julian date.
 */
static void julday_to_gregdate(
    unsigned long	julday,		/* Julian day number to convert */
    int			*year,		/* Gregorian year (out) */
    int			*month,		/* Gregorian month (1-12) (out) */
    int			*day		/* Gregorian day (1-31) (out) */
    )
{
    long	ja, jb;
    int		jc;
    long	jd;
    int		je, iday, imonth, iyear;
    double	xc;

    if (julday < 2299161)
	ja = julday;
    else
    {
	int	ia = (int) (((julday - 1867216) - 0.25) / 36524.25);

	ja = julday + 1 + ia - (int)(0.25 * ia);
    }

    jb = ja + 1524;
    xc = ((jb - 2439870) - 122.1) / 365.25;
    jc = (long) (6680.0 + xc);
    jd = 365 * jc + (int)(0.25 * jc);
    je = (int) ((jb - jd) / 30.6001);

    iday = (int) (jb - jd - (int)(30.6001 * je));

    imonth = je - 1;
    if (imonth > 12)
	imonth -= 12;

    iyear = jc - 4715;
    if (imonth > 2)
	iyear -= 1;
    if (iyear <= 0)
	iyear -= 1;

    *year = iyear;
    *month = imonth;
    *day = iday;
}


/*
 * Convert a Gregorian/Julian date to a Julian day number.
 *
 * The Gregorian calendar was adopted midday, October 15, 1582.
 */
static unsigned long gregdate_to_julday(
    int		year,	/* Gregorian year */
    int		month,	/* Gregorian month (1-12) */
    int		day	/* Gregorian day (1-31) */
    )
{
    long		igreg = 15 + 31 * (10 + (12 * 1582));
    int			iy;	/* signed, origin 0 year */
    int			ja;	/* Julian century */
    int			jm;	/* Julian month */
    int			jy;	/* Julian year */
    unsigned long	julday;	/* returned Julian day number */

    /*
     * Because there is no 0 BC or 0 AD, assume the user wants the start of 
     * the common era if they specify year 0.
     */
    if (year == 0)
	year = 1;

    iy = year;
    if (year < 0)
	iy++;
    if (month > 2)
    {
	jy = iy;
	jm = month + 1;
    }
    else
    {
	jy = iy - 1;
	jm = month + 13;
    }

    /*
     *  Note: SLIGHTLY STRANGE CONSTRUCTIONS REQUIRED TO AVOID PROBLEMS WITH
     *        OPTIMISATION OR GENERAL ERRORS UNDER VMS!
     */
    julday = day + (int)(30.6001 * jm);
    if (jy >= 0)
    {
	julday += 365 * jy;
	julday += jy / 4;
    }
    else
    {
	double		xi = 365.25 * jy;

	if ((int)xi != xi)
	    xi -= 1;
	julday += (int)xi;
    }
    julday += 1720995;

    if (day + (31* (month + (12 * iy))) >= igreg)
    {
	ja = jy/100;
	julday -= ja;
	julday += 2;
	julday += ja/4;
    }

    return julday;
}

/*
 * Encode a date as a double-precision value.
 */
static double utencdate(int year, int month, int day)
{
    return ((long)gregdate_to_julday(year, month, day) - 
	    (long)gregdate_to_julday(2001, 1, 1)) * 86400.0;
}


/*
 * Encode a time as a double-precision value.
 */
static double utencclock(int hours, int minutes, double seconds)
{
    return (hours*60 + minutes)*60 + seconds;
}



/*
 * Decompose a value into a set of values accounting for uncertainty.
 */
static void decomp(
       double	value,
       double	uncer,		/* >= 0 */
       int	nbasis,
       double	*basis,		/* all values > 0 */
       double	*count
    )
{
    int		i;

    for (i = 0; i < nbasis; i++)
    {
	double	r = fmod(value, basis[i]);	/* remainder */

	/* Adjust remainder to minimum magnitude. */
	if (ABS(2*r) > basis[i])
	    r += r > 0
		    ? -basis[i]
		    :  basis[i];

	if (ABS(r) <= uncer)
	{
	    /* The value equals a basis multiple within the uncertainty. */
	    double	half = value < 0 ? -basis[i]/2 : basis[i]/2;
	    modf((value+half)/basis[i], count+i);
	    break;
	}

	value = basis[i] * modf(value/basis[i], count+i);
    }

    for (i++; i < nbasis; i++)
	count[i] = 0;
}


/*
 * Decode a time from a double-precision value.
 */
static void dectime(
	double	value,
	int	*year,
	int	*month,
	int	*day,
	int	*hour,
	int	*minute,
	float	*second)
{
    long	days;
    int  	hours;
    int 	minutes;
    double	seconds;
    double	uncer;		/* uncertainty of input value */
    typedef union
    {
	double	    vec[7];
	struct
	{
	    double	days;
	    double	hours12;
	    double	hours;
	    double	minutes10;
	    double	minutes;
	    double	seconds10;
	    double	seconds;
	}	    ind;
    } Basis;
    Basis	counts;
    static Basis	basis;

    basis.ind.days      = 86400;
    basis.ind.hours12   = 43200;
    basis.ind.hours     = 3600;
    basis.ind.minutes10 = 600;
    basis.ind.minutes   = 60;
    basis.ind.seconds10 = 10;
    basis.ind.seconds   = 1;

    uncer = ldexp(value < 0 ? -value : value, -DBL_MANT_DIG);

    days = (long) floor(value/86400.0);
    value -= days * 86400.0;		/* make positive excess */

    decomp(value, uncer, sizeof(basis.vec)/sizeof(basis.vec[0]),
	   basis.vec, counts.vec);

    days   += (long) counts.ind.days;
    hours   = (int)counts.ind.hours12 * 12 + (int)counts.ind.hours;
    minutes = (int)counts.ind.minutes10 * 10 + (int)counts.ind.minutes;
    seconds = (int)counts.ind.seconds10 * 10 + counts.ind.seconds;

    *second = seconds;
    *minute = minutes;
    *hour = hours;
    julday_to_gregdate(gregdate_to_julday(2001, 1, 1) + days, year, month, day);
}

/*
 * Convert a Gregorian/Julian date and time into a temporal value.
 *
 */
void InvCalendar(int year, int month, int day, int hour,
		 int minute, double second, sUnit unit, double *value)
{
  *value = (utencdate(year, month, day) + 
	    utencclock(hour, minute, second) - unit.origin) /
            unit.factor;
}

/*
 * Convert a temporal value into UTC Gregorian/Julian date and time.
 *
 */
void Calendar(double value, sUnit unit, int *year, int *month,
	      int *day, int *hour, int *minute, double *second)
{
  float	sec;

  dectime(unit.origin + value*unit.factor, year, month, day, hour,
	  minute, &sec);
  *second = sec;
}
#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif



int _ExitOnError   = 1;	/* If set to 1, exit on error       */
int _Verbose = 1;	/* If set to 1, errors are reported */
int _Debug   = 0;       /* If set to 1, debugging           */


void SysError(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  printf("\n");
   fprintf(stderr, "Error (%s) : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  if ( errno )
    perror("System error message ");
	
  exit(1);
}

void Error(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  printf("\n");
   fprintf(stderr, "Error (%s) : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  if ( _ExitOnError ) exit(1);
}

void Warning(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  if ( _Verbose )
    {
       fprintf(stderr, "Warning (%s) : ", caller);
      vfprintf(stderr, fmt, args);
       fprintf(stderr, "\n");
    }

  va_end(args);
}

void Message(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

   fprintf(stdout, "%-18s : ", caller);
  vfprintf(stdout, fmt, args);
   fprintf(stdout, "\n");

  va_end(args);
}
#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif





#if ! defined (HAVE_MMAP)
#  if defined (__sun) || defined (__i386) || defined (__ia64)
#    define HAVE_MMAP
#  endif
#endif

#if defined (HAVE_MMAP)
#  include <sys/mman.h> /* mmap() is defined in this header */
#endif


#if ! defined   (FALSE)
#  define  FALSE  0
#endif

#if ! defined   (TRUE)
#  define  TRUE   1
#endif

typedef struct
{
  int             used;           /* table entry in use            */
  int             flag;           /* access and error flag         */
  int             eof;            /* end of file flag              */
  int             fd;             /* file descriptor used for read */
  FILE           *fp;             /* FILE pointer used for write   */
  int             mode;           /* file access mode              */
  char           *name;           /* file name                     */
  off_t           size;           /* file size                     */
  off_t           position;       /* file position                 */
  long            access;         /* file access                   */
  off_t           byteTrans;      /*                               */
  size_t          blockSize;      /* file block size               */
  int             bufferType;     /* buffer type ( 1:std 2:mmap )  */
  size_t          bufferSize;     /* file buffer size              */
  size_t          mappedSize;     /* mmap buffer size              */
  char           *buffer;         /* file buffer                   */
  long            bufferNumFill;  /* number of buffer fill         */
  char           *bufferPtr;      /* file buffer pointer           */
  off_t           bufferPos;
  off_t           bufferStart;
  off_t           bufferEnd;
  size_t          bufferCnt;
}
F_I_L_E;

enum F_I_L_E_Flags
  {
    FILE_READ  =  01,
    FILE_WRITE =  02,
    FILE_UNBUF =  04,
    FILE_EOF   = 010,
    FILE_ERROR = 020
  };

static F_I_L_E *fileTable;
static int  fileTableSize = 0;

#undef FILE_CHECKID

static int  FileInit = 0;
static int  FileInfo = 0;


#if ! defined (MIN_BUF_SIZE)
#  define  MIN_BUF_SIZE  131072L
#endif


static size_t FileBufferSizeMin = MIN_BUF_SIZE;
static long   FileBufferSizeEnv = -1;
static int    FileBufferTypeEnv =  0;

static int  FILE_Debug = 0;   /* If set to 1, debugging */


static void fileTablePrint(void);

/*
 * A version string.
 */

#undef   LIBVERSION
#define  LIBVERSION      1.3.2
#define  XSTRING(x)	 #x
#define  STRING(x) 	 XSTRING(x)
const char file_libvers[] = STRING(LIBVERSION) " of "__DATE__" "__TIME__;

/*
  21/05/2004  1.3.2 set min I/O Buffersize to 128k
 */


const char *
fileLibraryVersion(void)
{
  return (file_libvers);
}

#ifndef POSIXIO_DEFAULT_PAGESIZE
#define POSIXIO_DEFAULT_PAGESIZE 4096
#endif

static int
pagesize(void)
{
#if defined (HAVE_MMAP)
  return ((int) getpagesize());
#else
  return ((int) POSIXIO_DEFAULT_PAGESIZE);
#endif
}

void fileDebug(int debug)
{
  static char func[] = "fileDebug";

  FILE_Debug = debug;

  if ( FILE_Debug )
    Message(func, "Debug level %d", debug);
}

static void fileCheckID(const char *func, int fileID)
{
  if ( fileID < 0 || fileID >= fileTableSize )
    Error(func, "fileID %d undefined!", fileID);

  if ( ! fileTable[fileID].used )
    {
      if ( FILE_Debug )
	fileTablePrint();
	  
      Error(func, "fileID %d undefined!", fileID);
    }
}

int fileSetBufferType(int fileID, int type)
{
  static char func[] = "fileSetBufferType";
  int ret = 0;

  switch (type)
    {
    case FILE_TYPE_STD:
    case FILE_TYPE_MMAP:
      fileTable[fileID].bufferType = type;
      break;
    default:
      Error(func, "File type %d not implemented!", type);
    }

#if ! defined (HAVE_MMAP)
  if ( type == FILE_TYPE_MMAP ) ret = 1;
#endif

  return (ret);
}

int fileGetBufferType(int fileID)
{
  return (fileTable[fileID].bufferType);
}

int fileFlush(int fileID)
{
#if defined (FILE_CHECKID)
  static char func[] = "fileFlush";

  fileCheckID(func, fileID);
#endif

  return (fflush(fileTable[fileID].fp));
}

void fileClearerr(int fileID)
{
  if ( fileTable[fileID].mode != 'r' )
    clearerr(fileTable[fileID].fp);
}

int fileEOF(int fileID)
{
  return ((fileTable[fileID].flag & FILE_EOF) != 0 );
}

int fileError(int fileID)
{
  return ((fileTable[fileID].flag & FILE_ERROR) != 0 );
}

void fileRewind(int fileID)
{
  fileSetPos(fileID, (off_t) 0, SEEK_SET);
  fileClearerr(fileID);
}

off_t fileGetPos(int fileID)
{
  static char func[] = "fileGetPos";
  off_t filepos = 0;

#if defined (FILE_CHECKID)
  fileCheckID(func, fileID);
#endif

  if ( fileTable[fileID].mode == 'r' )
    filepos = fileTable[fileID].position;
  else
    filepos = ftell(fileTable[fileID].fp);

  if ( FILE_Debug ) Message(func, "Position %ld", filepos);

  return (filepos);
}

int fileSetPos(int fileID, off_t offset, int whence)
{
  static char func[] = "fileSetPos";
  int status = 0;
  off_t position;

#if defined (FILE_CHECKID)
  fileCheckID(func, fileID);
#endif

  if ( FILE_Debug ) Message(func, "Offset %8ld  Whence %3d", (long) offset, whence);

  switch (whence)
    {
    case SEEK_SET:
      if ( fileTable[fileID].mode == 'r' )
	{
	  position = offset;
	  fileTable[fileID].position = position;
	  if ( position < fileTable[fileID].bufferStart ||
	       position > fileTable[fileID].bufferEnd )
	    {
	      if ( fileTable[fileID].bufferType == FILE_TYPE_STD )
		fileTable[fileID].bufferPos = position;
	      else
		fileTable[fileID].bufferPos = position - position % pagesize();

	      fileTable[fileID].bufferCnt = 0;
	      fileTable[fileID].bufferPtr = NULL;
	    }
	  else
	    {
	      if ( fileTable[fileID].bufferPos != fileTable[fileID].bufferEnd + 1 )
		{
		  if ( FILE_Debug )
		    Message(func, "Reset buffer pos from %ld to %ld",
			    fileTable[fileID].bufferPos, fileTable[fileID].bufferEnd + 1);
			    
		  fileTable[fileID].bufferPos = fileTable[fileID].bufferEnd + 1;
		}
	      fileTable[fileID].bufferCnt = fileTable[fileID].bufferEnd - position + 1;
	      fileTable[fileID].bufferPtr = fileTable[fileID].buffer + position -
		                            fileTable[fileID].bufferStart;
	    }
	}
      else
	{
	  status = fseek(fileTable[fileID].fp, offset, whence);
	}
      break;
    case SEEK_CUR:
      if ( fileTable[fileID].mode == 'r' )
	{
	  fileTable[fileID].position += offset;
	  position = fileTable[fileID].position;
	  if ( position < fileTable[fileID].bufferStart ||
	       position > fileTable[fileID].bufferEnd )
	    {
	      if ( fileTable[fileID].bufferType == FILE_TYPE_STD )
		fileTable[fileID].bufferPos = position;
	      else
		fileTable[fileID].bufferPos = position - position % pagesize();

	      fileTable[fileID].bufferCnt = 0;
	      fileTable[fileID].bufferPtr = NULL;
	    }
	  else
	    {
	      if ( fileTable[fileID].bufferPos != fileTable[fileID].bufferEnd + 1 )
		{
		  if ( FILE_Debug )
		    Message(func, "Reset buffer pos from %ld to %ld",
			    fileTable[fileID].bufferPos, fileTable[fileID].bufferEnd + 1);
			    
		  fileTable[fileID].bufferPos = fileTable[fileID].bufferEnd + 1;
		}
	      fileTable[fileID].bufferCnt -= offset;
	      fileTable[fileID].bufferPtr += offset;
	    }
	}
      else
	{
	  status = fseek(fileTable[fileID].fp, offset, whence);
	}
      break;
    default:
      Error(func, "Whence = %d not implemented", whence);
    }

  if ( fileTable[fileID].position < fileTable[fileID].size )
    if ( (fileTable[fileID].flag & FILE_EOF) != 0 )
      fileTable[fileID].flag -= FILE_EOF;

  return (status);
}

static void fileTablePrint(void)
{
  int fileID;
  int lprintHeader = 1;

  for ( fileID = 0; fileID < fileTableSize; fileID++ )
    {
      if ( fileTable[fileID].used )
	{
	  if ( lprintHeader )
	    {
	      fprintf(stderr, "\nFile table:\n");
	      fprintf(stderr, "+-----+---------+");
	      fprintf(stderr, "----------------------------------------------------+\n");
	      fprintf(stderr, "|  ID |  Mode   |");
	      fprintf(stderr, "  Name                                              |\n");
	      fprintf(stderr, "+-----+---------+");
	      fprintf(stderr, "----------------------------------------------------+\n");
	      fprintf(stderr, "| %3d | ", fileID);
	      lprintHeader = 0;
	    }

	  switch ( fileTable[fileID].mode )
	    {
	    case 'r':
	      fprintf(stderr, "read   ");
	      break;
	    case 'w':
	      fprintf(stderr, "write  ");
	      break;
	    case 'a':
	      fprintf(stderr, "append ");
	      break;
	    default:
	      fprintf(stderr, "unknown");
	    }

          fprintf(stderr, " | %-51s|\n", fileTable[fileID].name);
	}
    }

  if ( lprintHeader == 0 )
    {
      fprintf(stderr, "+-----+---------+");
      fprintf(stderr, "----------------------------------------------------+\n");
    }
}

char *fileInqName(int fileID)
{
  return (fileTable[fileID].name);
}

int fileInqMode(int fileID)
{
#if defined (FILE_CHECKID)
  static char func[] = "fileInqMode";

  fileCheckID(func, fileID);
#endif

  return (fileTable[fileID].mode);
}

static void fileTableInitEntry(int fileID)
{
  static char func[] = "fileTableInitEntry";

  if ( fileID < 0 || fileID >= fileTableSize )
    Error(func, "fileID %d undefined!", fileID);


  fileTable[fileID].used          = FALSE;
  fileTable[fileID].flag          = 0;
  fileTable[fileID].fd            = -1;
  fileTable[fileID].fp            = NULL;
  fileTable[fileID].mode          = 0;
  fileTable[fileID].size          = 0;
  fileTable[fileID].name          = NULL;
  fileTable[fileID].access        = 0;
  fileTable[fileID].position      = 0;
  fileTable[fileID].byteTrans     = 0;
  fileTable[fileID].bufferType    = 0;
  fileTable[fileID].bufferSize    = 0;
  fileTable[fileID].mappedSize    = 0;
  fileTable[fileID].buffer        = NULL;
  fileTable[fileID].bufferNumFill = 0;
  fileTable[fileID].bufferStart   = 0;
  fileTable[fileID].bufferEnd     = -1;
  fileTable[fileID].bufferPos     = 0;
  fileTable[fileID].bufferCnt     = 0;
  fileTable[fileID].bufferPtr     = NULL;
}

static int fileTableNewEntry(void)
{
  static char func[] = "fileTableNewEntry";
  int fileID = 0;

  /*
    Look for a free slot in fileTable.
    (Create the table the first time through).
  */
  if ( !fileTableSize )
    {
      int i;

      fileTableSize = 2;
      fileTable = (F_I_L_E *) malloc(fileTableSize*sizeof(F_I_L_E));
      if( fileTable == NULL )
	SysError(func, "Allocation of F_I_L_E failed!");

      for( i = 0; i < fileTableSize; i++ )
	fileTableInitEntry(i);
    }
  else
    {
      while( fileID < fileTableSize )
	{
	  if ( ! fileTable[fileID].used ) break;
	  fileID++;
	}
    }
  /*
    If the table overflows, double its size.
  */
  if ( fileID == fileTableSize )
    {
      int i;

      fileTableSize = 2*fileTableSize;
      fileTable = (F_I_L_E *) realloc(fileTable, fileTableSize*sizeof(F_I_L_E));
      if( fileTable == NULL )
	SysError(func, "Reallocation of F_I_L_E failed!");

      for( i = fileID; i < fileTableSize; i++ )
	fileTableInitEntry(i);
    }

  fileTable[fileID].used = TRUE;

  return (fileID);
}

static void fileTableDelete(void)
{
  static char func[] = "fileTableDelete";
  int fileID;

  for (fileID = 0; fileID < fileTableSize; fileID++)
    if ( fileTable[fileID].used ) fileClose(fileID);

  if ( fileTable )
    {
      free(fileTable);
      fileTableSize = 0;
    }
}

static long fileGetenv(char *envName)
{
  static char func[] = "fileGetenv";
  char *envString;
  long envValue = -1;
  long fact = 1;

  envString = getenv(envName);

  if ( envString )
    {
      int loop;

      for ( loop = 0; loop < (int) strlen(envString); loop++ )
	{
	  if ( ! isdigit((int) envString[loop]) )
	    {
	      switch ( tolower((int) envString[loop]) )
		{
		case 'k':  fact = 1024;        break;
		case 'm':  fact = 1048576;     break;
		case 'g':  fact = 1073741824;  break;
		default:
		  fact = 0;
		  Message(func, "Invalid number string in %s: %s", envName, envString);
		  Warning(func, "%s must comprise only digits [0-9].",envName);
		}
	      break;
	    }
	}

      if ( fact ) envValue = fact*atol(envString);

      if ( FILE_Debug ) Message(func, "Set %s to %ld", envName, envValue);
    }

  return (envValue);
}

static void fileInitialize(void)
{
  static char func[] = "fileInitialize";
  long value;

  FileInit = 1;

  value = fileGetenv("FILE_DEBUG");
  if ( value >= 0 ) FILE_Debug = (int) value;

  FileInfo = (int) fileGetenv("FILE_INFO");

  value  = fileGetenv("FILE_BUFSIZE");
  if ( value >= 0 ) FileBufferSizeEnv = value;

  value = fileGetenv("FILE_TYPE");
#if ! defined (HAVE_MMAP)
  if ( value == FILE_TYPE_MMAP )
    {
      Warning(func, "MMAP not available!");
      value = 0;
    }
#endif
  if ( value > 0 )
    {
      switch (value)
	{
	case FILE_TYPE_STD:
	case FILE_TYPE_MMAP:
	  FileBufferTypeEnv = value;
	  break;
	default:
	  Warning(func, "File type %d not implemented!", value);
	}
    }

  if ( FILE_Debug ) atexit(fileTablePrint);

  atexit(fileTableDelete);
}

static void fileSetBuffer(int fileID)
{
  static char func[] = "fileSetBuffer";
  size_t buffersize = 0;

#if defined (FILE_CHECKID)
  fileCheckID(func, fileID);
#endif

  if ( fileTable[fileID].mode == 'r' )
    {
      if ( FileBufferTypeEnv )
	fileTable[fileID].bufferType = FileBufferTypeEnv;
      else if ( fileTable[fileID].bufferType == 0 )
	fileTable[fileID].bufferType = FILE_TYPE_STD;

      if ( FileBufferSizeEnv >= 0 )
	buffersize = FileBufferSizeEnv;
      else if ( fileTable[fileID].bufferSize > 0 )
	buffersize = fileTable[fileID].bufferSize;
      else
	{
	  buffersize = fileTable[fileID].blockSize * 4;
	  if ( buffersize < FileBufferSizeMin ) buffersize = FileBufferSizeMin;
	}

      if ( (size_t) fileTable[fileID].size < buffersize )
	buffersize = fileTable[fileID].size;

      if ( fileTable[fileID].bufferType == FILE_TYPE_MMAP )
	{
	  size_t blocksize = pagesize();
	  size_t minblocksize = 4 * blocksize;
	  buffersize = buffersize - buffersize % minblocksize;

	  if ( buffersize < (size_t) fileTable[fileID].size && buffersize < minblocksize )
	    buffersize = minblocksize;
	}
    }
  else
    {
      fileTable[fileID].bufferType = FILE_TYPE_STD;

      if ( FileBufferSizeEnv >= 0 )
	buffersize = FileBufferSizeEnv;
      else if ( fileTable[fileID].bufferSize > 0 )
	buffersize = fileTable[fileID].bufferSize;
      else
	{
	  buffersize = fileTable[fileID].blockSize * 4;
	  if ( buffersize < FileBufferSizeMin ) buffersize = FileBufferSizeMin;
	}
    }

  if ( buffersize == 0 ) buffersize = 1;

  if ( fileTable[fileID].bufferType == FILE_TYPE_STD )
    {
      fileTable[fileID].buffer = (char *) malloc(buffersize);
      if ( fileTable[fileID].buffer == NULL )
	SysError(func, "Allocation of file buffer failed!");
    }	

  if ( fileTable[fileID].mode != 'r' )
    if ( setvbuf(fileTable[fileID].fp, fileTable[fileID].buffer, _IOFBF, buffersize) )
      SysError(func, "setvbuf failed");

  fileTable[fileID].bufferSize = buffersize;
}

static int fileFillBuffer(int fileID)
{
  static char func[] = "fileFillBuffer";
  long nread;
  int fd;
  int ret;
  long offset = 0;
  off_t retseek;
  
  if ( FILE_Debug )
    Message(func, "fileID = %d  Cnt = %ld", fileID, fileTable[fileID].bufferCnt);

  if ( (fileTable[fileID].flag & FILE_EOF) != 0 ) return (EOF);

  if ( fileTable[fileID].buffer == NULL ) fileSetBuffer(fileID);
  
  if ( fileTable[fileID].bufferSize == 0 ) return (EOF);

  fd = fileTable[fileID].fd;

#if defined (HAVE_MMAP)
  if ( fileTable[fileID].bufferType == FILE_TYPE_MMAP )
    {
      if ( fileTable[fileID].bufferPos >= fileTable[fileID].size )
	{
	  nread = 0;
	}
      else
	{
	  nread = fileTable[fileID].bufferSize;
	  if ( (nread + fileTable[fileID].bufferPos) > fileTable[fileID].size )
	    nread = fileTable[fileID].size - fileTable[fileID].bufferPos;

	  if ( fileTable[fileID].buffer )
	    {
	      ret = munmap(fileTable[fileID].buffer, fileTable[fileID].mappedSize);
	      if ( ret == -1 )
		SysError(func, "munmap error for read %s", fileTable[fileID].name);
	      fileTable[fileID].buffer = NULL;
	    }

	  fileTable[fileID].mappedSize = (size_t) nread;

	  fileTable[fileID].buffer =
            (char *) mmap(0, (size_t) nread, PROT_READ, MAP_SHARED, fd, fileTable[fileID].bufferPos);

	  if ( (int) fileTable[fileID].buffer == -1 )
	    SysError(func, "mmap error for read %s", fileTable[fileID].name);

	  offset = fileTable[fileID].position - fileTable[fileID].bufferPos;
	}
    }
  else
#endif
    {
      retseek = lseek(fileTable[fileID].fd, fileTable[fileID].bufferPos, SEEK_SET);
      if ( retseek == (off_t)-1 )
	SysError(func, "lseek error at pos %ld file %s", (long) fileTable[fileID].bufferPos, fileTable[fileID].name);
	
      nread = (long) read(fd, fileTable[fileID].buffer, fileTable[fileID].bufferSize);
    }

  if ( nread <= 0 )
    {
      if ( nread == 0 )
	fileTable[fileID].flag |= FILE_EOF;
      else
	fileTable[fileID].flag |= FILE_ERROR;

      fileTable[fileID].bufferCnt = 0;
      return (EOF);
    }

  fileTable[fileID].bufferPtr = fileTable[fileID].buffer;
  fileTable[fileID].bufferCnt = nread;

  fileTable[fileID].bufferStart = fileTable[fileID].bufferPos;
  fileTable[fileID].bufferPos  += nread;
  fileTable[fileID].bufferEnd   = fileTable[fileID].bufferPos - 1;

  if ( FILE_Debug )
    {
      Message(func, "fileID = %d  Val     = %d",  fileID, (int) fileTable[fileID].buffer[0]);
      Message(func, "fileID = %d  Start   = %ld", fileID, fileTable[fileID].bufferStart);
      Message(func, "fileID = %d  End     = %ld", fileID, fileTable[fileID].bufferEnd);
      Message(func, "fileID = %d  nread   = %ld", fileID, nread);
      Message(func, "fileID = %d  offset  = %ld", fileID, offset);
      Message(func, "fileID = %d  Pos     = %ld", fileID, fileTable[fileID].bufferPos);
      Message(func, "fileID = %d  postion = %ld", fileID, fileTable[fileID].position);
    }

  if ( offset > 0 )
    {
      if ( offset > nread )
	Error(func, "Internal problem with buffer handling. nread = %d offset = %d", nread, offset);

      fileTable[fileID].bufferPtr += offset;
      fileTable[fileID].bufferCnt -= offset;
    }

  fileTable[fileID].bufferNumFill++;

  return ((unsigned char) *fileTable[fileID].bufferPtr);
}

static void fileCopyFromBuffer(int fileID, void *ptr, size_t size)
{
  static char func[] = "fileCopyFromBuffer";

  if ( FILE_Debug )
    Message(func, "size = %ld  Cnt = %ld", size, fileTable[fileID].bufferCnt);

  if ( fileTable[fileID].bufferCnt < size )
    Error(func, "Buffer to small. bufferCnt = %d", fileTable[fileID].bufferCnt);

  if ( size == 1 )
    {
      ((char *)ptr)[0] = fileTable[fileID].bufferPtr[0];

      fileTable[fileID].bufferPtr++;
      fileTable[fileID].bufferCnt--;
    }
  else
    {
      memcpy(ptr, fileTable[fileID].bufferPtr, size);

      fileTable[fileID].bufferPtr += size;
      fileTable[fileID].bufferCnt -= size;
    }
}

static size_t
fileReadFromBuffer(int fileID, void *ptr, size_t size)
{
  static char func[] = "fileReadFromBuffer";
  size_t nread, rsize;
  size_t offset = 0;

  if ( FILE_Debug )
    Message(func, "size = %ld  Cnt = %d", size, (int) fileTable[fileID].bufferCnt);

  if ( ((int)fileTable[fileID].bufferCnt) < 0 )
    Error(func, "Internal problem. bufferCnt = %d", (int) fileTable[fileID].bufferCnt);

  rsize = size;

  while ( fileTable[fileID].bufferCnt < rsize )
    {
      nread = fileTable[fileID].bufferCnt;
      /*
      fprintf(stderr, "rsize = %d nread = %d\n", (int) rsize, (int) nread);
      */
      if ( nread > (size_t) 0 )
	fileCopyFromBuffer(fileID, (char *)ptr+offset, nread);
      offset += nread;
      if ( nread < rsize )
	rsize -= nread;
      else
	rsize = 0;

      if ( fileFillBuffer(fileID) == -1 ) break;
    }

  nread = size - offset;

  if ( fileTable[fileID].bufferCnt < nread ) nread = fileTable[fileID].bufferCnt;

  if ( nread > (unsigned) 0 )
    fileCopyFromBuffer(fileID, (char *)ptr+offset, nread);

  return (nread+offset);
}

void fileSetBufferSize(int fileID, long buffersize)
{
  fileTable[fileID].bufferSize = buffersize;
}

/* 
 *   Open a file. Returns file ID, or -1 on error
 */
int fileOpen(const char *filename, const char *mode)
{
  static char func[] = "fileOpen";
  FILE *fp = NULL;    /* file pointer    (used for write) */
  int fd = -1;        /* file descriptor (used for read)  */
  int fileID = FILE_UNDEFID;
  int fmode = 0;

  if ( ! FileInit ) fileInitialize();

  fmode = tolower((int) mode[0]);

  switch ( fmode )
    {
    case 'r':  fd =  open(filename, O_RDONLY); break;
    case 'w':  fp = fopen(filename, "w");      break;
    case 'a':  fp = fopen(filename, "a");      break;
    default:   Error(func, "Mode %c unexpected\n", fmode);
    }

  if ( FILE_Debug )
    if ( fp == NULL && fd == -1 )
      Message(func, "Open failed on %s mode %c", filename, fmode);

  if ( fp )
    {
      struct stat filestat;

      fileID = fileTableNewEntry();

      fileTable[fileID].fp   = fp;
      fileTable[fileID].mode = fmode;
      fileTable[fileID].name = strdup(filename);

      if ( stat(fileTable[fileID].name, &filestat) != 0 )
	SysError(func, fileTable[fileID].name);

      fileTable[fileID].blockSize = filestat.st_blksize;

      if ( fmode == 'r' )
	fileTable[fileID].size = filestat.st_size;

      if ( FILE_Debug )
	Message(func, "File %s opened with ID %d", filename, fileID);
    }
  else if ( fd >= 0 )
    {
      struct stat filestat;

      fileID = fileTableNewEntry();

      fileTable[fileID].fd   = fd;
      fileTable[fileID].mode = fmode;
      fileTable[fileID].name = strdup(filename);
 
      if ( fstat(fd, &filestat) != 0 ) SysError(func, filename);

      fileTable[fileID].blockSize = filestat.st_blksize;

      if ( fmode == 'r' ) fileTable[fileID].size = filestat.st_size;

      if ( FILE_Debug )
	Message(func, "File %s opened with ID %d", filename, fileID);      
    }

  return (fileID);
}

/* 
 *   Close a file.
 */
void fileClose(int fileID)
{
  static char func[] = "fileClose";
  char *name;
  int ret;
  char *ftname[] = {"unknown", "standard", "mmap"};

  fileCheckID(func, fileID);

  name = fileTable[fileID].name;

  if ( FILE_Debug )
    Message(func, "fileID = %d  filename = %s", fileID, name);

  if ( FileInfo > 0 )
    {
      fprintf(stderr, "____________________________________________\n");
      fprintf(stderr, " file ID          : %d\n",  fileID);
      fprintf(stderr, " file name        : %s\n",  fileTable[fileID].name);
      if ( fileTable[fileID].mode == 'r' )
	fprintf(stderr, " file descriptor  : %d\n",  fileTable[fileID].fd);
      else
	fprintf(stderr, " file pointer     : %p\n",  (void *) fileTable[fileID].fp);

      fprintf(stderr, " file mode        : %c\n",  fileTable[fileID].mode);

      if ( sizeof(off_t) > sizeof(long) )
	{
	  fprintf(stderr, " file size        : %lld\n", (long long) fileTable[fileID].size);
	  fprintf(stderr, " file position    : %lld\n", (long long) fileTable[fileID].position);
	  fprintf(stderr, " bytes transfered : %lld\n", (long long) fileTable[fileID].byteTrans);
	}
      else
	{
	  fprintf(stderr, " file size        : %ld\n", (long) fileTable[fileID].size);
	  fprintf(stderr, " file position    : %ld\n", (long) fileTable[fileID].position);
	  fprintf(stderr, " bytes transfered : %ld\n", (long) fileTable[fileID].byteTrans);
	}

      fprintf(stderr, " file access      : %ld\n", fileTable[fileID].access);
      fprintf(stderr, " buffer type      : %d (%s)\n",  fileTable[fileID].bufferType, ftname[fileTable[fileID].bufferType]);
      fprintf(stderr, " num buffer fill  : %ld\n", fileTable[fileID].bufferNumFill);
      fprintf(stderr, " buffer size      : %lu\n", (unsigned long) fileTable[fileID].bufferSize);
      fprintf(stderr, " block size       : %lu\n", (unsigned long) fileTable[fileID].blockSize);
#if defined (HAVE_MMAP)
      fprintf(stderr, " page size        : %d\n",  pagesize());
#endif
      fprintf(stderr, "--------------------------------------------\n");
    }

  if ( fileTable[fileID].mode == 'r' )
    {
#if defined (HAVE_MMAP)
      if ( fileTable[fileID].buffer && fileTable[fileID].mappedSize )
	{
	  ret = munmap(fileTable[fileID].buffer, fileTable[fileID].mappedSize);
	  if ( ret == -1 )
	    SysError(func, "munmap error for close %s", fileTable[fileID].name);
	  fileTable[fileID].buffer = NULL;
	}
#endif
      ret = close(fileTable[fileID].fd);
      if ( ret == -1 )
	SysError(func, "EOF returned for close of %s\n", name);
    }
  else
    {
      ret = fclose(fileTable[fileID].fp);
      if ( ret == EOF )
	SysError(func, "EOF returned for close of %s\n", name);
    }

  if ( fileTable[fileID].name )    free((void*) fileTable[fileID].name);
  if ( fileTable[fileID].buffer )  free((void*) fileTable[fileID].buffer);

  fileTableInitEntry(fileID);
}

int fileGetc(int fileID)
{
  int ivalue = -1;
  int fillret = 0;

  if ( fileTable[fileID].bufferCnt == 0 ) fillret = fileFillBuffer(fileID);

  if ( fillret >= 0 )
    {
      ivalue = (unsigned char) *fileTable[fileID].bufferPtr++;
      fileTable[fileID].bufferCnt--;

      fileTable[fileID].position++;
      fileTable[fileID].byteTrans++;
      fileTable[fileID].access++;
    }

  return (ivalue);
}

size_t fileRead(int fileID, void *ptr, size_t size)
{
  static char func[] = "fileRead";
  size_t nread;
  
#if defined (FILE_CHECKID)
  fileCheckID(func, fileID);
#endif

  nread = fileReadFromBuffer(fileID, ptr, size);

  fileTable[fileID].position  += nread;
  fileTable[fileID].byteTrans += nread;
  fileTable[fileID].access++;

  if ( FILE_Debug ) Message(func, "size %ld  nread %ld", size, nread);

  return (nread);
}

size_t fileWrite(int fileID, const void *ptr, size_t size)
{
#if defined (FILE_CHECKID)
  static char func[] = "fileWrite";
#endif
  size_t nwrite;
  FILE *fp;

#if defined (FILE_CHECKID)
  fileCheckID(func, fileID);
#endif

  if ( fileTable[fileID].buffer == NULL ) fileSetBuffer(fileID);

  fp = fileTable[fileID].fp;

  nwrite = fwrite(ptr, 1, size, fp);

  fileTable[fileID].position  += nwrite;
  fileTable[fileID].byteTrans += nwrite;
  fileTable[fileID].access++;

  return (nwrite);
}
/*
// pbio.c
*/

#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif




#define BUFFLEN 4096

/*
// Default buffer size for I/O operations (set via setvbuf)
*/
#define SIZE BUFSIZ

static long size = SIZE;
static int sizeSet = 0;
static char  *envSize;

/*
// Debug flags.
*/
#define DEBUGOFF 1
#define DEBUG (debugSet > DEBUGOFF)
static char *debugLevel;
static int debugSet = 0;

static int  rawSize = 0;
static int *rawFlag = NULL;

/*
//------------------------------------------------------------------------
// pbopen - Open file
//------------------------------------------------------------------------
*/
void pbOpen(int *unit, const char *name, const char *mode, int *iret)
{
  /*
  // Purpose:
  //  Opens file, returns the index of a UNIX FILE pointer held in 
  //  an internal table (fileTable).
  //
  // First time through, reads value in environment variable PBIO_BUFSIZE
  // (if it is set) and uses it as the size to be used for internal file
  // buffers; the value is passed to setvbuf. If PBIO_BUFSIZE is not set,
  // a default value is used.
  //
  // Function  accepts:
  //    name = filename
  //    mode = r, w
  //
  // Function returns:
  //   INTEGER iret:
  //     -1 = Could not open file.
  //     -2 = Invalid file name.
  //     -3 = Invalid open mode specified
  //      0 = OK.
  */

  /*
  // See if DEBUG switched on.
  */
  if ( ! debugSet )
    {
      debugLevel = getenv("PBIO_DEBUG");
      if ( debugLevel == NULL )
        debugSet = DEBUGOFF;              /* off */
      else
	{
	  int loop;
	  for ( loop = 0; loop < (int) strlen(debugLevel) ; loop++ )
	    {
	      if ( ! isdigit((int) debugLevel[loop]) )
		{
		  printf("Invalid number string in PBIO_DEBUG: %s\n", debugLevel);
		  printf("PBIO_DEBUG must comprise only digits [0-9].\n");
		  debugSet = DEBUGOFF;
		}
	    }
	  debugSet = DEBUGOFF + atoi(debugLevel);
	}
      if ( DEBUG ) printf("PBIO_PBOPEN: debug switched on\n");
    }

  *unit = 0;
  *iret = 0;

  if ( DEBUG ) printf("PBIO_PBOPEN: file name = %s\n", name);

  switch (*mode)
    {
    case 'a':
    case 'A':
      break;
    case 'c':
    case 'C':
    case 'w':
    case 'W':
      break;
    case 'r':
    case 'R':
      break;
    default:
      *iret = -3;
      return;
    }

  if ( DEBUG ) printf("PBIO_PBOPEN: file open mode = %c\n", *mode);

  /*
  // Now allocate a buffer for the file, if necessary.
  */
  if ( ! sizeSet )
    {
      envSize = getenv("PBIO_BUFSIZE");
      if ( envSize == NULL )
        size = SIZE;             /* default */
      else
	{
	  int loop;
	  for ( loop = 0; loop < (int) strlen(envSize) ; loop++ )
	    {
	      if ( ! isdigit((int) envSize[loop]) )
		{
		  printf("Invalid number string in PBIO_BUFSIZE: %s\n", envSize);
		  printf("PBIO_BUFSIZE must comprise only digits [0-9].\n");
		  exit(1);
		}
	    }
	  size = atol(envSize);
	}
      if ( size <= 0 )
	{
	  printf("Invalid buffer size in PBIO_BUFSIZE: %s\n", envSize);
	  printf("Buffer size defined by PBIO_BUFSIZE must be positive.\n");
	  exit(1);
	}
      sizeSet = 1;
    }

  *unit = fileOpen(name, mode);

  if ( *unit == -1 ) *iret = -1;
  else
    {
      if ( size >= 0 ) fileSetBufferSize(*unit, size);

      if ( rawSize == 0 )
	{
	  rawSize = 8;
	  rawFlag = (int *) malloc(rawSize*sizeof(int));
	}
      if ( *unit >= rawSize )
	{
	  rawSize = *unit*2;
	  rawFlag = (int *) realloc(rawFlag, rawSize*sizeof(int));
	}

      rawFlag[*unit] = 1;
    }

  if ( DEBUG ) printf("PBIO_PBOPEN: file ID = %d\n", *unit);

  if ( DEBUG ) printf("PBIO_PBOPEN: file buffer size = %ld\n", size);
}

/*
//------------------------------------------------------------------------
// pbseek - Seek
//------------------------------------------------------------------------
*/
void pbSeek(int unit, int offset, int whence, int *iret)
{
  /*
  //
  // Purpose:
  //   Seeks to a specified location in file.
  //
  //  Function  accepts:
  //    unit = the index of a UNIX FILE pointer held in
  //           an internal table (fileTable).
  //
  //    offset = byte count
  //
  //    whence  = 0, from start of file
  //            = 1, from current position
  //            = 2, from end of file.  
  //
  //  Returns:
  //    iret:
  //      -2 = error in handling file,
  //      -1 = end-of-file
  //      otherwise,  = byte offset from start of file.
  */
  int my_offset = (int) offset;
  int my_whence = (int) whence;

  /*
  // Must use negative offset if working from end-of-file
  */
  if ( DEBUG )
    { 
      printf("PBIO_PBSEEK: file ID = %d\n", unit);
      printf("PBIO_PBSEEK: Offset = %d\n", my_offset);
      printf("PBIO_PBSEEK: Type of offset = %d\n", my_whence);
    }

  if ( my_whence == 2 ) my_offset = - abs(my_offset);

  *iret = (int) fileGetPos(unit);
  if ( DEBUG ) printf("PBIO_PBSEEK: current position = %d\n", *iret);
  if ( *iret == my_offset && my_whence == 0 )
    *iret = 0;
  else
    *iret = fileSetPos(unit, (off_t) my_offset, my_whence);

  if ( DEBUG ) printf("PBIO_PBSEEK: fileSetPos return code = %d\n", *iret);

  if ( *iret != 0 )
    {
      if ( ! fileEOF(unit) )
	{
	  *iret = -2;             /* error in file-handling */
	  perror("pbseek");
	}
      else
        *iret = -1;             /* end-of-file  */

      fileClearerr(unit);
      return;
    }

  /*
  // Return the byte offset from start of file
  */
  *iret = (int) fileGetPos(unit);

  if ( DEBUG )
    printf("PBIO_PBSEEK: byte offset from start of file = %d\n", *iret);

  return;
}

/*
//------------------------------------------------------------------------
// PBTELL - Tells current file position (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbTell(int unit, int *iret)
{
  /*
  //
  // Purpose:
  //   Tells current byte offset in file.
  //
  //  Function  accepts:
  //    unit = the index of a UNIX FILE pointer held in
  //           an internal table (fileTable).
  //
  //  Returns:
  //    iret:
  //      -2 = error in handling file,
  //      otherwise,  = byte offset from start of file.
  */


  /*
  // Return the byte offset from start of file
  */
  *iret = (int) fileGetPos(unit);

  if ( *iret < 0 )
    {
      if ( DEBUG )
	{           /* error in file-handling */
	  printf("PBIO_PBTELL: file ID = %d. ", unit);
	  printf("Error status = %d\n", *iret);
	}
      perror("pbtell");
      *iret = -2;
    }

  if ( DEBUG )
    {
      printf("PBIO_PBTELL: file ID = %d. ", unit);
      printf("Byte offset from start of file = %d\n",*iret);
    }

  return;
}
		
/*
//------------------------------------------------------------------------
//  pbread - Read
//------------------------------------------------------------------------
*/
void pbRead(int unit, void *buffer, int nbytes, int *iret)
{
  /*
  // Purpose:
  //  Reads a block of bytes from a file..
  //
  //  Function  accepts:
  //    unit = the index of a UNIX FILE pointer held in
  //           an internal table (fileTable).
  //
  //    nbytes = number of bytes to read.
  //
  //  Returns:
  //    iret:
  //      -2 = error in reading file,
  //      -1 = end-of-file,
  //      otherwise, = number of bytes read.
  */

  if ( DEBUG )
    {
      printf("PBIO_READ: file ID = %d. ", unit);
      printf("Number of bytes to read = %d\n", nbytes);
    }

  if ( (*iret = fileRead(unit, buffer, nbytes) ) != nbytes )
    {
      if ( ! fileEOF(unit) )
	{
	  *iret = -2;             /*  error in file-handling  */
	  perror("pbread");
	  fileClearerr(unit);
	  return;
	}
      else
	{
	  *iret = -1;             /*  end-of-file */
	  fileClearerr(unit);
	}
    }

    if ( DEBUG )
      {
	printf("PBIO_READ: file ID = %d. ", unit);
	printf("Number of bytes read = %d\n", nbytes);
      }

    return;
}
 
/*
//------------------------------------------------------------------------
//  pbwrite - Write (from FORTRAN)
//------------------------------------------------------------------------
*/
void pbWrite(int unit, void *buffer, int nbytes, int *iret)
{
  /*
  // Purpose:
  //  Writes a block of bytes to a file.
  //
  //  Function  accepts:
  //    unit = the index of a UNIX FILE pointer held in
  //           an internal table (fileTable).
  //
  //    nbytes = number of bytes to write.
  //
  //  Returns:
  //    iret:
  //      -1 = Could not write to file.
  //     >=0 = Number of bytes written.
  */
  unsigned char cbuf[4];

  if ( DEBUG )
    {
      printf("PBIO_WRITE: file ID = %d. ", unit);
      printf("Number of bytes to write = %d\n", nbytes);
    }

  if ( rawFlag[unit] == 0 )
    {
      /* If we are at the beginning of the file, write 4 zero bytes */
      if ( fileGetPos(unit) == 0 )
	{
	  cbuf[0] = cbuf[1] = cbuf[2] = cbuf[3] = 0;
	  if ( fileWrite(unit, cbuf, 4) != 4 ) perror("pbwrite");
	}

      cbuf[0] = (nbytes>>24) & 0xff;
      cbuf[1] = (nbytes>>16) & 0xff;
      cbuf[2] = (nbytes>> 8) & 0xff;
      cbuf[3] = (nbytes    ) & 0xff;

      /* 4 header Bytes with the record length (BIG Endian) */
      if ( fileWrite(unit, cbuf, 4) != 4 ) perror("pbwrite");
    }

  if ( (*iret = fileWrite(unit, buffer, nbytes) ) != nbytes)
    {
      perror("pbwrite");
      *iret = -1;
    }

  if ( rawFlag[unit] == 0 )
    {
      /* 4 trailing bytes */
      if ( fileWrite(unit, cbuf, 4) != 4 ) perror("pbwrite");
    }

  if ( DEBUG )
    {
      printf("PBIO_WRITE: file ID = %d. ", unit);
      printf("PBIO_WRITE: number of bytes written = %d\n", *iret);
    }

  return;
}

/*
//------------------------------------------------------------------------
//   pbclose - close
//------------------------------------------------------------------------
*/
void pbClose(int unit, int *iret)
{
  /*
  // Purpose:
  //  Closes file.
  //
  //  Function  accepts:
  //    unit = the index of a UNIX FILE pointer held in
  //           an internal table (fileTable).
  ////  Returns:
  //    iret:
  //      0 = OK.
  //      otherwise = error in handling file.
  */

  if( DEBUG )
    printf("PBIO_CLOSE: file ID = %d\n", unit);

  if ( rawFlag[unit] == 0 )
    {
      unsigned char cbuf[4];

      /* Write 4 trailing 0 Bytes */

      cbuf[0] = cbuf[1] = cbuf[2] = cbuf[3] = 0;
      fileWrite(unit, cbuf, 4);
    }

  *iret = 0;

  fileClose(unit);

  return;
}

/*
//------------------------------------------------------------------------
//  pbflush - flush
//------------------------------------------------------------------------
*/
void pbFlush(int unit)
{
  /*
  // Purpose:	Flushes file.
  */

  if ( DEBUG )
    printf("PBIO_FLUSH: file ID = %d\n", unit);

  fileFlush(unit);
}

/*
//------------------------------------------------------------------------
//  GRIBREAD
//------------------------------------------------------------------------
*/
void gribread(void *buffer, int buffsize, int *readsize, int *status, int unit)
{
  /*
  //  Called as a FORTRAN subroutine:
  //
  //    CALL GRIBREAD( KARRAY, KINLEN, KOUTLEN, IRET, KUNIT )
  //
  */
  size_t holdsize = buffsize;

  /*
  // Read GRIB product
  */

  *status = gribRead(unit, (unsigned char *) buffer, &holdsize);

  *readsize = (int) holdsize;

  if ( DEBUG )
    {
      printf("PBIO_GRIBREAD: file ID = %d. ", unit);
      printf("Number of bytes read = %d\n", *readsize);
    }

  return;
}

/*
//------------------------------------------------------------------------
//  PBSIZE
//------------------------------------------------------------------------
*/
void pbSize(int unit, int *plen)
{
  /*
  //  Returns the size in bytes of the next GRIB product.
  //
  //  Called from FORTRAN:
  //      CALL PBSIZE(KUNIT, LENGTH)
  //
  //  unit  = file id returned from PBOPEN.
  //  plen  = size in bytes of the next product.
  //        = -2 if error allocating memory for internal buffer.
  //
  //  The input file is left positioned where it started.
  */
  off_t recpos;
  long offset;
  int ierr;
  int recsize;

  /*
  //  Use a smallish buffer for processing; this should suffice for all cases
  //  except versions -1 and 0 of GRIB and large BUFR products
  */

  recpos = fileGetPos(unit);

  if ( DEBUG )
    {
      printf("PBIO_SIZE: file ID = %d. ", unit);
      printf("Current file position = %ld\n", (long) recpos);
    }

  ierr = gribFileSeek(unit, &offset);
  if ( ierr > 0 ) printf("GRIB record not found\n");

  if ( ierr == 0 )
    {
      recsize = gribReadSize(unit);
      *plen = recsize;
    }
  else
    *plen = 0;

  /*
  //  Put the file pointer back where it started
  */
  if ( DEBUG )
    {
      printf("PBIO_SIZE: file pointer set back to: %ld\n", (long) recpos);
      printf("PBIO_SIZE: Product size = %d\n", *plen);
    }

  (void) fileSetPos(unit, recpos, SEEK_SET);

  return ;
}

void pbGrib(int unit, void *buffer, int inlen, int *outlen, int *iret)
{
  /*
    PBGRIB

    PURPOSE
    _______

    Reads next GRIB product from a file.

    INTERFACE
    _________

    CALL PBGRIB(KUNIT,KARRAY,KINLEN,KOUTLEN,KRET)


    Input parameters
    ________________

    KUNIT  - 	Unit number for the file returned from PBOPEN
    KARRAY -	FORTRAN array big enough to hold the GRIB product
    KINLEN -	Size in BYTES of the FORTRAN array


    Output parameters
    ________________

    KOUTLEN -	Size in BYTES of the GRIB product read	
    KRET    -	 0  if a GRIB product has been successfully read
                -1  if end-of-file is hit before a GRIB product is read
		-3  if the size of KARRAY is not sufficient for the 
                   GRIB product


    Common block usage
    __________________

    None.


    Method
    ______

    Calls GRIBREAD.


    Externals
    _________

    GRIBREAD - Read next GRIB product.


    AUTHOR
    ______

    J.D.Chambers       ECMWF


    MODIFICATIONS
    _____________

    None.
  */
  int nread;
  int retval;

  /*
    Get the GRIB product
  */
  gribread(buffer, inlen, &nread, &retval, unit);

  /*
    Escape if the user buffer is too small to hold even the early sections of
    the product or EOF encountered
  */
  if ( retval == -4 )
    {
      *outlen = nread;
      *iret = -1;
      return;
    }
  /*
    Escape if no GRIB product is found in the file
  */
  if ( retval == -1 )
    {
      *outlen = 0;
      *iret = -1;
      return;
    }
  /*
    Check if the array is big enough for the GRIB product
  */
  if ( retval == -3 )
    {
      *outlen = nread;
      *iret = -3;
      return;
    }
  /*
    Set success code if product retrieved
  */
  if ( nread >= 0 )
    {
      *outlen = nread;
      *iret = 0;
    }

  return;
}

void pbSetRaw(int unit, int raw)
{
  rawFlag[unit] = raw;
}

#if ! defined   (_GRIBFORTRAN_H)
#  include "gribFortran.h"
#endif

#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined (HAVE_CF_INTERFACE)
#  include "cfortran.h"
#endif

/*
  Rename PBIO routines if RENAME_PB is set
*/
#if  defined  (RENAME_PB)
#  define  pbopen   pb_open
#  define  pbclose  pb_close
#  define  pbread   pb_read
#  define  pbwrite  pb_write
#  define  pbseek   pb_seek
#  define  pbtell   pb_tell
#  define  pbflush  pb_flush
#  define  pbsize   pb_size
#  define  pbgrib   pb_grib
#endif

/*
    Rename GRIB routines if RENAME_GR is set
*/
#if  defined  (RENAME_GR)
#  define  gribex   grib_ex
#  define  grprs0   gr_prs0
#  define  grprs1   gr_prs1
#  define  grprs2   gr_prs2
#  define  grprs3   gr_prs3
#  define  grprs4   gr_prs4
#  define  grsdbg   gr_sdbg
#  define  grsrnd   gr_srnd
#  define  grsref   gr_sref
#  define  grsvck   gr_svck
#endif

/*
  PDOUBLE on CRAY is 128bit
*/
#if  defined  (CRAY)
#  define  PDOUBLE  PFLOAT
#endif

#if defined (HAVE_CF_INTERFACE)
/*
  PBIO
*/
FCALLSCSUB4 (pbOpen,   PBOPEN,   pbopen,  PINT, STRING, STRING, PINT)
FCALLSCSUB2 (pbClose,  PBCLOSE,  pbclose, INT, PINT)
FCALLSCSUB4 (pbRead,   PBREAD,   pbread,  INT, PVOID, INT, PINT)
FCALLSCSUB4 (pbWrite,  PBWRITE,  pbwrite, INT, PVOID, INT, PINT)
FCALLSCSUB4 (pbSeek,   PBSEEK,   pbseek,  INT, INT, INT, PINT)
FCALLSCSUB2 (pbTell,   PBTELL,   pbtell,  INT, PINT)
FCALLSCSUB1 (pbFlush,  PBFLUSH,  pbflush, INT)
FCALLSCSUB2 (pbSize,   PBSIZE,   pbsize,  INT, PINT)
FCALLSCSUB5 (pbGrib,   PBGRIB,   pbgrib,  INT, PVOID, INT, PINT, PINT)
FCALLSCSUB2 (pbSetRaw, PBSETRAW, pbsetraw,  INT, INT)

/*
  GRIBEX Print
*/
FCALLSCSUB1 (gribprs0,    GRPRS0,    grprs0,    PINT)
FCALLSCSUB2 (gribprs1,    GRPRS1,    grprs1,    PINT, PINT)

FCALLSCSUB3 (gribprs2,    GRPRS2,    grprs2,    PINT, PINT, PVOID)
FCALLSCSUB3 (gribprs2dp,  GRPRS2DP,  grprs2dp,  PINT, PINT, PDOUBLE)
FCALLSCSUB3 (gribprs2sp,  GRPRS2SP,  grprs2sp,  PINT, PINT, PFLOAT)

FCALLSCSUB3 (gribprs3,    GRPRS3,    grprs3,    PINT, PINT, PVOID)
FCALLSCSUB3 (gribprs3dp,  GRPRS3DP,  grprs3dp,  PINT, PINT, PDOUBLE)
FCALLSCSUB3 (gribprs3sp,  GRPRS3SP,  grprs3sp,  PINT, PINT, PFLOAT)

FCALLSCSUB3 (gribprs4,    GRPRS4,    grprs4,    PINT, PINT, PVOID)
FCALLSCSUB3 (gribprs4dp,  GRPRS4DP,  grprs4dp,  PINT, PINT, PDOUBLE)
FCALLSCSUB3 (gribprs4sp,  GRPRS4SP,  grprs4sp,  PINT, PINT, PFLOAT)

/*
  GRIBEX
*/
FCALLSCSUB14 (gribEx,   GRIBEX,    gribex,    PINT, PINT, PINT, PVOID ,  PINT,
	       PVOID,   PINT, PVOID,   INT,   PINT,  INT, PINT, STRING,  PINT)
FCALLSCSUB14 (gribExDP, GRIBEXDP,  gribexdp,  PINT, PINT, PINT, PDOUBLE, PINT,
               PDOUBLE, PINT, PDOUBLE,  INT,  PINT,  INT, PINT, STRING,  PINT)
FCALLSCSUB14 (gribExSP, GRIBEXSP,  gribexsp,  PINT, PINT, PINT, PFLOAT,  PINT,
               PFLOAT,  PINT, PFLOAT,  INT,   PINT,  INT, PINT, STRING,  PINT)

/*
  GRIBEX defaults
*/
FCALLSCSUB1 (gribsdbg,    GRSDBG,    grsdbg,    INT)
FCALLSCSUB1 (gribsrnd,    GRSRND,    grsrnd,    INT)
FCALLSCSUB1 (gribsref,    GRSREF,    grsref,    PVOID)
FCALLSCSUB1 (gribsrefsp,  GRSREFSP,  grsrefsp,  FLOAT)
FCALLSCSUB1 (gribsrefdp,  GRSREFDP,  grsrefdp,  DOUBLE)
FCALLSCSUB1 (gribsvck,    GRSVCK,    grsvck,    INT)

#endif

void gribEx(int *isec0, int *isec1, int *isec2, void *fsec2, int *isec3,
            void *fsec3, int *isec4, void *fsec4, int klenp, int *kgrib,
            int kleng, int *kword, char *hoper, int *kret)
{
#if defined (FORT_REAL4)
  gribExSP(isec0, isec1, isec2, (float *) fsec2, isec3,
	   (float *) fsec3, isec4, (float *) fsec4, klenp, kgrib,
	   kleng, kword, hoper, kret);
#else
  gribExDP(isec0, isec1, isec2, (double *) fsec2, isec3,
	   (double *) fsec3, isec4, (double *) fsec4, klenp, kgrib,
	   kleng, kword, hoper, kret);
#endif
}

void gribprs2(int *isec0, int *isec2, void *fsec2)
{
#if defined (FORT_REAL4)
  gribprs2sp(isec0, isec2, (float *) fsec2);
#else
  gribprs2dp(isec0, isec2, (double *) fsec2);
#endif
}

void gribprs3(int *isec0, int *isec3, void *fsec3)
{
#if defined (FORT_REAL4)
  gribprs3sp(isec0, isec3, (float *) fsec3);
#else
  gribprs3dp(isec0, isec3, (double *) fsec3);
#endif
}

void gribprs4(int *isec0, int *isec4, void *fsec4)
{
#if defined (FORT_REAL4)
  gribprs4sp(isec0, isec4, (float *) fsec4);
#else
  gribprs4dp(isec0, isec4, (double *) fsec4);
#endif
}

void gribprs0(int *isec0)
{
  gribPrintSec0(isec0);
}

void gribprs1(int *isec0, int *isec1)
{
  gribPrintSec1(isec0, isec1);
}

void gribprs2dp(int *isec0, int *isec2, double *fsec2)
{
  gribPrintSec2DP(isec0, isec2, fsec2);
}

void gribprs2sp(int *isec0, int *isec2, float  *fsec2)
{
  gribPrintSec2SP(isec0, isec2, fsec2);
}

void gribprs3dp(int *isec0, int *isec3, double *fsec3)
{
  gribPrintSec3DP(isec0, isec3, fsec3);
}

void gribprs3sp(int *isec0, int *isec3, float  *fsec3)
{
  gribPrintSec3SP(isec0, isec3, fsec3);
}

void gribprs4dp(int *isec0, int *isec4, double *fsec4)
{
  gribPrintSec4DP(isec0, isec4, fsec4);
}

void gribprs4sp(int *isec0, int *isec4, float  *fsec4)
{
  gribPrintSec4SP(isec0, isec4, fsec4);
}

void gribprs4w(int *isec4)
{
  gribPrintSec4Wave(isec4);
}

void gribsdbg(int debug)
{
  gribSetDebug(debug);
}

void gribsrnd(int round)
{
  gribSetRound(round);
}

void gribsref(void *fref)
{
#if defined (FORT_REAL4)
  gribSetRefSP(*((float *) fref));
#else
  gribSetRefDP(*((double *) fref));
#endif
}

void gribsrefdp(double fref)
{
  gribSetRefDP(fref);
}

void gribsrefsp(float fref)
{
  gribSetRefSP(fref);
}

void gribsvck(int vcheck)
{
  gribSetValueCheck(vcheck);
}

 
 

 


 


 


 


 


 


 


 


 


 


 


 


 


 


 


 


 


 


 


 


 


 
 



 



char *LIBGRIB = "0.5.6"  " of ""Aug 16 2004"" ""13:55:18";


 




static const char grb_libvers[] = "0.5.6"  " of ""Aug 16 2004"" ""13:55:18";



const char *
gribLibraryVersion(void)
{
  return (grb_libvers);
}
