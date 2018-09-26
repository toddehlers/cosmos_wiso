MODULE mod_gauss
!
! -- gauss.h   29-09-95   Version 2.0   Author: Laurent Terray
!    *******   25-10-96   Version 2.1   Addition of amskred
!              16-03-99   Version 2.3   Addition of T213 and T319
!              26-03-99   Version 2.3   Changed troncature for number of 
!                                        latitude between equator and pole
!@
!@  Contents : variables related to gaussian troncature (if any spectral model)
!@  --------
!@
!@ -ninip16: number of longitudes for the 16 latitude circles (T21) (1D)
!@
!@ -ninip24: number of longitudes for the 24 latitude circles (T31) (1D)
!@
!@ -ninip32: number of longitudes for the 32 latitude circles (T42) (1D)
!@
!@ -ninip48: number of longitudes for the 48 latitude circles (T63) (1D)
!@
!@ -ninip80: number of longitudes for the 80 latitude circles (T106) (1D)
!@
!@ -ninip160: number of longitudes for the 160 latitude circles (T213-319) (1D)
!@
!@ -- nredu16 : number of points on T21 reduced gaussian grid
!@
!@ -- nredu24 : number of points on T31 reduced gaussian grid
!2
!@ -- nredu32 : number of points on T42 reduced gaussian grid
!@
!@ -- nredu48 : number of points on T63 reduced gaussian grid
!@
!@ -- nredu80 : number of points on T106 reduced gaussian grid
!@
!@ -- nredu160 : number of points on T213-T319 reduced gaussian grids
!@
!@ -- amskred : mask value for reduced grid
!@
!     -------------------------------------------------------------------
!
  USE mod_kinds_oasis
! - T21
  INTEGER (kind=ip_intwp_p), dimension(16), PARAMETER :: ninip16 = &
     (/20,30,40,48,54,60,64,64,64,64,64,64,64,64,64,64/)
! - T31
  INTEGER (kind=ip_intwp_p), dimension(24), PARAMETER :: ninip24 = &
     (/20,30,40,48,54,60,64,64,64,64,64,64,64,64,64,64, &
     96,96,96,96,96,96,96,96/)
! - T42
  INTEGER (kind=ip_intwp_p), dimension(32), PARAMETER :: ninip32 = &
     (/20,30,40,48,48,54,54,64,72,80,80,90,96,100, &
     108,108,120,120,120,128,128,128,128,128, &
     128,128,128,128,128,128,128,128/)
! - T63
  INTEGER (kind=ip_intwp_p), dimension(48), PARAMETER :: ninip48 = &
     (/16, 16, 18, 24, 30, 36, 48, 48, 54, 60, 72, 72, 80, 90, 90 &
     ,96,100,108,120,120,128,128,144,144,144,144,150,160,160,160 &
     ,162,180,180,180,180,180,180,192,192,192,192,192,192,192,192 &
     ,192,192,192/)
! - T106
  INTEGER (kind=ip_intwp_p), dimension(80), PARAMETER :: ninip80 = &
     (/16, 16, 18, 24, 30, 36, 48, 50, 60, 64, 72, 80, 80, 90, 96 &
     ,100,108,120,120,128,128,144,144, 144,150,160,162,180,180,180 &
     ,192,192,192,200,216,216,216,216, 240,240,240,240,240,250,250 &
     ,250,256,270,270,270,270,288,288, 288,288,288,288,300,300,300 &
     ,300,300,320,320,320,320,320,320, 320,320,320,320,320,320,320 &
     ,320,320,320,320,320/)
!- T213- T319
  INTEGER (kind=ip_intwp_p), dimension(160), PARAMETER  :: ninip160 = &
     (/18, 25, 36, 40, 45, 50, 60, 64, 72, 72, 80, 90, 90, 96,108 &
     ,120,120,125,128,135,144,150,160,160,180,180,180,192,192,200 &
     ,216,216,225,225,240,240,243,250,256,270,270,288,288,288,300 &
     ,300,320,320,320,320,324,360,360,360,360,360,360,375,375,375 &
     ,384,384,400,400,400,405,432,432,432,432,432,450,450,450,450 &
     ,480,480,480,480,480,480,480,500,500,500,500,500,512,512,540 &
     ,540,540,540,540,540,540,540,576,576,576,576,576,576,576,576 &
     ,576,576,600,600,600,600,600,600,600,600,600,640,640,640,640 &
     ,640,640,640,640,640,640,640,640,640,640,640,640,640,640,640 &
     ,640,640,640,640,640,640,640,640,640,640,640,640,640,640,640 &
     ,640,640,640,640,640,640,640,640,640,640/)
  INTEGER (kind=ip_intwp_p), PARAMETER  :: nredu16 = 1784
  INTEGER (kind=ip_intwp_p), PARAMETER  :: nredu24 = 3320 
  INTEGER (kind=ip_intwp_p), PARAMETER  :: nredu32 = 6232
  INTEGER (kind=ip_intwp_p), PARAMETER  :: nredu48 = 12228 
  INTEGER (kind=ip_intwp_p), PARAMETER  :: nredu80 = 23444 
  INTEGER (kind=ip_intwp_p), PARAMETER  :: nredu160 = 138346
  REAL (kind=ip_realwp_p), PARAMETER  :: amskred = 9999999.0
!
END MODULE mod_gauss
!     -------------------------------------------------------------------




