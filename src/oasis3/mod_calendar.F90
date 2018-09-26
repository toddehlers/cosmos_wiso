MODULE mod_calendar
!
! -- calendar.h   01-11-95   Version 2.0   Author: Laurent Terray
!    **********
!@
!@  Contents : variables related to the coupler calendar
!@  --------
!@
!@ -- nddeb : beginning date of the simulation (yyyymmdd)
!@            (always 00000000 in blkdata.f, not present in namcouple)
!@
!@ -- nadeb : beginning year of the simulation (yy)
!@            (always 00 as ndded always 00000000)
!@
!@ -- nmdeb : beginning month of the simulation (mm)
!@            (always 00 as ndded always 00000000)
!@
!@ -- njdeb : beginning day of the simulation (dd)
!@            (always 00 as ndded always 00000000)
!@
!@ -- ndate : initial date (yyyymmdd)
!@
!@ -- njini : initial day (dd)
!@
!@ -- nmini : initial month (mm)
!@
!@ -- naini : initial year (yy)
!@
!@ -- njnow : current day (dd)
!@ 
!@ -- njone : inf. day limit for linear time interpolation (dd)
!@              'sea <<<--->>> sea' case
!@
!@ -- njtwo : sup. day limit for linear time interpolation (dd)
!@              'sea <<<--->>> sea' case
!@
!@ -- ndone : inf. day limit for linear time interpolation (dd)
!@              'sea <<<--->>> ice' case
!@
!@ -- ndtwo : sup. day limit for linear time interpolation (dd)
!@              'sea <<<--->>> ice' case
!@
!@ -- nmnow : current month (mm)
!@
!@ -- nmone : inf. month limit for linear time interpolation (mm)
!@
!@ -- nmtwo : sup. month limit for linear time interpolation (mm)
!@
!@ -- nanow : current year (yy)
!@
!@ -- ndinc : day increment for each coupler time step
!@
!@ -- nsrec : number of records to be skipped in climatology (SE) file
!@
!@ -- nmrec : number of records to be skipped in interannual (MO) file
!@
!@ -- ncaltype : calendar type
!@                 0      = 365 day calendar (no leap years)
!@                 1      = Gregorian calendar
!@                 n (>1) = n day month caledar
!@
!@     -------------------------------------------------------------------
  USE mod_kinds_oasis
  INTEGER(kind=ip_intwp_p) :: ndate, njini, nmini, naini, njnow,njone, njtwo
  INTEGER(kind=ip_intwp_p) :: ndtwo, nmnow,nmone, nmtwo, nanow, ndinc, nsrec
  INTEGER(kind=ip_intwp_p) :: ndone, nmrec
  INTEGER(kind=ip_intwp_p) :: ncaltype, nddeb, nadeb, nmdeb, njdeb
  INTEGER(kind=ip_intwp_p) :: ig_date(6)
!
END MODULE mod_calendar
!     -------------------------------------------------------------------
