      SUBROUTINE SUCHIJ(abrei,alaij,k,ipos,jpos,dist,wetref)
!
!     subroutine to determine i and j indices of nearest point
!     from lat and long
!
! RJ: i and j are GLOBAL indices!!!
!
!     input:
!     abrei      latitude (deg.)
!     alaij      longitude (deg.)
!     k          level (integer)
!     wetref     if land shall be considered too : 0.
!                ocean points only                 1.
!
!     output:
!     ipos       i index
!     jpos       j index
!     dist       distance from point in m
!
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS
      IMPLICIT NONE

      REAL abrei,alaij,dist,wetref
      INTEGER k,ipos,jpos

      REAL pipi,grabog,aphi,alam,xx,yy,zz,ax,ay,az,d
      INTEGER i,j
!
      pipi=2.*asin(1.)
      grabog=pipi/180.
!
      aphi=grabog*abrei
      alam=grabog*alaij
!
      xx=cos(alam)*cos(aphi)
      yy=sin(alam)*cos(aphi)
      zz=sin(aphi)
!
      dist = 1.e20
      ipos = 0
      jpos = 0
!
!CDIR NOLOOPCHG
      do j=2,je_g-1
!CDIR NOVECTOR
       do i=2,ie_g-1
        if(weto_g(i,j,k).gt.wetref-0.1)then
          ax=cos(giph_g(2*i,2*j))*cos(gila_g(2*i,2*j))
          ay=cos(giph_g(2*i,2*j))*sin(gila_g(2*i,2*j))
          az=sin(giph_g(2*i,2*j))
          d=(xx-ax)**2+(yy-ay)**2+(zz-az)**2
          if (d<dist) then
            dist = d
            ipos = i
            jpos = j
          endif
        endif
       enddo
      enddo

      dist=sqrt(dist)*6350000.
      if (p_pe == p_io) then
        write(6,*)'in suchij end: ',ipos,jpos,k,weto_g(ipos,jpos,k)
      endif
      return
      end
