!*    NAMELIST *JSBGRID_CTL* To control the grid in jsbach stand alone runs
!                                                                 
!      
!      NAME            TYPE        PURPOSE                                       
!      ----            ----        -------                                       
!
!     *LON_EAST*      REAL        Eastern boundary of region
!     *LON_WEST*      REAL        Western boundary of region
!     *LAT_NORTH*     REAL        Northern boundary of region
!     *LAT_SOUTH*     REAL        Southern boundary of region
!     *NPROCA*        INTEGER     Number of processors in longitude direction
!     *NPROCB*        INTEGER     Number of processors in latitude direction
!     *NPEDIM*        INTEGER     Working dimension for blocks in each domain.
!                                 This could e.g. be set to the vector length
!                                 on a vector machine. The default (-1):
!                                 each domain is processed in one call

NAMELIST /JSBGRID_CTL/ &
      LON_EAST, LON_WEST, LAT_NORTH, LAT_SOUTH, NPROCA, NPROCB, NPEDIM
!                                                                       
!     ------------------------------------------------------------------
