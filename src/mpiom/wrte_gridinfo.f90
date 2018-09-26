      Subroutine wrte_gridinfo(kky,kkm,kkd)
!
!c**** *wrte_gridinfo* - save grid information.
!c
!c     ch,    *mpi-met, hh*    10.04.01
!
!     modified
!     --------
!     s.legutke,        *mpi-mad, hh*    01.10.01
!     - separate routine extracted from ollie (main)
!     - netcdf version possible (with cond.comp. pnetcdfo)
!
!     purpose
!     -------
!     
!
!     method
!     -------
!     
!
!**   interface.
!     ----------
!
!     *call*       *wrte_gridinfo(kky,kkm,kkd)*
!
!     *parameter*  *param1.h*     - grid size parameters for ocean model.
!     *common*     *commo1.h*     - ocean/sediment tracer arrays.
!     *common*     *units.h*      - std i/o logical units.
!
!**   interface to calling routine (parameter list):
!     ----------------------------------------------
!
!     *integer* *kky*   - actual year.
!     *integer* *kkm*   - actual month.
!     *integer* *kkd*    - actual day.
!
!
!     externals
!     ---------
!     none.
!
!**************************************************************************

      Use mo_param1,Only :    ie_g,je_g,ke
      Use mo_parallel
      Use mo_commo1, Only : dlxp,dlyp,dlxu,dlyu,dlxv,dlyv,depto,deuto,weto
      Use mo_units
      Use mo_kind


      Integer(kind=i4) i4i1,i4i2,i4i3,i4i4,i4i5
      Integer kky,kkm,kkd

! write to disk (extra format)
!
          i4i1=((kky*10000)+(kkm*100)+kkd)
          i4i3=0
          i4i4=(ie_g*je_g)
          i4i2=85

          If (p_pe == p_io) Then 
             Write(io_ou_dlxp)i4i1,i4i2,i4i3,i4i4
          Endif
          Call write_slice_sp(io_ou_dlxp,dlxp)

          i4i2=86
          If (p_pe == p_io) Then 
             Write(io_ou_dlyp)i4i1,i4i2,i4i3,i4i4
          Endif
          Call write_slice_sp(io_ou_dlyp,dlyp)

          i4i2=185
          If (p_pe == p_io) Then 
             Write(io_ou_dlxu)i4i1,i4i2,i4i3,i4i4
          Endif
          Call write_slice_sp(io_ou_dlxu,dlxu)

          i4i2=186
          If (p_pe == p_io) Then 
          Write(io_ou_dlyu)i4i1,i4i2,i4i3,i4i4
          Endif
          Call write_slice_sp(io_ou_dlyu,dlyu)

          i4i2=184
          If (p_pe == p_io) Then 
          Write(io_ou_deut)i4i1,i4i2,i4i3,i4i4
          Endif
          Call write_slice_sp(io_ou_deut,deuto)

          i4i2=84
          If (p_pe == p_io) Then 
          Write(io_ou_dept)i4i1,i4i2,i4i3,i4i4
          Endif
          Call write_slice_sp(io_ou_dept,depto)

          Do k=1,ke
          i4i2=172
          i4i3=k
          If (p_pe == p_io) Then 
          Write(io_ou_weto)i4i1,i4i2,i4i3,i4i4
          Endif
          Call write_slice_sp(io_ou_weto,weto(:,:,k))
          Enddo


       Return
       End
