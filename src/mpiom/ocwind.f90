SUBROUTINE ocwind

!sv**** *ocwind* -  update ocean velocities with the surface stress.
!sv
!sv     modified
!sv     --------
!sv     stephan venzke,    mpi, 1999
!sv       26.08.99
!sv       - included fluxes option (i.e. distinguish bewteen tau over water and ice)
!sv         here only tau over water is considered (i.e. txo is replaced by aofltxwo)
!sv         tau over ice is used in sbr ocice
!sv       02.09.99
!sv       - included common-block mo_fluxes1 for forcing with atmospheric fluxes
!sv         *common*    *fluxes1* - atmospheric fields on the ocean grid.
!sv ------------------------------------------------------------------------------
!
!        uwe  08.03.00
!           include tauwat, ocean ice stress from ocice
  USE mo_param1
  USE mo_parallel
  USE mo_commo1
  USE mo_commoau1

#ifdef __coupled
  USE mo_fluxes1
#endif


#ifdef bounds_exch_tp
  REAL   pqx,pqy,tye_g(ie_g,je_g),tauwatv_g(ie_g,je_g)
#endif


!$OMP PARALLEL PRIVATE(I,J,UR,VR)
  dtdect=dt/dzw(1)
!#ifdef bounds_exch_save
!$OMP SINGLE
  call bounds_exch('p',sicomo,'ocwind 1')
!$OMP END SINGLE
!#endif

#ifdef bounds_exch_tp
!$OMP SINGLE
!!  call bounds_exch('vf',tye,'ocwind 2')
!!  call bounds_exch('vf',tauwatv,'ocwind 3')
!$OMP END SINGLE

!$OMP SINGLE
      call gather_arr(tye,tye_g,0)
      call gather_arr(tauwatv,tauwatv_g,0)
      call p_bcast(tye_g,0)
      call p_bcast(tauwatv_g,0)
      do i=2,ie_g/2-1
      iop=ie_g+1-i
      pqx=0.5*(tye_g(i,2)-tye_g(iop,2))
      pqy=0.5*(tauwatv_g(i,2)-tauwatv_g(iop,2))
      tye_g(i,2)=pqx
       tauwatv_g(i,2)=pqy
      tye_g(iop,2)=-pqx
      tauwatv_g(iop,2)=-pqy
      enddo
      call scatter_arr(tye_g,tye,p_io)
      call scatter_arr(tauwatv_g,tauwatv,p_io)
!$OMP END SINGLE
#endif


!$OMP DO
  DO j=2,je-1
     DO i=2,ie-1
!hh use sea ice concentration from previous timestep to by consistent to sea ice dynamics 
        ur=1.-0.5*(sicomp(i,j)+sicomp(i+1,j))    
        vr=1.-0.5*(sicomp(i,j)+sicomp(i,j+1))

!uwe include ice stress on water
!sv 26.8.99 included fluxes option (i.e. distinguish bewteen tau over water and ice)
!sv         here only tau over water is considered (i.e. txo is replaced by aofltxwo)
!sv         tau over ice is used in sbr ocice

#ifdef __coupled
!sv 26.8.99 achtung: bei stefanie sind die windstresse ueber wasser schon mit der 
!sv                  mittleren compactness (gemittelt ueber ganze atmospaeren-zelle) gewichtet !!! 
!sv         d.h. hier muesste ur eigentlich durch ur/(gesamtcompcness in atmospaeren-zelle) ersetzt werden!

        uoo(i,j,1)=uoo(i,j,1)+dtdect*amsuo(i,j,1)*(aofltxwo(i,j)*ur       &
                  +(1.-ur)*amsuo(i,j,1)*tauwatu(i,j))

        voe(i,j,1)=voe(i,j,1)+dtdect*amsue(i,j,1)*(aofltywe(i,j)*vr       &
                  +(1.-vr)*amsue(i,j,1)*tauwatv(i,j))
#else

        uoo(i,j,1)=uoo(i,j,1)+dtdect*amsuo(i,j,1)*(txo(i,j)*ur            &
                  +(1.-ur)*amsuo(i,j,1)*tauwatu(i,j))

        voe(i,j,1)=voe(i,j,1)+dtdect*amsue(i,j,1)*(tye(i,j)*vr            &
             +(1.-vr)*amsue(i,j,1)*tauwatv(i,j))
#endif
     END DO
  END DO
!$OMP END DO
!#ifdef bounds_exch_save
!$OMP SINGLE
  CALL bounds_exch('u',uoo,'ocwind 2')
  CALL bounds_exch('v',voe,'ocwind 3')
!$OMP END SINGLE
!#endif
!$OMP END PARALLEL

END SUBROUTINE ocwind



