subroutine tcheck(l,feld)
  USE MO_MPI
  USE MO_PARAM1
  USE MO_PARALLEL
  USE MO_COMMO1
  USE MO_UNITS

  IMPLICIT NONE
  real helfie_g(ie_g,je_g),vq(ie_g),feld(ie,je,ke)
  integer i,k,ir,l,id1(1),id2(1)
   if(l.le.50)return
    vq=0.
    call gather_arr(feld(:,:,1),helfie_g(:,:),p_io)
      IF (p_pe==p_io) THEN

    do i=2,ie_g/2 - 1
    ir=ie_g+1-i
    vq(i)=helfie_g(i,2)-helfie_g(ir,3)
    enddo
!    write(0,*)'linie',(i,vq(i),i=1,ie_g/2)
     id1=minloc(vq)
     id2=maxloc(vq)
    i=id1(1)
    k=id2(1)
    write(0,*)'vcheck',l,i,k,minval(vq),maxval(vq),helfie_g(i,2)&
   &, helfie_g(k,2)
    ENDIF
    return
    end
