      module mod_flux_01y2
!
!     calculate total and AB=01(v,j,Omega) resolved fluxes
!
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      use mod_absorcion_01y2
      use mod_photoini_01y2
      use mod_colini_01y2
     
      implicit none

      save
* paraS2
      real*8 :: photonorm
      real*8 :: r1flux_ang,ekinmin_ev,ekinmax_ev
      real*8 :: r1flux,ekinmin,ekinmax
      integer :: netot,ncontfile,icanref
! total flux
      integer :: ir4flux        ! r1 distance index at which total flux is evaluated
      complex*16,allocatable ::  zAR1flux(:) !vector to evaluate derivatives in r1 at r1(ir4flux)
      complex*16,allocatable :: zCkcheby(:)
      complex*16,allocatable :: zfR1flux(:,:,:),zR(:),zdR(:)      
      integer,allocatable :: ir1flux(:)
      integer,allocatable :: indangflux(:),indcanflux(:),indr2flux(:)
      complex*16, allocatable :: zCR(:,:), zCdR(:,:)
      real*8,allocatable :: S2prodtot(:),S2prod(:),reacfct(:)
      real*8,allocatable :: apaqini(:)
      
! flux in individual channels in the 01  + 2 reactant channel
      
      integer :: ir2balint,nflux,nfluxprocdim
      real*8, allocatable :: rpopi(:,:,:), funreac(:,:,:)
      real*8, allocatable :: S2(:,:,:)
      
      complex*16, allocatable :: zfR2(:,:,:) ,zfR1(:,:,:)
      complex*16,allocatable :: zS2(:,:,:,:)


      real*8,allocatable :: etotS2(:),S2factor(:,:,:,:)
     &           ,Cvj(:,:,:),Cvjproc(:,:,:)
      integer,allocatable :: nfluxproc(:)

      contains
********************************************
*   functions of   mod_flux_01y2.f    *
********************************************
!=======================================================

!--------------------------------------------------
      subroutine ini_flux
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2

      implicit none

      include "mpif.h"

      integer :: ir1,ir2,iang,iang_proc,ican,ican_proc,ir,ip
      integer :: ielec,iflux,j,iv,ie
      real*8 :: r1,abstot,r2,xk,delta2r,box1
      integer :: ierror,iloop0,indt,it0,mcan,m,k
      integer*8 :: imem

************************************************************
      namelist /inputflux/r1flux_ang,netot,ekinmin_eV,ekinmax_eV
     &                   ,ncontfile
************************************************************

!  reading data in namelist

      r1flux=absr1-1.d0
      ncontfile=0
           
      write(6,'(40("="),/,10x,"flux_mod",/,40("="))')
      write(6,*)
      write(6,*)'  flux data'
      write(6,*)'  ---------'
      open(10,file='input.dat',status='old')
         read(10,nml = inputflux)
         write(6,nml = inputflux)
         call flush(6)
      close(10)

      ekinmin=ekinmin_eV*ev2cm*conve1
      ekinmax=ekinmax_eV*ev2cm*conve1
! for reactants
      allocate(rpopi(nvini:nvmax,jini:jmax,ncanmax)
     &             ,funreac(nanguprocdim,npun1,ncanmax)  
     &             ,S2prodtot(netot),S2prod(netot)
     &             ,reacfct(netot)
     &             ,zfR2(npun1,nanguprocdim,ncanmax)
     &             ,zCkcheby(netot),apaqini(netot)
     &     , stat=ierror)

!     determining ir2balint for 01(v,j,Omega) flux

      do ir2=1,npun2
         r2=rmis2+dble(ir2-1)*ah2
         if(r2.lt.absr2.and.ir2.gt.ir2balint)ir2balint=ir2
      enddo

      if(npun1.gt.1)then
         r1flux=r1flux_ang
         ir4flux=0
         do ir1=1,npun1
            r1=rmis1+dble(ir1-1)*ah1
            if(r1.le.r1flux.and.ir1.gt.ir4flux)ir4flux=ir1
         enddo
         if(ir4flux.lt.0.and.ir4flux.gt.npun1)then
            write(6,*)'  ir4flux = ',ir4flux,' out of range '
            call flush(6)
            call MPI_BARRIER(MPI_COMM_WORLD, ierror)
            stop
         endif
         write(6,*)' total flux evaluated at index ir4flux= ',ir4flux
      
!     calculating vector to perform derivative d/dr1 with fft at ir4flux
      
         allocate(zAR1flux(npun1)
     &       , stat=ierror)
         norealproc_mem=norealproc_mem
     &    +npun1*2+netot*2

         box1=dble(npun1-1)*ah1
         call fftmom(box1,npun1,npun1,xm1red,hbr,xpr1,xp2r1)
         do m=1,npun1
            delta2r=dble(ir4flux-m)*ah1
            zAR1flux(m)=dcmplx(0.d0,0.d0)
            do k=1,npun1
               xk=xpr1(k)
               zAR1flux(m) = zAR1flux(m)
     &               + dcmplx(0.d0,xk)*cdexp(dcmplx(0.d0,xk*delta2r))    
            enddo
            zAR1flux(m) = zAR1flux(m)/dble(2*npun1)
         enddo
      
! setting ir1flux(iang) at which flux is analized

         allocate(ir1flux(nangu),nfluxproc(0:nproc-1)
     &       , stat=ierror)
         nointegerproc_mem=nointegerproc_mem+nangu+nproc
     
         ir1flux(:)=0

         nflux=0
         do ip=0,nproc-1
            nfluxproc(ip)=0
            do iang_proc=1,nanguproc
               iang=indangreal(iang_proc,ip)
           
               do ican_proc=1,ncanproc(ip)
                  ican=ibasproc(ican_proc,ip)
! ir1flux: distances in r1 at which flux is evaluated
                  do ir=1,npunreal(iang)
                     ir1=indr1(ir,iang)
                     r1=dble(ir1-1)*ah1+rmis1
                     if(r1.le.r1flux.and.ir1.gt.ir1flux(iang))
     &                                      ir1flux(iang)=ir1

                  enddo ! ir

! nfluxproc(ip) and nflux
                  do ir=1,npunreal(iang)
                     ir1=indr1(ir,iang)
                     if(ir1.eq.ir1flux(iang).and.npun1.gt.1)then
                        nflux=nflux+1
                        nfluxproc(ip)=nfluxproc(ip)+1
                     endif
                  enddo
               enddo ! ican_proc
            enddo ! iang_proc
!          write(6,*)' for proc ',ip,'  nfluxproc= ',nfluxproc(ip)
         enddo ! ip

!  allocating matrices

         nfluxprocdim=nfluxproc(idproc)
     
          allocate( zCR(netot,nfluxprocdim)
     &             ,zCdR(netot,nfluxprocdim)
     &             ,zfR1(npun2,nanguprocdim,ncanmax)
     &             ,zfR1flux(npun2,nanguprocdim,ncanmax)
     &             ,zR(npun2),zdR(npun2)
     &             ,indangflux(nfluxprocdim)
     &             ,indcanflux(nfluxprocdim)
     &             ,indr2flux(nfluxprocdim)
     &       , stat=ierror)
         if(ierror.ne.0)then
           write(*,*)" error in initmem for flux "
           stop
         else
           nointegerproc_mem=nointegerproc_mem+nfluxprocdim*3
           imem=
     &      netot*nfluxprocdim*2*2  ! complex number 
     &       + (npun1+2*npun2)*nanguprocdim*ncanmax*2 ! complex number
     &       + npun2*2*2 ! complex number
     &       +netot*3
     &       +(nvmax-nvini+1)*(jmax-jini+1)*ncanmax
     &       +nanguprocdim*npun1*ncanmax       
            norealproc_mem=norealproc_mem+imem
            write(6,*)' zCR,zCdR,etc in iniflux, proc= '
     &          ,idproc
     &          ,' imem=',imem
     &          ,' memory(Gb)= '
     &          ,dble(imem*8)*1.d-9
            call flush(6)
         endif
      
         indangflux(:)=0
         indcanflux(:)=0
         indr2flux(:)=0
       
         zCR(:,:)=dcmplx(0.d0,0.d0)
         zCdR(:,:)=dcmplx(0.d0,0.d0)
         S2prodtot(:)=0.d0
         reacfct(:)=0.d0

         zfR2(:,:,:)=dcmplx(0.d0,0.d0)
         zfR1(:,:,:)=dcmplx(0.d0,0.d0)
         zfR1flux(:,:,:)=dcmplx(0.d0,0.d0)
    
!     determining indexes to evaluate flux

         do ip=idproc,idproc
            iflux=0
            do iang_proc=1,nanguproc
               iang=indangreal(iang_proc,ip)
           
               do ican_proc=1,ncanproc(ip)
                  ican=ibasproc(ican_proc,ip)
                  do ir=1,npunreal(iang)
                     ir1=indr1(ir,iang)
                     ir2=indr2(ir,iang)
                     if(ir1.eq.ir1flux(iang).and.npun1.gt.1)then
                        iflux=iflux+1
                        indangflux(iflux)=iang
                        indcanflux(iflux)=ican
                        indr2flux(iflux)=ir2
                     endif
                  enddo
               enddo ! ican_proc
            enddo ! iang_proc
             write(6,*)' for proc ',ip,'  nfluxproc= ',nfluxproc(ip)
         enddo ! ip

      endif ! npun1>1

!     starting the allocation of S2 paraS2
         if(iprod.eq.0.or.npun1.eq.1)then
            nviniprod=nvini
            nvmaxprod=nvmax
            jiniprod=jini
            jmaxprod=jmax
            iomminprod=iommin
            iommaxprod=iommax
         endif
         if(iprod.ge.0)then
            ierror=0
         
            allocate(etotS2(netot)
     &     , S2factor(netot,nvini:nvmax,jini:jmax,nelecmax)
     &     , S2(nvini:nvmax,jini:jmax,ncanmax)
     &     , zS2(netot,nvini:nvmax,jini:jmax,ncanmax)
     &     , Cvj(nvini:nvmax,jini:jmax,ncanmax)
     &     , Cvjproc(nvini:nvmax,jini:jmax,ncanmax)
     &       , stat=ierror)
            if(ierror.ne.0)then
              write(*,*)" error in initmem for paraS2reac "
              call flush(6)
              stop
           else
              norealproc_mem=norealproc_mem
     &          +netot*(1+(nvmax-nvini+1)*(jmax-jini+1)*nelecmax)
     &          +(nvmax-nvini+1)*(jmax-jini+1)*nelecmax*3
     &          +netot*(nvmax-nvini+1)*(jmax-jini+1)*nelecmax*2              
               write(6,*)'norealproc_mem= ',norealproc_mem
               write(6,*)'nointegerproc_mem= ',nointegerproc_mem
               call flush(6)
            endif
        
         endif ! iprod>0

**>> reading information for continuing propagation

      it0=0
      indt=0
      iloop0=0
      open(3,file='cont.data',status='old',err=1)
      read(3,*)it0,indt,iloop0
 1    continue
      
      if(it0.eq.0)then
         abstot=0.d0
         if(npun1.gt.1)zS2(:,:,:,:)=dcmplx(0.d0,0.d0)
      else
         read(3,*)abstot,nflux,mcan
         close(3)
         if(ncontfile.eq.1.and.npun1.gt.1)then
            write(name,'("cont",i2.2,".S2")')idproc
            open(4,file=name,status='old',form='unformatted')
            read(4)zS2
            write(6,*)' zS2 read in ',name
            read(4)zCR
            write(6,*)' zCR read in ',name
            read(4)zCdR
            write(6,*)' zCdR read in ',name
            close(4)
         endif
      endif

      call Smatini

      write(6,*)' ending ini_flux in proc= ',idproc
     &     ,' norealproc_mem=',norealproc_mem
     &     ,' nointegerproc_mem=',nointegerproc_mem
     &     ,' memory(Gb)= '
     &   ,dble(nointegerproc_mem*4+norealproc_mem*8)*1.d-9
      call flush(6)
      
      return
      end subroutine ini_flux

!--------------------------------------------------

      subroutine Smatini
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2

      implicit none

      include "mpif.h"
      integer :: ican,ierror

       icanref=0
       do ican=1,ncan
          if(nelebas(ican).eq.ielecref
     &             .and.iombas(ican).eq.iomref)icanref=ican
       enddo
       if(icanref.ge.1.and.icanref.le.ncan)then
          if(idproc.eq.0)write(6,*) '  --> reference channel is ican= '
     &                                ,icanref
          if(idproc.eq.0)write(6,*) '        for ielec ',ielecref
     &                                ,'and Omega= ',iomref
       else
          write(6,*)'  icanref can not be set'
          call flush(6)
          call MPI_BARRIER(MPI_COMM_WORLD, ierror)
          stop
       endif

! For Iphoto=0

       if(iphoto.eq.0)then
          write(6,'(/,80("-"),/,4x
     &        ,"Building Smatrix auxiliar for collisions"
     &        ,/,80("-"),//)')
          call Smatini_col

! For Iphoto > 0
       elseif(iphoto.gt.0)then
          write(6,'(/,80("-"),4x
     &        ,"Building Smatrix auxiliar for photoinitiated process"
     &        ,/,80("-"),//)')
          call Smatini_photo
      endif     
      
      return
      end subroutine Smatini

!--------------------------------------------------

      subroutine Smatini_col
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2

      implicit none

      include "mpif.h"
      integer, parameter :: npunt=1000
      real*8 :: rgaussr(npunt)
      real*8 :: ekinfin,estep,etotmax,etotmin,pfin,ahgauss,arg,ekinini
      real*8 :: df,dg,f,g,paqini,pepe,pini,r,r2,xnorm,xr2,abstot,erot
      real*8 :: rfin,wro
      integer :: ie,iele,iv,j,ir2,key
      complex*16 :: zexpo,zfft

*> Emin,Emax kinetic energies to determine interval for S2 elements evaluation

       etotmin=ekinmin
       etotmax=ekinmax
       estep=(etotmax-etotmin)/dble(netot-1)

**> For Iphoto=0
       photonorm=1.d0
       if(idproc.eq.0)then
             write(6,*)' Kin. Ener. range to evaluate S^2 (eV):'
     &             ,etotmin/conve1/8065.5d0,etotmax/conve1/8065.5d0
       endif


*1              building initial gaussian
*a)  estimated closest l:

       if(idproc.eq.0)open(43,file='gaussR',status='unknown')
       ahgauss=(rfin2-rmis2)/dble(npunt-1)
       do ir2=1,npunt
          r2=rmis2+dble(ir2-1)*ahgauss
          call  rgauscolini(xr2,r2,rcolini,alpha0,xk0,factor0,il0)
          rgaussr(ir2)=xr2/dsqrt(2.d0)
       if(idproc.eq.0)write(43,*)r2,rgaussr(ir2)*rgaussr(ir2)
       enddo
       if(idproc.eq.0)close(43)

***>> S^2 matrix magnitudes: 
*   imposing bessel behaviour for the closest l on a single Omega channel:
*    consistent with a CS approach in which a single Omega is considered
*    The more general problem of considering body-fixed bessel function"
*                   (see JCP, 113,(2000),1781 and refs. therein)
*            will be considered later

      key=-1
      if(idproc.eq.0)write(6,"(/,10x,'For Jtot=',i3
     &         ,' j= ',i3,' Omg= ',i3,' v=',i2,
     &     /,15x,' the closest value of l= ',i3,//
     &      )")Jtot,jref,iomref,nvref,il0

*2)  calculating initial energy distribution

      if(idproc.eq.0)open(43,file='gaussE',status='unknown')
      xnorm=0.d0
      pepe=dble(il0*(il0+1))
      rfin=rcolini+1.d0
      do ie=1,netot
         etotS2(ie)=etotmin+dble(ie-1)*estep
         ekinini=etotS2(ie)
         if(ekinini.gt.0.d0)then
            pini=dsqrt(ekinini*xmasa0*2.d0/hbr/hbr)
            zfft=dcmplx(0.d0,0.d0)
            do ir2=1,npunt
               r=rmis2+dble(ir2-1)*ahgauss
               erot=hbr*hbr*pepe*0.5d0/(rcolini*rcolini*xmasa0)
!               erot=hbr*hbr*pepe*0.5d0/(r*r*xmasa0)

               if(ekinini.gt.erot)then
                 arg=r*pini
                 CALL BESPH2(F,DF,G,DG,PEPE,ARG,KEY,0)
                 zexpo=dcmplx(-g,f)
                 zfft=zfft+dcmplx(rgaussr(ir2),0.d0)*zexpo
               endif
            enddo
            zfft=zfft*ahgauss/2.d0/pi
            paqini=dreal(zfft*dconjg(zfft))
         else
            pini=0.d0
            paqini=1.d0
         endif
         apaqini(ie)=paqini
** reactants
         do iele=1,nelec
         do iv=nvini,nvmax
         do j=jini,jmax
            ekinfin=etotS2(ie)-ediat(iv,j,iele)
            if(ekinfin.gt.0.d0)then
               pfin=dsqrt(ekinfin*xm2reac*2.d0/hbr/hbr)
            else
               pfin=0.d0
            endif
            if(paqini.lt.1.d-15)then
               S2factor(ie,iv,j,iele)=0.d0
            else
               S2factor(ie,iv,j,iele)=
     &              hbr*hbr*pfin*pini/(2.d0*paqini)/xm2reac/xmasa0
            endif
         enddo
         enddo
         enddo
      

** total
         xnorm=xnorm+paqini
         if(paqini.lt.1.d-15)then
            reacfct(ie)=0.d0
         else
            reacfct(ie)=hbr*hbr*pini/paqini/xmasa0/xm1reac/4.d0/pi/pi
         endif
         
         if(idproc.eq.0)write(43,*) etotS2(ie)/conve1/8065.5d0,paqini
     &         ,(etotS2(ie)-ediat(nvref,jref,ielecref))/conve1/8065.5d0
      enddo
      if(idproc.eq.0)close(43)

      xnorm=xnorm*estep
      if(idproc.eq.0)write(6,*)' Norm ini. wvp in energy= ',xnorm
      if(idproc.eq.0)call flush(6)

      
      return
      end subroutine Smatini_col
!--------------------------------------------------
      subroutine Smatini_photo
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2

      implicit none

      include "mpif.h"
      integer, parameter :: npunt=1000
      real*8 :: rgaussr(npunt)
      real*8 :: ekinfin,estep,etotmax,etotmin,pfin,ahgauss,arg,ekinini
      real*8 :: df,dg,f,g
      integer :: ie,iele,iv,j

      photonorm=1.d0
*> Emin,Emax kinetic energies to determine interval for S2 elements evaluation
     
       etotmin=ekinmin
       etotmax=ekinmax
       estep=(etotmax-etotmin)/dble(netot-1)


         if(idproc.eq.0)then
                  write(6,*)' Total ener. for S^2 in products (eV): '
     &             ,etotmin/conve1/8065.5d0,etotmax/conve1/8065.5d0
         endif
               
         do ie=1,netot
            etotS2(ie)=etotmin+dble(ie-1)*estep
!pru            reacfct(ie)=0.5d0/pi/xm1reac
            reacfct(ie)=0.25d0/pi/xm1reac
** reactants
            do iele=1,nelec
            do iv=nvini,nvmax
            do j=jini,jmax

               ekinfin=etotS2(ie)-ediat(iv,j,iele)
               if(ekinfin.gt.0.d0)then
                  pfin=dsqrt(ekinfin*xm2reac*2.d0/hbr/hbr)
               else
                  pfin=0.d0
               endif
!pru               S2factor(ie,iv,j,iele)=2.d0*pi*pfin/xm2reac
               S2factor(ie,iv,j,iele)=4.d0*pfin/(pi*xm2reac)
            enddo
            enddo
            enddo
         enddo
c      endif

      return
      end subroutine Smatini_photo
      
!--------------------------------------------------

      subroutine totalflux_k(ikcheb,kminelastic)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      implicit none
      include "mpif.h"
      complex*16 :: zzz,zexpo
      integer :: i,icanp,ielec,iom,iangp,ir,ir1,ir2,ican,iang,itotp
      integer :: iflux,ie
      integer :: kminelastic,ikcheb
      
!     evaluates total flux for r1=r1flux
!      add the contribution at the present k Chevishev iteration

      zfR1flux(:,:,:)=dcmplx(0.d0,0.d0)
      if(npun1.gt.1)then

         zzz=dcmplx(0.d0,hbr/(xm1red*ah1))
         do i=1,ntotproc(idproc)
            call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
            zfR1flux(ir2,iangp,icanp)=zfR1flux(ir2,iangp,icanp)
     &                      +zAR1flux(ir1)*dcmplx(rpaqproc(i),0.d0)
         enddo
         do icanp=1,ncanproc(idproc)
            ican=ibasproc(icanp,idproc)

            ielec=nelebas(ican)
            do iangp=1,nanguproc
               iang=indangreal(iangp,idproc)
               zR(:)=dcmplx(0.d0,0.d0)
               zdR(:)=dcmplx(0.d0,0.d0)
           
               do ir=1,npunreal(iang)
                  ir1=indr1(ir,iang)

                  if(ir1.eq.ir1flux(iang))then
                     ir2=indr2(ir,iang)
                     itotp=indtotproc(ir,icanp,iangp,idproc)
                     zR(ir2)=dcmplx(rpaqproc(itotp),0.d0)
                     zdR(ir2)=zfR1flux(ir2,iangp,icanp)
                  endif
               enddo  !  loop in ir

               do iflux=1,nfluxproc(idproc)
                  if(indangflux(iflux).eq.iang)then
                  if(indcanflux(iflux).eq.ican)then
                     ir2=indr2flux(iflux)
                     do ie=1,netot
                        zexpo=zCkcheby(iE)
                        if(iphoto.eq.0.and.iprod.eq.1)then
                     if(ikcheb.lt.kminelastic)zexpo=dcmplx(0.d0,0.d0)
                        endif
                        zCR(ie,iflux)=zCR(ie,iflux)+zexpo*zR(ir2)
                        zCdR(ie,iflux)=zCdR(ie,iflux)+zexpo*zdR(ir2)
                     enddo
                  endif
                  endif
               enddo
           
            enddo !  loop in iang
         enddo  ! loop in icanp

      endif ! npun1 > 1
      
      return
      end subroutine  totalflux_k
!--------------------------------------------------

      subroutine Cvjflux_k(ikcheb,kminelastic)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      implicit none
      include "mpif.h"
      complex*16 :: zzz,zexpo,zp,zder,zfluxvj
      real*8 :: xj,rori,xx,rpopu
      integer :: i,icanp,ielec,iom,iangp,ir,ir1,ir2,ican,iang,itotp
      integer :: maxnv,nsend,iv,j,m,ifail,ie,ierr
      integer :: kminelastic,ikcheb
      
!     evaluates flux on the individual v,j,e states in the 01 + 2 Jacobi coordinates
!      add the contribution at the present k Chevishev iteration

      rpopi(:,:,:)=0.d0
      funreac(:,:,:)=0.d0
      cvjproc(:,:,:)=0.d0
      cvj(:,:,:)=0.d0
      do itotp=1,ntotproc(idproc)
         call indiproc(itotp,icanp,ielec,iom,iangp,ir,ir1,ir2)
         if(ir2.eq.ir2abs)then
            funreac(iangp,ir1,icanp)=rpaqproc(itotp)/dsqrt(ah2)
         endif
      enddo

*     Projection of the wvp on the reactant eigenstates

      do icanp=1,ncanproc(idproc)
         ican=ibasproc(icanp,idproc)
         ielec=nelebas(ican)
         iom=iombas(ican)
         do ir1=1,npun1
            do j=j00(ican),jmax,inc
               xj=0.d0
               do iangp=1,nanguproc
                  iang=indangreal(iangp,idproc)
                  rori=funreac(iangp,ir1,icanp)
                  xj=xj+rori*Djmmp(iang,j,ielec,iom)
               enddo
  
               maxnv=min0(nvmax,noBCstates(j,ielec))
               do iv=nvini,maxnv
                  rpopi(iv,j,ican)=rpopi(iv,j,ican)
     &                +xj*fd(ir1,iv,j,ielec)
               enddo
            enddo ! j
         enddo ! ir1
      enddo  ! icanp

* flux on each individual channel

      do icanp=1,ncanproc(idproc)
         ican=ibasproc(icanp,idproc)
         ielec=nelebas(ican)
         do j=j00(ican),jmax,inc
            maxnv=min0(nvmax,noBCstates(j,ielec))
            do iv=nvini,maxnv
               Cvjproc(iv,j,ican)=rpopi(iv,j,ican)
            enddo
         enddo
      enddo

      nsend=(nvmax-nvini+1)*(jmax-jini+1)*ncanmax
      call MPI_REDUCE(Cvjproc,Cvj,nsend,MPI_REAL8,MPI_SUM
     &     ,0,MPI_COMM_WORLD,ierr)

      
***>> and performing FFT for Balint-Kurti method

      if(idproc.eq.0)then
         do ican=1,ncan
            ielec=nelebas(ican)
            do j=j00(ican),jmax,inc
               maxnv=min(nvmax,noBCstates(j,ielec))
               do iv=nvini,maxnv

                  ifail=0
                  if(ican.eq.icanref.and.iv.eq.nvref.and.j.eq.jref)then
                     if(iprod.eq.1)then
                  
                     else
                        if(iphoto.eq.0.and.ikcheb.lt.kminelastic)ifail=1
                     endif
                  endif
 
                  if(ifail.eq.0)then
                     zzz=dcmplx(Cvj(iv,j,ican),0.d0)
                     do ie=1,netot
                        zS2(ie,iv,j,ican)=zS2(ie,iv,j,ican)
     &                  +zzz*zCkcheby(iE)
                     enddo
                  endif

               enddo  ! iv
            enddo  ! j
         enddo  ! ican
      endif ! idproc = 0
       
      return
      end subroutine Cvjflux_k
!--------------------------------------------------

      subroutine totalflux
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      implicit none
      include "mpif.h"
      integer :: ie,iflux,ierr
      
*  total flux to products

      do ie=1,netot
         S2prod(ie)=0.d0          
         S2prodtot(ie)=0.d0          
         do iflux=1,nfluxproc(idproc)
            S2prod(ie)=S2prod(ie)
     &                +dimag(dconjg(zCR(ie,iflux))*zCdR(ie,iflux))
         enddo
         
         S2prod(ie)=S2prod(ie)*reacfct(ie)/ah1
      enddo

      call MPI_REDUCE(S2prod,S2prodtot,netot,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
      return
      end subroutine  totalflux
!--------------------------------------------------
!=======================================================
      end module mod_flux_01y2
