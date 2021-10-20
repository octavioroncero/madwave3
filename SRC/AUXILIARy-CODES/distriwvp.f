      program distriparawvp  ! version v6  in progress
*
* This program extract the state-2-state S^2 matrix elements
*         and writes them in distriXXXX files  
*     to be used to calculate Cumulative Reaction probabilities (CRP)
*                             integral cross sections (ICS)
*                             rates
*        and the "Cvjprod" previously generated by madwave3 propagations
*

      implicit real*8(a-h,o-y)
      implicit complex*16(z)
      character*30 name,viejo,posicion,system,frpini,fRgini,fcoefini
*
* energies in eV (input/output)
*
      parameter (nmulgrid=16,ntreverse=1000,dt=0.05) ! data for backward propagation

*********************************************************
      namelist /inputtime/ntimes, nloop,kminelastic
      namelist /inputcol/Rcolini_ang,ecol_ev,deltaE_ev
      namelist /inputflux/r1flux_ang,netot,ekinmin_eV,ekinmax_eV
     &                   ,ncontfile
      namelist /inputgridbase/npun1,rmis1,rfin1,npun1min
     &                       ,npun2,rmis2,rfin2
     &                       ,nangu,nangplot
     &     ,Jtot,iparity,inc,nelecmax,iommin,iommax,j0
     &     ,jini,jmax,nvini,nvmax
     &                       ,nvref,jref,iomref,ielecref
      namelist /inputprod/iprod
     &                   ,nviniprod,nvmaxprod
     &                   ,jiniprod,jmaxprod
     &     ,iomminprod,iommaxprod
     &     ,Rbalinprod,n2prod0,n2prod1,nangproj0,nangproj1
      namelist /inputpotmass/system,xm1,xm0,xm2
     &                      ,VcutmaxeV,radcutmaxeV,rotcutmaxeV
*********************************************************
*Dimensions
      parameter (npunt=1024,npuntott=npunt*nmulgrid)
      dimension zgaussr(npunt)
      dimension zgaussr0(npuntott),zgaussr1(npuntott)
      dimension pr(npuntott),p2r(npuntott),vcentri(npuntott)
      dimension work(2*npuntott)

* preparing dimensions depending on Jacobi coordinates used in propagation (iprod=1 -products- or iprod=2 reactants) 

      real*8,allocatable,dimension(:,:,:) :: Cvjprod
      real*8, allocatable, dimension(:,:,:,:,:) :: prodCvj
      real*8,allocatable,dimension(:,:,:,:) :: S2pro,Cvj
      real*8,allocatable,dimension(:,:,:) :: S2prodfac,ediat,ediatprod
      real*8,allocatable,dimension(:) :: vibprod,rotprod,vibrot
      complex*16,allocatable,dimension(:,:,:,:) :: zS2prod
      integer :: nloopreal,ifail

**>> constants

      zero = dcmplx(0.d0,0.d0)
      zeye = dcmplx(0.d0,1.d0)
      pi = dacos(-1.d0)
      conve1 = 1.197d-4
      hbr = 0.063533625d0
      CONVL=.52917726D0
      CONVE=4.55633538D-6
      CONVM=.182288853D4
      nelec=nelecmax
      ceVcm=1.23984245d-4

      
!     reading data
         iommax=Jtot
         open(10,file='input.dat',status='old')
         read(10,nml = inputgridbase)
         close(10)
         iommaxprod=Jtot
         open(10,file='input.dat',status='old')
         read(10,nml = inputtime)
         close(10)
         open(10,file='input.dat',status='old')
         read(10,nml = inputcol)
         close(10)
         open(10,file='input.dat',status='old')
         read(10,nml = inputflux)
         close(10)
         open(10,file='input.dat',status='old')
         read(10,nml = inputprod)
         close(10)
         open(10,file='input.dat',status='old')
         read(10,nml = inputpotmass)
         close(10)

         write(6,*)'  --- input data ---'
         write(6,nml = inputgridbase)
         write(6,nml = inputtime)
         write(6,nml = inputcol)
         write(6,nml = inputflux)
         write(6,nml = inputprod)         
         write(6,nml = inputpotmass)
         write(6,*)'  --- end input data ---'

         emindistri=ekinmin_eV
         emaxdistri=ekinmax_eV
         nener=netot
         nloop0=nloop

**> reading initial input

      open(16,file='../pot/cont.pot',status='old')
      read(16,*)nelec
      do ie=1,nelec
         read(16,*)
      enddo
      read(16,*,err=11)vmin,ijklm,vmaxtot
      go to 22
 11    continue
          vmaxtot=vcutmaxeV
          close(16)
          open(16,file='../pot/cont.pot',status='old')
          read(16,*)nelec
          do ie=1,nelec
             read(16,*)
          enddo
          write(6,*)' old verion of mod_pot: setting vmaxtot=vcutmax'
          read(16,*)vmin
          
 22    continue
      close(16)
      vmintot=vmin
      vmaxtoteV=vmaxtot/(8065.5d0*conve1)
      emaxtot=vmaxtoteV+radcutmaxeV+3.d0*rotcutmaxeV
      emaxtot=emaxtot*8065.5d0*conve1
      delta2=0.5d0*(emaxtot-vmintot)
      emindlt=vmin+delta2

      write(6,*)' -- Energy interval  in the calculation (eV)'
      write(6,*)'    Emin= ',vmin/conve1/8065.5d0
     &         ,'     Emax= ',emaxtot/conve1/8065.5d0
      write(6,*)vmin/conve1/8065.5d0
      write(6,*)emindlt/conve1/8065.5d0,delta2/conve1/8065.5d0

**> preparing magnitudes for iprod=1 or 2

      xmtot=xm0+xm1+xm2
      if(iprod.eq.1)then              ! products Jacobi coordinates
         iom0=iommin
         iom1=iommax
         j00=jini
         j1=jmax
         nv0=nvini
         nv1=nvmax
         nelecprod=nelec
         xmasa=(xm1*(xm0+xm2))/xmtot
         xmasaprod=(xm2*(xm0+xm1))/xmtot
         allocate(
     &      ediat(nviniprod:nvmaxprod,jiniprod:jmaxprod,nelecmax)
     &     ,ediatprod(nvini:nvmax,jini:jmax,nelecmax)
     &          , stat=ierror)
         open(16,file='../func/bcwf',status='old')
         do ielec=1,nelecprod
         do j=jini,jmax
            read(16,*)iielec,jj,noBCstates
            do iv=nvini,noBCstates
               read(16,*)iiv,ediatprod(iv,j,ielec)
               do ir1=1,npun1
                  read(16,*)
               enddo
            enddo
         enddo
         enddo
         close(16)
         open(16,file='../func/prodwf',status='old')
         do j=jiniprod,jmaxprod
            read(16,*)jj,noBCprod,nn2prod1,nnelec
            do iv=nviniprod,noBCprod
               read(16,*)iiv,ediat(iv,j,1)
            
               do ir2=1,nn2prod1
                  read(16,*)
               enddo
            enddo
         enddo
         close(16)
         allocate(Cvj(nv0:nv1,j00:j1,iom0:iom1,nelecmax))
      elseif(iprod.eq.2)then          ! reactants Jacobi coordinates
         iom0=iomminprod
         iom1=iommaxprod
         j00=jiniprod
         j1=jmaxprod
         nv0=nviniprod
         nv1=nvmaxprod
         nelecprod=1
         allocate(Cvjprod(nv0:nv1,j00:j1,iom0:iom1))

         xmasa=(xm2*(xm0+xm1))/xmtot
         xmasaprod=(xm1*(xm0+xm2))/xmtot
         allocate(
     &      ediatprod(nviniprod:nvmaxprod,jiniprod:jmaxprod,nelecmax)
     &     ,ediat(nvini:nvmax,jini:jmax,nelecmax)
     &          , stat=ierror)
         open(16,file='../func/bcwf',status='old')
         do ielec=1,nelecmax
         do j=jini,jmax
            read(16,*)iielec,jj,noBCstates
            do iv=nvini,noBCstates
               read(16,*)iiv,ediat(iv,j,ielec)
               do ir1=1,npun1
                  read(16,*)
               enddo
            enddo
         enddo
         enddo
         close(16)
         open(16,file='../func/prodwf',status='old')
         do j=jiniprod,jmaxprod
            read(16,*)jj,noBCprod,nn2prod1,nnelec
            do iv=nviniprod,noBCprod
               read(16,*)iiv,ediatprod(iv,j,nelecprod)
            
               do ir2=1,nn2prod1
                  read(16,*)
               enddo
            enddo
         enddo
         close(16)
         write(6,*)' masses= ',xmasa,xmasaprod,xm1,xm0,xm2
      endif

      do ielec=1,nelec
      do iv=nvini,nvmax
         do j=jini,jmax
            write(26,*)ielec,iv,j,ediat(iv,j,ielec)/(conve1*8065.5d0)
         enddo
      enddo
      enddo
      do iv=nviniprod,nvmaxprod
         do j=jiniprod,jmaxprod
            ielec=1
          write(36,*)ielec,iv,j,ediatprod(iv,j,ielec)/(conve1*8065.5d0)
         enddo
      enddo
**> allocating matrices

      allocate(prodCvj(ntimes*nloop0,nv0:nv1,j00:j1,nelecmax,iom0:iom1)
     &          ,zS2prod(nv0:nv1,j00:j1,nelecmax,iom0:iom1)
     &          ,S2prodfac(nv0:nv1,j00:j1,nelecmax)
     &          ,S2pro(nv0:nv1,j00:j1,nelecmax,iom0:iom1)
     &          ,vibprod(nv0:nv1)
     &          ,rotprod(j00:j1)
     &          ,vibrot(j00:j1)
     &          , stat=ierror)

**> Reading Cvj coefficients

      kcheby=0
      ifail=1
      if(iprod.eq.1)then
         do iloop=1,nloop0
            write(name,'("Cvj.",i4.4)')iloop
            write(6,*)name
            open(16,file=name,status='old',form='unformatted',err=1)
            do ik=1,ntimes
               kcheby=kcheby+1
               read(16,err=2)it,Cvj
               do Iom=iom0,iom1
               do ielec=1,nelecmax
               do j=j00,j1
                  do iv=nv0,nv1
                  prodCvj(kcheby,iv,j,ielec,iom)=Cvj(iv,j,iom,ielec)
               enddo
               enddo
               enddo
               enddo
            enddo
           
            close(16)
         enddo  ! end iloop
      else
        do iloop=1,nloop0
            write(name,'("Cvjprod.",i4.4)')iloop
            write(6,*)name
            open(16,file=name,status='old',form='unformatted',err=1)
            do ik=1,ntimes
               kcheby=kcheby+1
               read(16,err=3)it,Cvjprod
               do Iom=iom0,iom1
               do ielec=1,nelecprod
               do j=j00,j1
               do iv=nv0,nv1
               prodCvj(kcheby,iv,j,ielec,iom)=Cvjprod(iv,j,iom)
               enddo
               enddo
               enddo
               enddo
            enddo
            
            close(16)
         enddo
      endif
      ifail=0
      nloopreal=nloop0
 1    continue
      if(ifail.eq.1)then
         nloopreal=iloop-1
      endif
      write(6,*)'Reading coefficients up to nloop = ',loopreal

**> Initial wavepacket

         photonorm=1.d0
         r2col = rcolini_ang
         ecol=ecol_ev*8065.5*conve1
         deltae=deltae_ev*8065.5*conve1
         ediatref=ediat(nvref,jref,ielecref)
         write(6,*)' ediatref= ',ediatref/conve1

         xk0=dsqrt(2.d0*xmasa*ecol)/hbr
         alpha0=dsqrt(2.d0*ecol/xmasa)*hbr/deltae
         factor=(2.d0/pi/alpha0/alpha0)**(0.25d0)
         ahgauss=(rfin2-rmis2)/dble(npunt-1)
            write(6,*)'  **  R2 initial conditions **'
            write(6,*)'      colision energy(eV)= ',ecol/conve1/8065.5d0
            write(6,*)'      energy width(eV)= ',deltae/conve1/8065.5d0
            write(6,*)'  r20 =',r2col,'   k= ',xk0,'   alpha= ',alpha0
            write(6,*)
            write(6,*)' redmass= ',xmasa
         DET=DBLE(Jtot*(Jtot+1)+jref*(jref+1)-2*iomref*iomref)
         IL=0.5D0*(-1.D0+DSQRT(1.D0+4.d0*DET))+0.5D0
         pepe=dble(il*(il+1))
         key=-1
         
         do ir2=1,npunt
            r2=rmis2+dble(ir2-1)*ahgauss
             call  rgauscolini(xr2,r2,r2col,alpha0,xk0,factor,il)
            zgaussr(ir2)=dcmplx(xr2/dsqrt(2.d0),0.d0)
            erot=hbr*hbr*pepe*0.5d0/(r2*r2*xmasa)
            write(27,*)r2,dreal(zgaussr(ir2)*dconjg(zgaussr(ir2))),erot
         enddo

**>> Main Loop in energy

      open(15,file='distriS2prod.elec',status='unknown')
      ifile=36
      ivfile=15
      iii=ifile
      iiiv=ivfile
      
      do ielec=1,nelecprod
      do ivprod=nv0,nv1
          iiiv=iiiv+1
          write(name,"('distriS2prod.v',i2.2,'.e',i1.1)")
     &                 ivprod,ielec
          open(iiiv,file=name,status='unknown')
!          do iomprod=iom0,iom1
!             iii=iii+1
!             write(name,"('distriS2prod.v',i2.2,'.Omg',i2.2,'.e',i1.1)")
!     &                 ivprod,iomprod,ielec
!             open(iii,file=name,status='unknown')
!         enddo
      enddo
      enddo
      emin=emindistri*8065.5d0*conve1
      emax=emaxdistri*8065.5d0*conve1
      delta=(emax-emin)/dble(nener-1)
      r=rcolini+1.d0
      erot=hbr*hbr*pepe*0.5d0/(r*r*xmasa)
      do ie=1,nener
         e=emin+dble(ie-1)*delta
         ekinini=e!-ediatref

**> initial flux
         pini=0.d0
         paqini=0.d0
         if(ekinini.gt.0.d0)then
            pini=dsqrt(ekinini*xmasa*2.d0/hbr/hbr)
            zfft=0.d0
            do ir2=1,npunt
               r=rmis2+dble(ir2-1)*ahgauss
               arg=r*pini
               
               CALL BESPH2(F,DF,G,DG,PEPE,ARG,KEY,0)
                  
               if(dabs(f).gt.1.d-20)then
                  zexpo=dcmplx(-g,f)
                  zfft=zfft+zgaussr(ir2)*zexpo
               endif
            enddo
            zfft=zfft*ahgauss/2.d0/pi
            paqini=dreal(zfft*dconjg(zfft))
         endif
         write(6,'(10(1x,e15.7))')e/(conve1*8065.5),paqini,pini
**> products properties
         do iele=1,nelec
         do iv=nv0,nv1
         do j=j00,j1
c            if(iprod.ne.1)then
               ekinfin=e-ediatprod(iv,j,iele)+ediatref
c            elseif(iprod.eq.1)then
c               ekinfin=e-ediat(iv,j,iele)
c            endif

            if(ekinfin.gt.0.d0)then
               pfin=dsqrt(ekinfin*xmasaprod*2.d0/hbr/hbr)
            else
               pfin=0.d0
            endif
            if(paqini.lt.1.d-15)then
               S2prodfac(iv,j,iele)=0.d0
            else
               S2prodfac(iv,j,iele)=
     &              hbr*hbr*pfin*pini/(2.d0*paqini)/xmasa/xmasaprod
            endif
         enddo
         enddo
         enddo

**> initialization

         do ielec=1,nelecprod
         do iomprod=iom0,iom1
         do j=j00,j1
         do iv=nv0,nv1
             zS2prod(iv,j,ielec,iomprod)=dcmplx(0.d0,0.d0)
         enddo
         enddo
         enddo
         enddo
**> sum in Chevishev iterations

         zfactor=dcmplx(2.d0*hbr/delta2,0.d0)
         do ikcheb=1,ntimes*nloopreal
            Es=(E-emindlt)/delta2
            expo=-dble(ikcheb)*dacos(Es)
            zexpo=cdexp(dcmplx(0.d0,expo))
            deno=1.d0-Es*Es
            deno=dsqrt(deno)
            zCkcheby=zfactor*zexpo/dcmplx(deno,0.d0)
            zexpo=zCkcheby
            do ielec=1,nelecprod
            do iomprod=iom0,iom1
            do j=j00,j1
            do iv=nv0,nv1
               xxx=prodCvj(ikcheb,iv,j,ielec,iomprod)
               
               zS2prod(iv,j,ielec,iomprod)=
     &                  zS2prod(iv,j,ielec,iomprod)
     &                  +dcmplx(xxx,0.d0)*zexpo
            enddo
            enddo
            enddo
            enddo

         enddo  ! end loop in cheb iterations   

**> Finishing
            S2no=0.d0
            iii=ifile
            iiiv=ivfile
            do iv=nv0,nv1
               vibprod(iv)=0.d0
            enddo
        
            do ielec=1,nelecprod
            do iv=nv0,nv1
               iiiv=iiiv+1
               do j=j00,j1
                  vibrot(j)=0.d0
               enddo
                     
               do iom=iom0,iom1
                  iii=iii+1
                  do j=j00,j1
                     rotprod(j)=0.d0
                  enddo
                  if(iabs(iom).le.j1)then
                     jjjj0=max0(iabs(iom),j00)
                     do j=jjjj0,j1
                        zzz=zS2prod(iv,j,ielec,iom)/2.d0/pi
                        Av=dreal(zzz*dconjg(zzz))
                        S2pro(iv,j,ielec,iom)=Av*S2prodfac(iv,j,ielec)
                        S2no=S2no+S2pro(iv,j,ielec,iom)
                        vibprod(iv)=vibprod(iv)+S2pro(iv,j,ielec,iom)
                        rotprod(j)=rotprod(j)+S2pro(iv,j,ielec,iom)
                        vibrot(j)=vibrot(j)+rotprod(j)      
                     enddo
                  endif
!                  write(iii,'(1000(1x,e15.7))') e/(conve1*8065.5d0)
!     &                  ,(rotprod(j)*photonorm,j=j00,j1) 
               enddo  ! iom
               write(iiiv,'(1000(1x,e15.7))') e/(conve1*8065.5d0)
     &                  ,(vibrot(j)*photonorm,j=j00,j1) 
 
            enddo  ! iv
            enddo  ! ielec
      
            write(15,"(1000(1x,e15.7))")e/(conve1*8065.5d0)
     &             ,S2no*photonorm
     &            ,(vibprod(iv)*photonorm,iv=nv0,nv1)
      enddo   ! end loop in energy

**> deallocation

      deallocate(prodCvj
     &          ,Cvjprod
     &          ,zS2prod
     &          ,S2prodfac
     &          ,S2pro
     &          ,vibprod
     &          ,rotprod
     &          ,vibrot
     &          , stat=ierror)
      stop

 2    continue
      write(6,*)' Cvj for iloop=',iloop,' not complete,ik= ',ik
      call flush(6)
      stop
 3    continue
      write(6,*)' Cvjprod for iloop=',iloop,' not complete,ik= ',ik
      call flush(6)
      stop
      
      end


****************
      subroutine propatras(zgaussr0,zgaussr1,work,vcentri,p2r,hbr,dt
     &                 ,rmis2,ahgauss,ntreverse,npuntot)

      implicit complex*16(z)
      implicit real*8(a-h,o-y)
      dimension zgaussr0(npuntot),zgaussr1(npuntot)
      dimension vcentri(npuntot),p2r(npuntot)
      dimension work(2*npuntot)

      conve1 = 1.197d-4

      xnorm=0.d0
      do ir=1,npuntot
         xnorm=xnorm+dreal(zgaussr0(ir)*dconjg(zgaussr0(ir)))
      enddo
      xnorm=xnorm*ahgauss

      open(33,file='paq0',status='unknown')
      do ir=1,npuntot
         r=rmis2+dble(ir-1)*ahgauss
         xxx=dreal(zgaussr0(ir)*dconjg(zgaussr0(ir)))
         if(xxx.lt.1.d-30)xxx=0.d0
         write(33,'(3(1x,e15.7))')r,xxx,vcentri(ir)/(conve1*8065.5d0)                                   
      enddo
      close(33)

      do it=1,ntreverse
* cdexp(+iV dt /2 hbr)
         do ir=1,npuntot
            arg=vcentri(ir)*dt*0.5d0/hbr
            zarg=dcmplx(0.d0,arg)
            zgaussr1(ir)=zgaussr0(ir)*cdexp(zarg)
         enddo
* cdexp(+iT dt / hbr)
          call four1(zgaussr1,work,npuntot,+1)
          do ir=1,npuntot
             arg=p2r(ir)*dt/hbr
             zarg=cdexp(dcmplx(0.d0,arg))
             zgaussr1(ir)=zgaussr1(ir)*zarg
          enddo
          call four1(zgaussr1,work,npuntot,-1)

* cdexp(+iV dt /2 hbr)
         do ir=1,npuntot
            arg=vcentri(ir)*dt*0.5d0/hbr
            zarg=dcmplx(0.d0,arg)
            zgaussr0(ir)=zgaussr1(ir)*cdexp(zarg)
         enddo
         
         xnorm=0.d0
         do ir=1,npuntot
            xnorm=xnorm+dreal(zgaussr0(ir)*dconjg(zgaussr0(ir)))
         enddo
         xnorm=xnorm*ahgauss
         write(6,*)it,xnorm
      enddo ! it
      open(33,file='paq1',status='unknown')
      do ir=1,npuntot
         r=rmis2+dble(ir-1)*ahgauss
         xxx=dreal(zgaussr0(ir)*dconjg(zgaussr0(ir)))
         if(xxx.lt.1.d-30)xxx=0.d0
         write(33,'(3(1x,e15.7))')r,xxx
     &                                   
      enddo
      close(33)

      return
      end

c******************************************************************

      SUBROUTINE four1(zpaq,data,nn,isign)
      implicit real*8(a-h,o-y)
      implicit complex*16(z)
      dimension zpaq(nn),data(2*nn)      
       
      n=2*nn

      k=0
      do i=1,n,2
         k=k+1
         data(i)=dreal(zpaq(k))
         data(i+1)=dimag(zpaq(k))
      enddo

      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif

      pp=1.d0
      if(isign.eq.-1)pp=dble(nn)
      k=0
      do i=1,n,2
         k=k+1
         zpaq(k)=dcmplx(data(i)/pp,data(i+1)/pp)
      enddo


      return
      END


