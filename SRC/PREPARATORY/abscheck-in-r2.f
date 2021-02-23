      implicit real*8(a-h,o-y)
      implicit complex*16(z)
      character*30 :: name,system,frpini,fRgini,fcoefini

      parameter(nexact=10,rpeq=0.7,npunaux=10000)

      real*8,allocatable,dimension(:) :: rpaq,rHpaq
      real*8,allocatable,dimension(:) :: rpaqexact,rHpaqexact
      real*8,allocatable,dimension(:) :: rflanz,rflanzexact
      real*8,allocatable,dimension(:) :: pot,absfR2,absfR2exact
      real*8,allocatable,dimension(:,:) ::  potmat
      real*8,allocatable,dimension(:) :: sigdiat,sigatom
      real*8,allocatable,dimension(:) :: pr2,p2r2
      real*8,allocatable,dimension(:) :: pr2exact,p2r2exact
      real*8,allocatable,dimension(:) :: flux0,flux,fluxexact,eflux
      complex*16,allocatable,dimension(:) :: zCkcheby,zS2,zS2exact
      integer,allocatable,dimension(:) :: iomdiat,iomatom
*********************************************************
      namelist /inputgridbase/npun1,rmis1,rfin1,npun1min
     &                       ,npun2,rmis2,rfin2
     &                       ,nangu,nangplot
     &     ,Jtot,iparity,inc,nelecmax,iommin,iommax,j0
     &     ,jini,jmax,nvini,nvmax
     &     ,nvref,jref,iomref,ielecref
*********************************************************
      namelist /inputabsorcion/absr1,absalp1,n1expo,ninterv1
     &                        ,absr2,absalp2,n2expo,ninterv2     
*********************************************************
      namelist /inputpotmass/system,xm1,xm0,xm2
     &                      ,VcutmaxeV,radcutmaxeV,rotcutmaxeV
*********************************************************
      namelist /inputcol/Rcolini_ang,ecol_ev,deltaE_ev
*********************************************************
      namelist /inputtime/ntimes, nloop,kminelastic
*********************************************************
      namelist /inputflux/r1flux_ang,netot,ekinmin_eV,ekinmax_eV
     &                   ,ncontfile
************************************************************


**> constants

      zero = dcmplx(0.d0,0.d0)
      zeye = dcmplx(0.d0,1.d0)
      pi = dacos(-1.d0)
      conve1 = 1.197d-4
      CONVL=.52917726D0
      CONVE=4.55633538D-6
      CONVM=.182288853D4
      hbr=convl/dsqrt(convm*conve/conve1)
      eV2cm = 8065.5d0
      
!     reading input data
      
      npun1min=32
      iomref=0
      ielecref=1
      iommin=0
      
      write(6,*)
      write(6,*)'  grid and basis data'
      write(6,*)'  -------------------'
         open(10,file='input.dat',status='old')
         read(10,nml = inputgridbase)
         write(6,nml = inputgridbase)
         call flush(6)
         close(10)
      write(6,*)
      write(6,*)'  absorption data'
      write(6,*)'  -------------------'
         open(10,file='input.dat',status='old')
         read(10,nml = inputabsorcion)
         write(6,nml = inputabsorcion)
         call flush(6)
         close(10)
      write(6,'(80("-"),/,10x
     &      ,"Mass and pot determination for 01+2= ",a50
     &      ,/80("-"),/)')system
         open(10,file='input.dat',status='old')
         read(10,nml = inputpotmass)
         write(6,nml = inputpotmass)
         close(10)
      write(6,'(40("="),/,10x,"colini_mod",/,40("="))')
      write(6,*)
      write(6,*)'  colini data'
      write(6,*)'  -----------'
      open(10,file='input.dat',status='old')
         read(10,nml = inputcol)
         write(6,nml = inputcol)
         call flush(6)
      close(10)
      Rcolini=Rcolini_ang
      ecol=ecol_ev*ev2cm*conve1
      deltaE=deltaE_ev*ev2cm*conve1
      write(6,*)
      write(6,*)'  propagation data'
      write(6,*)'  ----------------'
      open(10,file='input.dat',status='old')
         read(10,nml = inputtime)
         write(6,nml = inputtime)
         call flush(6)
      close(10)
      write(6,*)
      write(6,*)'  flux data'
      write(6,*)'  ---------'
      open(10,file='input.dat',status='old')
         read(10,nml = inputflux)
         write(6,nml = inputflux)
         call flush(6)
      close(10)

      r1flux=r1flux_ang
      ekinmin=ekinmin_eV*ev2cm*conve1
      ekinmax=ekinmax_eV*ev2cm*conve1

      nelec=nelecmax
!**>> masses
!* reduced masses (amu are the mass units using zots, in this program

      xmtot=xm0+xm1+xm2

* Reactant Jacobi
      xm1reac = (xm0*xm1)/(xm0+xm1)
      xm2reac = (xm2*(xm0+xm1))/xmtot

**>> maximum energy to cut potential
     
      vcutmax=vcutmaxeV*ev2cm*conve1
      radcutmax=radcutmaxeV*ev2cm*conve1
      rotcutmax=rotcutmaxeV*ev2cm*conve1

!! allocating matrices

      allocate(rpaq(npun2),rHpaq(npun2),absfR2(npun2)
     &       ,pot(npun2*nexact),absfR2exact(npun2*nexact)       
     &       ,rpaqexact(npun2*nexact),rHpaqexact(npun2*nexact)
     &       ,rflanz(npun2),rflanzexact(npun2*nexact)
     &       ,potmat(nelecmax,nelecmax),sigdiat(nelecmax)
     &       ,sigatom(nelecmax),iomdiat(nelecmax),iomatom(nelecmax)
     &       ,pr2(npun2),p2r2(npun2)
     &       ,flux0(netot),flux(netot),fluxexact(netot),eflux(netot)
     &       ,zCkcheby(netot),zS2(netot),zS2exact(netot)
     &       ,pr2exact(npun2*nexact),p2r2exact(npun2*nexact)
     &      ,stat=ierror)
      if(ierror.ne.0)then
        write(*,*)" error in initmem for zpaq and zHpaq, etc "
        stop
      endif

* Potential and barrier

      ah2 = (rfin2-rmis2)/dble(npun2-1)
      n2=npun2*nexact
      call setxbcpotele(iomdiat,iomatom,sigdiat,sigatom,nelec,nelecmax)
      r1=rpeq
      r1au=rpeq/convl
      ctet=0.d0
      vmin=1.d10
      open(10,file='potR2',status='unknown')
      write(6,*)' nelec,nelecmax,ielecref= ', nelec,nelecmax,ielecref
      do ir2=1,npun2*nexact
         r2=rmis2+dble(ir2-1)*ah2
         r2au=r2/convl
         call potelebond(r1au,r2au,ctet,potmat,nelec,nelecmax)

         pot(ir2)=potmat(ielecref,ielecref)/conve*conve1

         xl2=dble(Jtot*(Jtot+1))+dble(jref*(jref+1))
     &          -2.d0*dble(iomref*iomref)

         rot=xl2*hbr*hbr*0.5d0/(xm2reac*r2*r2)
         vinfty=pot(ir2)
         pot(ir2)=pot(ir2)+rot
         if(pot(ir2).gt.vcutmax+rotcutmax)pot(ir2)=vcutmax+rotcutmax
      enddo
      do ir2=1,npun2*nexact
         pot(ir2)=pot(ir2)-vinfty
         r2=rmis2+dble(ir2-1)*ah2
         if(pot(ir2).lt.vmin)vmin=pot(ir2)
         write(10,*)r2,pot(ir2)/conve1/8065.5d0,rot/conve1/8065.5d0
      enddo
      close(10)

* Absorption function

      ir2abs=npun2
      do ir2=1,npun2
         r2=rmis2+dble(ir2-1)*ah2
         f2=1.d0
         potabs2=0.d0
         if(r2.gt.absr2)then
            if(ir2.lt.ir2abs)ir2abs=ir2
            delta=(r2-absr2)/dble(ninterv2)
            deltan=(delta)**n2expo
            expo = absalp2*deltan
            f2 = dexp(-expo)
         endif
         absfr2(ir2)=f2
      enddo

      ir2absexact=n2
      do ir2=1,n2
         r2=rmis2+dble(ir2-1)*ah2
         f2=1.d0
         potabs2=0.d0
         if(r2.gt.rfin2)then
            if(ir2.lt.ir2absexact)ir2absexact=ir2
            delta=(r2-rfin2)/dble(ninterv2*nexact)
            deltan=(delta)**n2expo
            expo = absalp2*deltan/(dble(nexact*n2expo))
            f2 = dexp(-expo)
         endif
         absfr2exact(ir2)=f2
      enddo

* radial kinetic initialization

      box=dble(npun2-1)*ah2
      call sinmom(box,npun2,npun2,xm2reac,hbr,pr2,p2r2)
      xk2max=0.d0
      do ir2=1,npun2
         if(dabs(p2r2(ir2)).gt.xk2max)xk2max=dabs(p2r2(ir2))
      enddo

      n2=npun2*nexact
      box=dble(n2-1)*ah2
      call sinmom(box,n2,n2,xm2reac,hbr,pr2exact,p2r2exact)
      do ir2=1,n2
         if(dabs(p2r2exact(ir2)).gt.xk2max)xk2max=dabs(p2r2exact(ir2))
      enddo

* parameters for Chebyshev propagation

      emax=vcutmax+rotcutmax+xk2max
      emin=vmin
         write(6,*)'     Emin= ',vmin/conve1/8065.5d0
     &               ,'  Emax= ',emax/conve1/8065.5d0
         write(6,*)
      delta=emax-vmin
      delta2=delta*0.5d0
      emindlt=vmin+delta2
      write(6,*)' Energy scaling: delta2,emindlt=',delta2,emindlt

* initial wvp

      r2col = rcolini
      xk0=dsqrt(2.d0*xm2reac*ecol)/hbr
      alpha0=dsqrt(2.d0*ecol/xm2reac)*hbr/deltae
      factor=(2.d0/pi/alpha0/alpha0)**(0.25d0)
      DET=DBLE(Jtot*(Jtot+1)+jref*(jref+1)-2*iomref*iomref)
      IL=0.5D0*(-1.D0+DSQRT(1.D0+4*DET))+0.5D0
      PEPE=dble(il*(il+1))
      xxx=0.d0
      do ir2=1,npun2*nexact
         r2=rmis2+dble(ir2-1)*ah2
         call  rgauscolini(xr2,r2,r2col,alpha0,xk0,factor,il)
         rpaqexact(ir2)=xr2
         if(ir2.le.npun2)rpaq(ir2)=xr2
       
         xxx=xxx+xr2*xr2
      enddo
      xxx=xxx*ah2
      call Tr2sinFFT(rpaq,rHpaq,p2r2,npun2)

      it=0
         write(name,"('rpaq.k',i6.6)")it
         open(10,file=name,status='unknown')
         do ir2=1,npun2
            r2=rmis2+dble(ir2-1)*ah2
            write(10,*)
     &         r2,rpaq(ir2),rpaqexact(ir2),rHpaq(ir2)
         enddo
         close(10)

* flux initialization

      destep=(ekinmax-ekinmin)/dble(netot-1)
      do ie=1,netot
         e=ekinmin+dble(ie-1)*destep
         eflux(ie)=e
                  es=(E-emindlt)/delta2
                  expo=-dble(it)*dacos(Es)
                  write(69,*)e,es,expo,emindlt,delta2
         zS2(ie)=dcmplx(0.d0,0.d0)
         zS2exact(ie)=dcmplx(0.d0,0.d0)
         ekinini=eflux(ie)
         if(ekinini.gt.0.d0)then
            pini=dsqrt(ekinini*xm2reac*2.d0/hbr/hbr)
            zfft=0.d0
            ah2aux=(rfin2-rmis2)/dble(npunaux-1)
            do ir2=1,npunaux
               r=rmis2+dble(ir2-1)*ah2aux
               arg=r*pini
               CALL BESPH2(F,DF,G,DG,PEPE,ARG,KEY,0)
               zexpo=dcmplx(-g,f)
         call  rgauscolini(xr2,r,r2col,alpha0,xk0,factor,il)
               zfft=zfft+dcmplx(xr2,0.d0)*zexpo
            enddo
            zfft=zfft*ah2aux/2.d0/pi
            flux0(ie)=dreal(zfft*dconjg(zfft))
         endif
      enddo

* initializing propagation

      do ir2=1,npun2
        rflanz(ir2)=rpaq(ir2)
      enddo
      call Tr2sinFFT(rpaq,rHpaq,p2r2,npun2)
      do ir2=1,npun2
          rHpaq(ir2)=rHpaq(ir2)  + pot(ir2)*rpaq(ir2)           
      enddo
      do ir2=1,npun2
          rpaq(ir2)= (rHpaq(ir2)-emindlt*rflanz(ir2))/delta2  
      enddo
**
      do ir2=1,n2
        rflanzexact(ir2)=rpaqexact(ir2)
      enddo
      call Tr2sinFFT(rpaqexact,rHpaqexact,p2r2exact,n2)
      do ir2=1,n2
          rHpaqexact(ir2)=rHpaqexact(ir2)  + pot(ir2)*rpaqexact(ir2)           
      enddo
      do ir2=1,n2
          rpaqexact(ir2)= (rHpaqexact(ir2)
     &           -emindlt*rflanzexact(ir2))/delta2  
      enddo


* main loop in time

      it=0
      do iloop=1,nloop
         do i=1,ntimes
            it=it+1
* wvp

            call Tr2sinFFT(rpaq,rHpaq,p2r2,npun2)
         
            xnorm=0.0d0
            do ir2=1,npun2
               rHpaq(ir2)=rHpaq(ir2)  + pot(ir2)*rpaq(ir2)           
 
               fabs=absfr2(ir2)
               rHnorm=(rHpaq(ir2)-emindlt*rpaq(ir2))/delta2
               xxx=(2.d0*rHnorm-fabs*rflanz(ir2))*fabs
               rflanz(ir2)=rpaq(ir2)
               rpaq(ir2)=xxx
               xnorm=xnorm+rpaq(ir2)*rpaq(ir2)
            enddo
* wvp exact

            n2=npun2*nexact
            call Tr2sinFFT(rpaqexact,rHpaqexact,p2r2exact,n2)
    
            xnorm2=0.d0
            do ir2=1,n2
               rHpaqexact(ir2)=rHpaqexact(ir2)+ pot(ir2)*rpaqexact(ir2)           
               
               fabs=absfR2exact(ir2)
               rHnorm=(rHpaqexact(ir2)-emindlt*rpaqexact(ir2))/delta2
               xxx=(2.d0*rHnorm-fabs*rflanzexact(ir2))*fabs
               rflanzexact(ir2)=rpaqexact(ir2)
               rpaqexact(ir2)=xxx
               xnorm2=xnorm2+rpaqexact(ir2)*rpaqexact(ir2)
            enddo

            write(6,*)it,xnorm*ah2,xnorm2*ah2
* cheby coefficients and fluxes

            zCvj=dcmplx(rpaq(ir2abs),0.d0)
            zCvjexact=dcmplx(rpaqexact(ir2abs),0.d0)

            if(it.gt.kminelastic)then
               d=2.d0
               if(it.eq.0)d=1.d0
               zfactor=dcmplx(d*hbr/delta2,0.d0)
               do ie=1,nEtot
                  E=eflux(ie)
                  Es=(E-emindlt)/delta2
                  expo=-dble(it)*dacos(Es)
                  zexpo=cdexp(dcmplx(0.d0,expo))
                  deno=1.d0-Es*Es
                  deno=dsqrt(deno)
                  zCkcheby(iE)=zfactor*zexpo/dcmplx(deno,0.d0)

                  zS2(ie)=zS2(ie)+zCvj*zCkcheby(iE) 
                  zS2exact(ie)=zS2exact(ie)+zCvjexact*zCkcheby(iE)
               enddo
            endif

         enddo
*printing wavepackets

         write(name,"('rpaq.k',i6.6)")it
         open(10,file=name,status='unknown')
         do ir2=1,n2
            r2=rmis2+dble(ir2-1)*ah2
c            if(r2.lt.absr2)then
             f=0.d0
             if(ir2.le.npun2)f=rpaq(ir2)
            write(10,*)r2,f,rpaqexact(ir2)
c            endif
         enddo
         close(10)

         write(name,"('flux.k',i6.6)")it
         open(10,file=name,status='unknown')
         do ie=1,netot
            e=eflux(ie)
            pini=dsqrt(e*xm2reac*2.d0/hbr/hbr)
            S2=dreal(zS2(ie)*dconjg(zS2(ie)))*0.25d0/(pi*pi)
            S2exact=dreal(zS2exact(ie)*dconjg(zS2exact(ie)))
     &            *0.25d0/(pi*pi)

            write(10,*)e/conve1/8065.5d0
     &             ,flux0(ie)
     &             ,S2*((hbr*pini/xm2reac)**2)
     &             ,S2exact*((hbr*pini/xm2reac)**2)
     &             ,((hbr*pini/xm2reac)**2)

         enddo
         close(10)
      enddo

      stop
      end

**************************** TR2sinFFT ******************************

      subroutine Tr2sinFFT(rpaq,rHpaq,p2,npun)
      implicit real*8(a-h,o-y)
      implicit complex*16(z)
  
      INTEGER FFTW_R2HC
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_DFT_R2HC_ICKY
      PARAMETER (FFTW_DFT_R2HC_ICKY=512)
      INTEGER FFTW_NONTHREADED_ICKY
      PARAMETER (FFTW_NONTHREADED_ICKY=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)
      integer*8 plan,flags
      INTEGER req   ! , status(MPI_STATUS_SIZE)
 
       dimension rpaq(npun),rHpaq(npun),p2(npun)

      flags=FFTW_ESTIMATE


      n1=npun
     
      divi=0.5d0/dble(n1-1)
             do ir1=1,n1
                rHpaq(ir1)=rpaq(ir1)
             enddo
     

            call dfftw_plan_r2r_1d(plan,n1,rHpaq,rHpaq
     &                    ,FFTW_RODFT00
     &                     ,flags)
             call dfftw_execute(plan)

             call dfftw_destroy_plan(plan)

       
             do ir1=1,n1
                rHpaq(ir1)=rHpaq(ir1)*p2(ir1)*divi
             enddo
             call dfftw_plan_r2r_1d(plan,n1,rHpaq,rHpaq
     &                    ,FFTW_RODFT00
     &                     ,flags)
             call dfftw_execute(plan)
             call dfftw_destroy_plan(plan)

    
       return
       end
