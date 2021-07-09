      program sigmaFromS2prod
      implicit none
!-------------------------------------------------------------------------------------------------------------
!     Calculates the reaction cross section
!     by adding all partial waves total reaction probabilities as
!      
!     sigma= pi/{2j+1 k_vj^2}  sum_J=0^Jtotmax P_J(E) (2J+1)
!
!     where the P_J(E) correspond to the 2nd column of S2prod files
!
!     Three approximations can be used used:
!     a) J-shifting:  only needs P_J=0 and a good estimate of Be
!     b) interpolation based on J-shifting: read NJcalc values of J
!        to be read,
!        and uses J-shifting interpolation por intermediate values
!     c) Exact: reading all J up to Jtotmax 
!
!     Input  look at a,b,c
!     -----
!
!    a)  values read in namelist/sigmaS2prod ( kept in input.dat files used for madwave3 propagations):
!        -Jtotmax :(maximum value of Jtot in the partial wave expansion)
!        -Ncalc: Number of P_J to be read
!        -nESigma,EminSigma_eV,EmaxSigma_Ev: energy grid used (EminSigma_eV and EmaxSigma_eV must be inside
!        ekinmin_eV,ekinmax_eV values defined in inputflux namelist)
! 
!    b)  The  particular J values for which P_J is read are defined in file: "CalculatedJ.dat"   read in namelist/sigmaFs2prod
!     kept in input.dat files used for madwave3 propagations.
!         values of J and recommened B_J rotational constants
!
!    c)  the S2prod files for each J to be read are in directories ../JXYZ for J=XYZ
!
!     Output files:
!     ------
!     sigma.dat:  with 2 columns with energy (in eV) and cross-section in Angstrom^2
!     PJ.dat file: with columns energy (in eV) and Jtotmax+1 columns with PJ(E) from J=0 to J=Jtotmax
!
!     Convert magnitudes internaly to atomic units
!
!---------------------------------------------------------------------------------------------------------------------------
      integer:: npun1,npun2,nangu,nangplot,Jtot,iparity,inc
     &     ,nelecmax,iommin,iommax,j0,jini,jmax,nvini,nvmax,nvref
     &     ,jref,iomref,ielecref,Jtotmax,nCalc,nESigma,netot
     &     ,ncontfile,npun1min,ierror,i,iE,iJ,iJ0,iJ1,iEE
     &     ,iprod,nviniprod,nvmaxprod,jiniprod,jmaxprod
     &     ,iomminprod,iommaxprod,n2prod0,n2prod1,nangproj0,nangproj1
     &     ,iom0max,iom
      double precision :: r1flux_ang
     &     ,ekinmin_eV,ekinmax_eV,rmis1,rfin1,rmis2,rfin2
     &     ,xm1,xm0,xm2,EminSigma_eV,EmaxSigma_Ev
     &     ,VcutmaxeV,radcutmaxeV,rotcutmaxeV,pi,conve1,hbr,conve,convm
     &     ,ceVcm,EminSigma,EmaxSigma,xmtot,xmasa,xK2,Rbalinprod,de
     &     ,fac,sigma,convl,E,E_eV,sigmaprod,xx
      character*50 :: system,name
      double precision, allocatable :: PJcalc(:,:,:),Ecalc(:,:,:)
      double precision, allocatable :: PJcalcProd(:,:,:),PJprod(:,:,:)
     &                               ,PJ(:,:,:),BeJcalc(:)
      integer, allocatable :: nEcalc(:,:),Jtotcalc(:)
      double precision, allocatable :: x(:),f(:,:)
      
      namelist /inputflux/r1flux_ang,netot,ekinmin_eV,ekinmax_eV
     &                   ,ncontfile
      namelist /inputgridbase/npun1,rmis1,rfin1,npun1min
     &                ,npun2,rmis2,rfin2,nangu,nangplot
     &              ,Jtot,iparity,inc,nelecmax,iommin,iommax,j0
     &               ,jini,jmax,nvini,nvmax,nvref,jref,iomref,ielecref
       
      namelist /inputprod/iprod
     &                   ,nviniprod,nvmaxprod
     &                   ,jiniprod,jmaxprod
     &     ,iomminprod,iommaxprod
     &     ,Rbalinprod,n2prod0,n2prod1,nangproj0,nangproj1
      namelist /inputpotmass/system,xm1,xm0,xm2
     &                      ,VcutmaxeV,radcutmaxeV,rotcutmaxeV
      namelist/sigmaS2prod/Jtotmax,nCalc
     &     ,nESigma,EminSigma_eV,EmaxSigma_Ev

      write(6,'(40("-"))')
      write(6,'(10x,"Program sigmaS2prod")')
      write(6,'(40("-"))')

!**>> constants

      pi = dacos(-1.d0)
      conve1 = 1.197d-4
      hbr = 0.063533625d0
      CONVL=.52917726D0
      CONVE=4.55633538D-6
      CONVM=.182288853D4
      ceVcm=1.23984245d-4

!**>>     reading data

      open(10,file='input.dat',status='old')
      read(10,nml = inputgridbase)
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
      open(10,file='input.dat',status='old')
      read(10,nml = sigmaS2prod)
      close(10)

      write(6,*)'  --- input data ---'
      write(6,nml = inputgridbase)
      write(6,nml = inputflux)
      write(6,nml = inputprod)
      write(6,nml = inputpotmass)
      write(6,nml = sigmaS2prod)
      write(6,*)'  --- end input data ---'
!**>> reduced mass
      xm0=xm0*convm
      xm1=xm1*convm
      xm2=xm2*convm
      xmtot=xm0+xm1+xm2
      if(iprod.eq.1)then
        xmasa=(xm1*(xm0+xm2))/xmtot
        write(6,'("--> Using product Jacobi coordinates",/
     &         ,10x,"with BC mases= ",2(F12.8,1x))')xm0/convm,xm2/convm
      else               
        xmasa=(xm2*(xm0+xm1))/xmtot
        write(6,'("--> Using reactant Jacobi coordinates",/
     &         ,10x,"with BC mases= ",2(F12.8,1x))')xm0/convm,xm1/convm 
      endif
!**>> energy grids
      if(EminSigma_eV.ge.ekinmin_eV.and.EmaxSigma_Ev.le.ekinmax_eV)then
         EminSigma=EminSigma_eV*conve/ceVcm
         EmaxSigma=EmaxSigma_eV*conve/ceVcm
      else
         write(6,*)' EminSigma_eV, EmaxSigma_eV out of the interval '
         write(6,*)'    defined in flux in ekinmin_eV,ekinmax_eV'
         call flush(6)
         stop
      endif

!**   >> allocating matrices

      allocate(PJcalc(netot,nCalc,-jref:jref)
     &     ,Ecalc(netot,nCalc,-jref:jref),nEcalc(nCalc,-jref:jref)
     &     ,PJ(nEsigma,0:Jtotmax,-jref:jref),Jtotcalc(ncalc)
     &     ,PJprod(nEsigma,0:Jtotmax,-jref:jref)
     &     ,PJcalcProd(netot,nCalc,-jref:jref)
     &     ,x(netot),f(netot,2),BeJcalc(ncalc)
     &     ,stat=ierror)

      PJ(:,:,:)=0.d0
      PJprod(:,:,:)=0.d0
!**   >>  Reading calculated PJ

      open(10,file="CalculatedJ.dat",status="old")
      do i=1,nCalc
         read(10,*)Jtotcalc(i),BeJcalc(i)
         BeJcalc(i)=BeJcalc(i)*conve/ceVcm
      enddo
      close(10)

      do i=1,nCalc
         iom0max=jref
         if(iom0max.gt.Jtotcalc(i))iom0max=Jtotcalc(i)
         
         do iom=-iom0max,iom0max

            if(iom.eq.0)then
               write(name,'("../Omg",i1,"/J",i3.3,"/S2prod")')
     &              iom,Jtotcalc(i)
            elseif(iom.gt.0)then
               write(name,'("../Omg",i1,"/p/J",i3.3,"/S2prod")')
     &              iabs(iom),Jtotcalc(i)
            elseif(iom.lt.0)then
               write(name,'("../Omg",i1,"/m/J",i3.3,"/S2prod")')
     &              iabs(iom),Jtotcalc(i)
            endif
            write(6,*)' Reading S2prod for J,Omg=',Jtotcalc(i),iom
     &                                              ,' in ',name
            open(10,file=name,status='old')
            iEE=0
            do iE=1,netot
               if(iprod.lt.2)then
                  read(10,*,end=1)Ecalc(ie,i,iom),PJcalc(ie,i,iom)
                  PJcalcProd(ie,i,iom)=0.d0
               else
                  read(10,*,end=1)Ecalc(ie,i,iom),PJcalc(ie,i,iom)
     &              ,xx,PJcalcProd(ie,i,iom)
               endif
               if(PJcalc(ie,i,iom).gt.1.d0)PJcalc(ie,i,iom)=1.d0
               if(PJcalcProd(ie,i,iom).gt.1.d0)PJcalcProd(ie,i,iom)=1.d0
               
               iEE=iEE+1
               Ecalc(ie,i,iom)=Ecalc(ie,i,iom)*conve/ceVcm
            enddo
 1          continue
            nEcalc(i,iom)=iEE
            close(10)
         enddo  ! iom 
      enddo  ! Jtotcalc

!**   >> Interpolating/extrapolating to get PJ(0,...Jtotmax)

      call JshiftInterPol(PJ,PJcalc,Ecalc,BeJCalc,x,f
     &     ,EminSigma,EmaxSigma
     &     ,nEcalc,Jtotcalc,nESigma,Jtotmax,nCalc,netot,jref)
      call JshiftInterPol(PJprod,PJcalcProd,Ecalc,BeJCalc,x,f
     &     ,EminSigma,EmaxSigma
     &     ,nEcalc,Jtotcalc,nESigma,Jtotmax,nCalc,netot,jref)

      dE=(EmaxSigma-EminSigma)/dble(nESigma)
      do iom=-jref,jref
         if(iom.eq.0)then
            write(name,'("PJOmg",i1,".dat")')iabs(iom)
         elseif(iom.gt.0)then
            write(name,'("PJOmg",i1,"p.dat")')iabs(iom)
         elseif(iom.lt.0)then
            write(name,'("PJOmg",i1,"m.dat")')iabs(iom)
         endif
         open(10,file=name,status='unknown')
         do ie=1,nESigma
            E=EminSigma+dble(ie-1)*dE
            E_eV=E*ceVcm/conve
            do Jtot=0,Jtotmax
               if(PJ(ie,Jtot,iom).lt.1.d-50)then
                       PJ(ie,Jtot,iom)=0.d0
               endif
            enddo
            write(10,'(1000(1x,e15.7))')E_eV
     &                 ,(PJ(ie,Jtot,iom),Jtot=0,Jtotmax)
         enddo
         close(10)
         if(iom.eq.0)then
            write(name,'("PprodJOmg",i1,".dat")')iabs(iom)
         elseif(iom.gt.0)then
            write(name,'("PprodJOmg",i1,"p.dat")')iabs(iom)
         elseif(iom.lt.0)then
            write(name,'("PprodJOmg",i1,"m.dat")')iabs(iom)
         endif
         open(10,file=name,status='unknown')
         do ie=1,nESigma
            E=EminSigma+dble(ie-1)*dE
            E_eV=E*ceVcm/conve
            do Jtot=0,Jtotmax
               if(PJprod(ie,Jtot,iom).lt.1.d-50)then
                       PJprod(ie,Jtot,iom)=0.d0
               endif
            enddo
            write(10,'(1000(1x,e15.7))')E_eV
     &                 ,(PJprod(ie,Jtot,iom),Jtot=0,Jtotmax)
         enddo
         close(10)
      enddo
! total cross section
      open(10,file='sigma.dat',status='unknown')
      write(10,*)'# Energy (eV)   cross-sections (Angstrom^2)'
      do ie=1,nESigma
         E=EminSigma+dble(ie-1)*dE
         E_eV=E*ceVcm/conve
         xK2=2.d0*xmasa*E        ! all in a.u. 
         fac=pi/(xk2*dble(2*jref+1))
         fac=fac*convl*convl     ! converting to Angstrom^{-2}
         sigma=0.d0
         sigmaprod=0.d0
         do Jtot=0,Jtotmax
            iom0max=jref
            if(iom0max.gt.Jtot)iom0max=Jtot
            do iom=-iom0max,iom0max
            
               sigma=sigma+dble(2*Jtot+1)*PJ(ie,Jtot,iom)
               sigmaprod=sigmaprod+dble(2*Jtot+1)*PJprod(ie,Jtot,iom)
            enddo
         enddo
         write(10,'(300(1x,e15.7))')E_eV,sigma*fac,sigmaprod*fac
      enddo
      close(10)
! helicity cross section
      if(jref.gt.0)then
         do iom=-jref,jref
            if(iom.eq.0)then
               write(name,"('sigma-Omg0.dat')")
            elseif(iom.gt.0)then
               write(name,"('sigma-Omg',i1,'p.dat')")iom
            elseif(iom.lt.0)then
               write(name,"('sigma-Omg',i1,'m.dat')")iabs(iom)
            endif
            open(10,file=name,status='unknown')
           write(10,*)'# Energy (eV) helicity X-sections(Angstrom^2)'
            do ie=1,nESigma
               E=EminSigma+dble(ie-1)*dE
               E_eV=E*ceVcm/conve
               xK2=2.d0*xmasa*E  ! all in a.u. 
               fac=pi/(xk2)
               fac=fac*convl*convl     ! converting to Angstrom^{-2}
               sigma=0.d0
               sigmaprod=0.d0
               do Jtot=max(0,iabs(iom)),Jtotmax
            
                  sigma=sigma+dble(2*Jtot+1)*PJ(ie,Jtot,iom)
                  sigmaprod=sigmaprod+dble(2*Jtot+1)*PJprod(ie,Jtot,iom)
               enddo
               write(10,'(300(1x,e15.7))')E_eV,sigma*fac,sigmaprod*fac
            enddo
            close(10)
         enddo
      endif
      
      stop
      end

!-------------------------------------------------------------------
      subroutine JshiftInterPol(PJ,PJcalc,Ecalc,BeJCalc,x,f
     &     ,EminSigma,EmaxSigma
     &     ,nEcalc,Jtotcalc,nESigma,Jtotmax,nCalc,netot,jref)
      implicit none
      integer :: nESigma,Jtotmax,nCalc,netot,jref,iom,iom0max,iom0,iom1
      integer :: iparsign
      integer :: nEcalc(nCalc,-jref:jref),Jtotcalc(nCalc)
      double precision :: PJ(nEsigma,0:Jtotmax,-jref:jref)
     &                   ,BeJcalc(ncalc)
      double precision :: PJcalc(netot,nCalc,-jref:jref)
     &                    ,Ecalc(netot,nCalc,-jref:jref)
      double precision :: x(netot),f(netot,2),EminSigma,EmaxSigma
      double precision :: dE,spl,Esigma,shift,Eshifted,fac0,fac1
      integer :: Jtot,nesta,J0,J1,i0,i1,ne,ie,iold,ieSigma,i
      double precision :: PJ0(nESigma), PJ1(nESigma)
      

      dE=(EmaxSigma-EminSigma)/dble(nESigma)

      do Jtot=0,Jtotmax
         iom0max=jref
         if(iom0max.gt.Jtot)iom0max=Jtot
         
         do iom=-iom0max,iom0max
            nesta=0
            J0=-10
            J1=-10
            i1=-1
            i0=-1
            do i=1,nCalc
               if(Jtotcalc(i).eq.Jtot)then
                  nesta=1
                  J0=Jtotcalc(i)
                  i0=i
                  iom0=iom
                  goto 1
               elseif(i.gt.1)then
                  if(Jtotcalc(i-1).lt.Jtot.and.Jtotcalc(i).gt.Jtot)then
                     J0=Jtotcalc(i-1)
                     J1=Jtotcalc(i)
                     i0=i-1
                     i1=i
                  endif
               endif
            enddo

 1          continue
            iom0=iom
            iom1=iom
!     chossing parity so the presence of Omega is or not as in present case
!            if(nesta.eq.0)then
!            if(iom.eq.0)then
!               iom0=iom
!               iom1=iom
!            elseif(iom.gt.0)then
!               iparsign=(-1)**Jtot
!               iom0=iom
!               if((-1)**J0.ne.iparsign)iom0=-iom
!               iom1=iom
!               if((-1)**J1.ne.iparsign)iom1=-iom
!            endif
!            endif
! 
            if(nesta.eq.1)then
               ne=nEcalc(i0,iom)
               do ie=1,ne
                  x(ie)=Ecalc(ie,i0,iom)
                  f(ie,1)=PJcalc(ie,i0,iom)
                  f(ie,2)=0.d0
               enddo
               call splset(f,x,ne,netot)
               iold=2
               do iEsigma=1,nESigma
                  Esigma=EminSigma+dble(iEsigma-1)*dE
                  call splinqq(f,x,iold,ne,Esigma,netot,spl)
                  PJ(iEsigma,Jtot,iom)=spl
               enddo
            elseif(nesta.eq.0.and.J0.ge.0.and.J1.ge.0)then
               PJ0(:)=0.d0
               PJ1(:)=0.d0
           ! PJ0
               ne=nEcalc(i0,iom0)
               do ie=1,ne
                  x(ie)=Ecalc(ie,i0,iom0)
                  f(ie,1)=PJcalc(ie,i0,iom0)
                  f(ie,2)=0.d0
               enddo
               call splset(f,x,ne,netot)
               iold=2
               do iEsigma=1,nESigma
                  Esigma=EminSigma+dble(iEsigma-1)*dE
                  shift=BeJCalc(i0)*dble(Jtot*(Jtot+1))
     &                 - BeJCalc(i0)*dble(J0*(J0+1))
                  Eshifted=Esigma-shift
                  spl=0.d0
                  if(Eshifted.ge.x(1).and.Eshifted.le.x(ne))then
                     call splinqq(f,x,iold,ne,Eshifted,netot,spl)
                     PJ0(iEsigma)=spl
                  elseif(Eshifted.lt.x(1))then
                     PJ0(iEsigma)=0.d0
                  elseif(Eshifted.gt.x(ne))then
                     PJ0(iEsigma)=f(ne,1)
                  endif
               enddo
           ! PJ1
               ne=nEcalc(i1,iom1)
               do ie=1,ne
                  x(ie)=Ecalc(ie,i1,iom1)
                  f(ie,1)=PJcalc(ie,i1,iom1)
                  f(ie,2)=0.d0
               enddo
               call splset(f,x,ne,netot)
               iold=2
               do iEsigma=1,nESigma
                  Esigma=EminSigma+dble(iEsigma-1)*dE
                  shift=BeJCalc(i1)*dble(Jtot*(Jtot+1))
     &              - BeJCalc(i1)*dble(J1*(J1+1))
                  Eshifted=Esigma-shift
                  spl=0.d0
                  if(Eshifted.ge.x(1).and.Eshifted.le.x(ne))then
                     call splinqq(f,x,iold,ne,Eshifted,netot,spl)
                     PJ1(iEsigma)=spl
                  elseif(Eshifted.lt.x(1))then
                     PJ1(iEsigma)=0.d0
                  elseif(Eshifted.gt.x(ne))then
                     PJ1(iEsigma)=f(ne,1)
                  endif
               enddo
            ! average
               do iEsigma=1,nESigma
                  fac0=dble(J1-Jtot)/dble(J1-J0)
                  fac1=dble(Jtot-J0)/dble(J1-J0)
                  PJ(iEsigma,Jtot,iom)=fac0*PJ0(iESigma)
     &                                +fac1*PJ1(iESigma)
               enddo
            elseif(Jtot.gt.Jtotcalc(nCalc))then ! Jtot > maximum J calculated extrapolation
!     PJ0
               i0=nCalc
               J0=Jtotcalc(i0)
               ne=nEcalc(i0,iom0)
               do ie=1,ne
                  x(ie)=Ecalc(ie,i0,iom0)
                  f(ie,1)=PJcalc(ie,i0,iom0)
                  f(ie,2)=0.d0
               enddo
               call splset(f,x,ne,netot)
               iold=2
               do iEsigma=1,nESigma
                  Esigma=EminSigma+dble(iEsigma-1)*dE
                  shift=BeJCalc(i0)*dble(Jtot*(Jtot+1))
     &              - BeJCalc(i0)*dble(J0*(J0+1))
                  Eshifted=Esigma-shift
                  spl=0.d0
                  if(Eshifted.ge.x(1).and.Eshifted.le.x(ne))then
                     call splinqq(f,x,iold,ne,Eshifted,netot,spl)
                     PJ(iEsigma,Jtot,iom)=spl
                  elseif(Eshifted.lt.x(1))then
                     PJ(iEsigma,Jtot,iom)=0.d0
                  elseif(Eshifted.gt.x(ne))then
                     PJ(iEsigma,Jtot,iom)=f(ne,1)
                  endif
               enddo
            
            endif
         enddo ! iom
      enddo                  ! Jtot
      return 
      end subroutine JshiftInterPol
