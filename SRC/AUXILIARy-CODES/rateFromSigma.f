      program rateFromSigma
      implicit none
!-------------------------------------------------------------------------------------
!     Calculates the reaction rate constant
!           from a reaction cross section
!     by integrating over a Boltzmann distribution of collision energy
!      
!      
!      K(T) =    int dE   sigma(E) E e^{-E/T}
!
!
!     Input  look at a,b
!     -----
!
!    a)  read in unit 5: energy(eV), cross section (Angstrom^2) 
! 
!    b)  no. of energies, E_min and E_max  read in namelist/sigmaFs2prod
!        no. of temperaturesi, Tmin, T_max read in namelist/InputRateSigma
!     kept in input.dat files used for madwave3 propagations.
!
!
!     Output files:
!     ------
!     rate.dat: 2 columns with temperature (in K) and rate constant in cm^3/s
!
!     Convert magnitudes internaly to atomic units
!
!-------------------------------------------------------------------------------------

      integer:: ncontfile,npun1min,ierror,i,iE,iJ,iJ0,iJ1,iEE
     &     ,iprod,nviniprod,nvmaxprod,jiniprod,jmaxprod
     &     ,iomminprod,iommaxprod,n2prod0,n2prod1,nangproj0,nangproj1
     &     ,Jtotmax,nCalc,nEsigma,nTemp,iTemp,extrapolation,nEextra
      double precision :: 
     &     xm1,xm0,xm2,EminSigma_eV,EmaxSigma_Ev
     &     ,VcutmaxeV,radcutmaxeV,rotcutmaxeV,pi,conve1,hbr,conve,convm
     &     ,ceVcm,au2eV,cm2K,EminSigma,EmaxSigma,xmtot,xmasa,xK2
     &     ,Rbalinprod,de
     &     ,fac,convl,E,E_eV,xx,TempMin,TempMax,rate,deltaTemp,deltaE
     &     ,Temp,Temp_au,au2cm,au2s,au2fs,au2K
     &     ,Afactor_angstroms2eV,Epower,ee_au,Afactor
     &     ,Eminextra_eV,deltaEextra,rate0,ee,sigma,xmasa_au,sigma_au
      character*50 :: system,name
      double precision, allocatable :: Esigma(:),sigmaE(:)
      
      namelist /inputpotmass/system,xm1,xm0,xm2
     &                      ,VcutmaxeV,radcutmaxeV,rotcutmaxeV
      namelist /inputprod/iprod
     &                   ,nviniprod,nvmaxprod
     &                   ,jiniprod,jmaxprod
     &     ,iomminprod,iommaxprod
     &     ,Rbalinprod,n2prod0,n2prod1,nangproj0,nangproj1
      namelist/sigmaS2prod/Jtotmax,nCalc
     &     ,nESigma,EminSigma_eV,EmaxSigma_Ev
      namelist/InputRateSigma/TempMin,TempMax,nTemp
     &     ,extrapolation,Afactor_angstroms2eV,Epower
     &     ,Eminextra_eV,nEextra
      write(6,'(40("-"))')
      write(6,'(10x,"Program rateFromSigma")')
      write(6,'(40("-"))')

!**>> constants

      pi = dacos(-1.d0)
      conve1 = 1.197d-4
      hbr = 1.d0 ! 0.063533625d0
      CONVL=.52917726D0
      CONVE=4.55633538D-6
      CONVM=.182288853D4
      ceVcm=1.23984245d-4
      au2eV=27.2113957d0
      cm2K=1.d0/0.629d0
      au2fs=24.188843265d-3
      au2K=315773.213d0
      au2cm=convl*1.d-8
      au2s=au2fs*1.d-15

      write(6,*) ' au_L^3/au_t = ',au2cm*au2cm*au2cm/au2s

!**>>     reading data

      open(10,file='input.dat',status='old')
      read(10,nml = inputpotmass)
      close(10)
      open(10,file='input.dat',status='old')
      read(10,nml = inputprod)
      close(10)
      open(10,file='input.dat',status='old')
      read(10,nml = sigmaS2prod,err=2)
      close(10)
      extrapolation=0
      open(10,file='input.dat',status='old')
      read(10,nml = InputRateSigma,err=1)
      close(10)

      write(6,*)'  --- input data ---'
      write(6,nml = inputprod)
      write(6,nml = inputpotmass)
      write(6,nml = sigmaS2prod)
      write(6,nml = InputRateSigma)
      write(6,*)'  --- end input data ---'
      deltaTemp=(TempMax-TempMin)/dble(nTemp-1)

      write(6,*)'calculation of rate constants by numerical integration'
      if(extrapolation.eq.0)then
         write(6,*)'    with NO extrapolation at low energies'
      else
         write(6,*)'   extrapolating with the function: A e**alpha '
         write(6,*)'      with A (Ang^2/ eV**alpha) = '
     &        ,Afactor_angstroms2eV,'  alpha= ',Epower
      endif
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
      xmasa_au=xmasa
!**   >> reading cross section

      allocate(Esigma(nESigma),sigmaE(nESigma)
     &     ,stat=ierror)

      read(5,*)
      do ie=1,nESigma
         read(5,*,err=3)Esigma(ie),sigmaE(ie)
         sigmaE(ie)=sigmaE(ie)/(convl*convl)  ! converting sigma to bohr^2
         Esigma(ie)=Esigma(ie)/au2eV
      enddo
      deltaE=(Esigma(nEsigma)-Esigma(1))/dble(nEsigma-1)
      if(extrapolation.ne.0)then
         Afactor=sigmaE(1)/(Esigma(1)**Epower)
      endif
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
!     rate constant: calculation and writting
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      open(10,file='rate.dat',status='unknown')
      write(10,*)'# Temperature (K)    rate constant (cm^3/s)'

      do itemp=1,ntemp
         Temp=TempMin+dble(itemp-1)*deltaTemp
         Temp_au=Temp/au2K
         fac=dsqrt(8.d0/(pi*xmasa*((temp_au)**3)))
! this is a.u. so that
         fac = fac *(au2cm/au2s)*(au2cm*au2cm)
! in units of  Angstroms ^3 /fs
!           and to pass to cm ^3 /s
!     1 Angstroms ^3 /fs = (10^{-8})**3 / 10^{-15} cm^3/s
!     = 10^{-9} cm^3/s
         rate=0.d0
         rate0=0.d0
         if(extrapolation.ne.0)then
            deltaEextra=(Esigma(1)*au2eV-eminextra_eV)/dble(nEextra-1)
            do ie=1,nEextra
               ee=Eminextra_eV+dble(ie-1)*deltaEextra
               ee_au=ee/au2eV
               xk2=2.d0*xmasa_au*ee_au  ! k^2 in atomic units
               sigma_au=Afactor*(ee_au**Epower) !
               xx=sigma_au*ee_au
               rate0=rate0+xx*dexp(-ee_au/Temp_au)
               if(itemp.eq.1)then
                   write(69,*)ee,Afactor*(ee_au**Epower)
               endif
            enddo
            rate0=rate0*deltaEextra/au2eV
         endif
         do ie=1,nESigma
            xx=sigmaE(ie)*Esigma(ie)
            rate=rate+xx*dexp(-Esigma(ie)/Temp_au)
            if(itemp.eq.1)then
               write(69,*)Esigma(ie)*au2ev,sigmaE(ie)
            endif
         enddo
         rate=rate*deltaE
         rate=rate+rate0
         write(10,'(300(1x,e15.7))')Temp,rate*fac
      enddo
      close(10)

      open(10,file='EdistriTmax',status='unknown')
      do ie=1,nESigma
         write(10,*)Esigma(ie)*au2eV,dexp(-Esigma(ie)/Temp_au),Temp
      enddo
      close(10)

      
      stop
 1    continue
      write(6,*)'  namelist "InputRateSigma" is absent in input.dat'
      write(6,*)' &InputRateSigma '
      write(6,*)' TempMin, TempMax,nTemp '
      write(6,*)' extrapolation,Afactor_angstroms2eV,Epower'
      write(6,*)'     ,Eminextra_eV,nEextra '
      call flush(6)
      stop
 2    continue
      write(6,*)'  namelist "sigmaS2prod" is absent in input.dat'
      write(6,*)' &sigmaS2prod'
      write(6,*)' Jtotmax,nCalc,nESigma,EminSigma_eV,EmaxSigma_Ev'
      call flush(6)
      stop
 3    continue
      write(6,*)' end of ifile=5 sigma input'
      write(6,*)' check number of energies'
      

      end

