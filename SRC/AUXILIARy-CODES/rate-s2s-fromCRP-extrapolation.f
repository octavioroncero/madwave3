      program rates_from_CRP_Extrapolation

      implicit none
!-------------------------------------------------------------------------------------
!     Calculates the reaction state-2-state rate constant
!           from a Cumulative Reaction probability
!     by integrating over a Boltzmann distribution of collision energy
!      
!      
!      K(T) =    int dE   sigma(E) E e^{-E/T}
!
!
!     Input  look at a,b
!     -----
!
!     a)  read in unit 5 all individual state-2-state CRP files
!          from initial state (nvref,jref,ieleref) of reactants
!         to 
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

      real*8, allocatable :: CRP(:,:,:,:),rate(:),energy(:)
      real*8, allocatable :: xk2(:),ratetot(:),ratetotlow(:)
     &                      ,Afactor(:,:,:)
      character*40 name,system
      real*8 :: rmis1,rfin1,rmis2,rfin2,xm1,xm0,xm2
     &     ,VcutmaxeV,radcutmaxeV,rotcutmaxeV,Rbalinprod
      real*8 :: EminSigma_eV,EmaxSigma_Ev
      real*8 :: TempMin,TempMax,Afactor_angstroms2eV,Eminextra_eV
      integer :: npun1,npun1min,npun2,nangu,nangplot
     &     ,Jtot,iparity,inc,nelecmax,iommin,iommax,j0
     &     ,jini,jmax,nvini,nvmax
     &     ,nvref,jref,iomref,ielecref
      integer :: iprod,nviniprod,nvmaxprod,jiniprod,jmaxprod
     &     ,iomminprod,iommaxprod,n2prod0,n2prod1,nangproj0,nangproj1
      integer :: Jtotmax,nCalc,nESigma,nTemp
     &     ,extrapolation,nEextra
      integer :: nvinirate,nvmaxrate,jinirate,jmaxrate,iomminrate
     &          ,iommaxrate
      real*8 :: deltaE,deltaTemp,E,temp,ratevib,yyy,deltaEextra,xxk2
      real*8 :: xx,eexpo,fac,temp_au,Epower,e_au,sigma_ang2,sigma_au
      integer :: iTemp,ie,j,iv,ielec,ifile
      real*8 :: pi,conve1,hbr, CONVL,CONVE,CONVM,ceVcm,au2eV
     &         ,cm2K,au2fs,au2K,au2cm,au2s,xmasa,xmtot,xmasa_au
      
      namelist /inputgridbase/npun1,rmis1,rfin1,npun1min
     &                       ,npun2,rmis2,rfin2
     &                       ,nangu,nangplot
     &     ,Jtot,iparity,inc,nelecmax,iommin,iommax,j0
     &     ,jini,jmax,nvini,nvmax
     &                       ,nvref,jref,iomref,ielecref
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
      write(6,'(10x,"Program rate-s2s-FromCRP")')
      write(6,'(40("-"))')
      call flush(6)
**>> constants

      pi = dacos(-1.d0)
      conve1 = 1.197d-4
      hbr = 1.d0    !     0.063533625d0
      CONVL=.52917726D0
      CONVE=4.55633538D-6
      CONVM=.182288853D4
      ceVcm=1.23984245d-4
      au2eV=27.2113957d0
      
      ceVcm=1.23984245d-4
      au2eV=27.2113957d0
      cm2K=1.d0/0.629d0
      au2fs=24.188843265d-3
      au2K=315773.213d0
      au2cm=convl*1.d-8
      au2s=au2fs*1.d-15

!**>>     reading data

      open(10,file='input.dat',status='old')
      read(10,nml = inputgridbase)
      close(10)
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
      
!**>> reduced mass
      xm0=xm0*convm
      xm1=xm1*convm
      xm2=xm2*convm
      xmtot=xm0+xm1+xm2
      if(iprod.eq.1)then
        xmasa=(xm1*(xm0+xm2))/xmtot
        write(6,'("--> Using product Jacobi coordinates",/
     &,10x,"with BC mases= ",2(F12.8,1x))')xm0/convm,xm2/convm
        nvmaxrate=nvmax
        nvinirate=nvini
        jmaxrate=jmax
        jinirate=jini
      else               
        xmasa=(xm2*(xm0+xm1))/xmtot
        write(6,'("--> Using reactant Jacobi coordinates",/
     &         ,10x,"with BC mases= ",2(F12.8,1x))')xm0/convm,xm1/convm 
        nvmaxrate=nvmaxprod
        nvinirate=nviniprod
        jmaxrate=jmaxprod
        jinirate=jiniprod
      endif
      xmasa_au=xmasa
!     allocating matrices

      allocate(
     &  CRP(nEsigma,jinirate:jmaxrate,nvinirate:nvmaxrate,nelecmax)
     &     ,rate(jinirate:jmaxrate)
     &    ,Afactor(jinirate:jmaxrate,nvinirate:nvmaxrate,nelecmax)
     & , energy(nEsigma),xk2(nEsigma)
     & ,ratetot(nTemp),ratetotlow(nTemp))

!     reading files (energy in eV, cross section in Angstroms)
!        converting to au

      ifile=9
      do ielec=1,nelecmax
      do iv=nvinirate,nvmaxrate
        
         write(name,'("CRP.vf",i2.2,".ef",i1.1)')iv,ielec
         write(6,*)ielec,iv,nEsigma,' reading in ',name
         call flush(6)
         open(ifile,file=name,status='unknown')
         do ie=1,nEsigma
            read(ifile,*)energy(ie),xk2(ie)
     &                  ,(CRP(ie,j,iv,ielec),j=jinirate,jmaxrate)
!            write(6,*)energy(ie),xk2(ie)
!     &                  ,(CRP(ie,j,iv,ielec),j=jinirate,jmaxrate)
         enddo
         close(ifile)

      enddo
      enddo

!     converting energy to a.u.
      
      do ie=1,nEsigma
         energy(ie)=energy(ie)/au2eV
      enddo
      deltaE=energy(2)-energy(1)
   
!     preparing Afactor for extrapolation as sigma= A e**Epower

      Afactor(:,:,:)=0.d0
      deltaEextra=0.d0
      if(extrapolation.ne.0)then
         deltaEextra=(Energy(1)-eminextra_eV/au2eV)/dble(nEextra)
         xxk2=2.d0*xmasa_au*energy(1)
         do ielec=1,nelecmax
            do iv=nvinirate,nvmaxrate
               do j=jinirate,jmaxrate
                  if(CRP(1,j,iv,ielec).gt.1.d-10)then
                     fac=(energy(1)**(Epower))*energy(1)
                     Afactor(j,iv,ielec)=CRP(1,j,iv,ielec)/(fac)
                     write(6,*)j,iv,ielec,Afactor(j,iv,ielec)
                  endif
               enddo
            enddo
         enddo
      endif
      
* integration in energy

      ratetot(:)=0.d0
      ratetotlow(:)=0.d0
      
      do ielec=1,nelecmax
      do iv=nvinirate,nvmaxrate
         write(name,"('ratevib.e'i2.2,'.v',i2.2)")ielec,iv
         open(9,file=name,status='unknown')
         write(name,"('rate.e'i2.2,'.v',i2.2)")ielec,iv
         open(10,file=name,status='unknown')
     
         do iTemp=1,ntemp
            temp=TempMin+dble(itemp-1)*deltaTemp
            temp_au=temp/au2K
!            fac=dsqrt(8.d0/(pi*xmasa*((temp_au)**3)))
            fac=dsqrt(2.d0*pi/(((xmasa_au*temp_au)**3)))
! this is a.u. so that
            fac = fac *(au2cm/au2s)*(au2cm*au2cm)
! in units of  Angstroms ^3 /fs
!           and to pass to cm ^3 /s
!     1 Angstroms ^3 /fs = (10^{-8})**3 / 10^{-15} cm^3/s
!     = 10^{-9} cm^3/s
            
            ratevib=0.d0
            do j=jinirate,jmaxrate
               rate(j)=0.d0
                
               if(extrapolation.ne.0)then
                  yyy=0.d0
                  do ie=1,nEextra
                     e_au=Eminextra_eV/au2eV+dble(ie-1)*deltaEextra                     
                     xxk2=2.d0*xmasa_au*e_au
                     sigma_au=(e_au**Epower)                     
                     xx= sigma_au*e_au
                     yyy=yyy+xx*dexp(-e_au/Temp_au)
                  enddo
                  rate(j)=yyy*Afactor(j,iv,ielec)*deltaEextra
                  ratetotlow(itemp)=ratetotlow(itemp)
     &               +yyy*Afactor(j,iv,ielec)*deltaEextra
               endif
               
               yyy=0.d0
               do ie=1,nEsigma
                  eexpo=dexp(-energy(ie)/Temp_au)                 
                  yyy=yyy+eexpo*CRP(ie,j,iv,ielec)
               enddo
               rate(j)=rate(j)+yyy*deltaE
               
               ratevib=ratevib+rate(j)
            enddo
            ratetot(itemp)=ratetot(iTemp)+ratevib
            
            write(10,'(1000(1x,e15.7))')temp
     &            ,(rate(j)*fac,j=jinirate,jmaxrate)

            write(9,*)temp,ratevib*fac
         enddo  ! itemp

      enddo !j
      enddo !iv

      open(16,file='ratetot.res',status='unknown')
      do itemp=1,nTemp
         temp=TempMin+dble(itemp-1)*deltaTemp
         temp_au=temp/au2K
!         fac=dsqrt(8.d0/(pi*xmasa*((temp_au)**3)))
          fac=dsqrt(2.d0*pi/(((xmasa_au*temp_au)**3)))
         fac = fac *(au2cm/au2s)*(au2cm*au2cm)
         write(16,*)temp,ratetot(itemp)*fac,ratetotlow(itemp)*fac
      enddo
      close(16)
      
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
      stop
      end
