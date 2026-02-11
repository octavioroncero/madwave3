      program CrossPhoto
!     *********************************************************************
!     Initialize the potential, reactant and product wf, and dipole transition moment
!     to be read to proceed for several wave packet calculations witj MadWave3
!     for AB+C --> AC + C
!     for 01+2 --> 02 + 1
!     for different initial states
!     *********************************************************************

      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      implicit none

      real*8 :: r1flux_ang,ekinmin_eV,ekinmax_eV
      real*8 :: ekinmin,ekinmax,energy_ini
      integer :: netot,ncontfile
      real*8,allocatable :: etot(:),flux_total(:),photon_energy_au(:)
      real*8,allocatable :: flux_arr2y3(:,:),flux_diat_arr1(:,:)
      real*8,allocatable :: flux_diat_arr2(:,:),flux_diat_arr3(:,:)
      real*8,allocatable :: flux_AyByC(:,:)
      real*8 :: ephoton,lambda
      integer ::ielec,ie
 
************************************************************
      namelist /inputflux/r1flux_ang,netot,ekinmin_eV,ekinmax_eV
     &                   ,ncontfile
************************************************************
! initialization
      call input_grid
      call pot0
!      call  paralelizacion
      write(6,*)
      write(6,*)'  --- Initialization of the potential ---'
      write(6,*)

      nelec=nelecmax
      call setxbcpotele(iomdiat,iomatom,sigdiat,sigatom
     &                     ,nelec,nelecmax)
!     rovibrational states of reactants and products
         
       if( xm2.gt.1.d-3)then
          call radial_functions01_read

          call product_radialf_read
       endif
!     initialize total flux quantities

      open(10,file='input.dat',status='old')
         read(10,nml = inputflux)
         write(6,nml = inputflux)
         call flush(6)
      close(10)
      ekinmin=ekinmin_eV*ev2cm*conve1
      ekinmax=ekinmax_eV*ev2cm*conve1
      
      allocate(  etot(netot),flux_total(netot)
     &        ,photon_energy_au(netot)
     &        ,flux_arr2y3(netot,nelec)
     &        ,flux_diat_arr1(netot,nelec)
     &        ,flux_diat_arr2(netot,nelec)
     &        ,flux_diat_arr3(netot,nelec) 
     &        ,flux_AyByC(netot,nelec) )
!       
!     reading total fluxes, tranforming energies to a.u., flux(in zots) ---> cross section in cm^2
!       
      call read_S2prod(etot,flux_total,photon_energy_au,energy_ini
     &                ,netot,xm2)

! from here only for triatomic systems
      if (xm2 .gt. 1.d-3) then
!
!     Reading fluxes in rovibrational states of (reactants, iarr=1) AB + C

         call read_diat_reac_flux(flux_diat_arr1,energy_ini,nelec,netot)

!
!     Reading total fluxes on products rearrangements (iar=2 and 3)

         call read_products_flux(flux_arr2y3,energy_ini,nelec,netot)

!
!     Reading fluxes on rovibrational states of products in iarr=2

        call read_diat_prod2_flux(flux_diat_arr2,energy_ini,nelec,netot)
      
!     Infering flux of rearrangement 3 (iarr3)

         call diat_prod3_flux(flux_arr2y3,flux_diat_arr1
     &                    ,flux_diat_arr2,flux_diat_arr3
     &                     ,nelec,netot)
      
!     Infering flux for total fragmentantion

         call total_AyByC_flux(flux_arr2y3,flux_diat_arr2,flux_diat_arr3
     &     ,flux_AyByC,nelec,netot)

      end if ! xm2>1d-3
!     printing total fluxes

      write(name,"('xsection-total')")
      open(10,file=name,status='unknown')
         write(10,*)'# photon energy (a.u.) wavelength (Angstrom)'
     &            ,'  total xsection(Angstrom^2)'
         do ie=1,netot
            lambda=1.d8/(photon_energy_au(ie)*hartree2cm)
            write(10,'(10(1x,e20.7e3))')photon_energy_au(ie)
     &                   ,lambda          , flux_total(ie)        

         enddo 
      close(10)

!     printing electronic fluxes

      if (xm2 .gt. 1.d-3) then
         do ielec=1,nelec
            write(name,"('xsection-elec',i2.2)")ielec
            open(10,file=name,status='unknown')
            write(10,*)'# photon energy (a.u.) wavelength (Angstom)'
     &        ,'arr_i(i=1-3)-xsection  A+B+C xsection (Angstrom^2)'
               do ie=1,netot
                  lambda=1.d8/(photon_energy_au(ie)*hartree2cm)
                  write(10,'(10(1x,e20.7e3))')photon_energy_au(ie)
     &                             ,lambda
     &                             , flux_diat_arr1(ie,ielec)        
     &                             , flux_diat_arr2(ie,ielec)        
     &                             , flux_diat_arr3(ie,ielec)        
     &                             , flux_AyByC(ie,ielec)        
     &                             , flux_arr2y3(ie,ielec)
               enddo 
            close(10)
         end do
      end if

         
       stop
       end program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_S2prod(etot,flux_total,photon_energy_au
     &                       ,energy_ini,netot,xm2)
      use mod_gridYpara_01y2
      implicit none
      integer :: netot
      real*8 :: etot(netot),flux_total(netot)
      real*8 :: photon_energy_au(netot),xm2
      real*8 :: energy_ini_eV,energy_ini
      character*50 :: name
      integer :: ie
      real*8 :: ephoton,factor,xxx

      etot(:)=0.d0
      flux_total(:)=0.d0
      
      open(10,file='auto',status='old')
      read(10,*)xxx,energy_ini_eV
      close(10)

      energy_ini=energy_ini_eV/au2eV

      factor=(4.d0/(pi*hbr*zot2au))*CSconstant_au
     &     *(convl*convl)
      
      open(10,file='S2prod',status='old')
      do ie=1,netot
         if(xm2 .gt. 1.d-3)then
            read(10,*)etot(ie),xxx,flux_total(ie)
         else
            read(10,*)etot(ie),flux_total(ie)
         endif  
         ephoton=etot(ie)/au2eV-energy_ini
         flux_total(ie)=ephoton*factor*flux_total(ie)
         photon_energy_au(ie)=ephoton
      end do
      close(10)
     

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_diat_reac_flux(flux_diat_arr1,energy_ini
     &                                   ,nelec,netot)
      use mod_gridYpara_01y2
      implicit none
      integer :: nelec,netot
      real*8 :: flux_diat_arr1(netot,nelec)
      character*50 :: name
      integer :: ie,ielec
      real*8 :: ephoton,factor,xxx,eee,fff,energy_ini

      flux_diat_arr1(:,:)=0.d0
      
      factor=(4.d0/(pi*hbr*zot2au))*CSconstant_au
     &     *(convl*convl)

       do ielec=1,nelec
          write(name,'("distriS2reac.elec",i2.2)')ielec
          open (10,file=name,status='old')
          do ie=1,netot
             read(10,*)eee,fff
             ephoton=eee/au2eV-energy_ini
             flux_diat_arr1(ie,ielec)=ephoton*factor*fff
          end do
          close(10)
       end do
     
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine  read_products_flux(flux_arr2y3,energy_ini,nelec,netot)
      use mod_gridYpara_01y2
      implicit none
      integer :: nelec,netot
      real*8 :: flux_arr2y3(netot,nelec)
      real*8 :: energy_ini
      character*50 :: name
      integer :: ie,ielec
      real*8 :: ephoton,factor,xxx,eee,fff
      
      flux_arr2y3(:,:)=0.d0
      factor=(4.d0/(pi*hbr*zot2au))*CSconstant_au
     &     *(convl*convl)

       write(name,'("S2prodelec")')
       open (10,file=name,status='old')
       do ie=1,netot
          read(10,*)eee,(flux_arr2y3(ie,ielec),ielec=1,nelec)
          ephoton=eee/au2eV-energy_ini
          do ielec=1,nelec
             flux_arr2y3(ie,ielec)=ephoton*factor*flux_arr2y3(ie,ielec)
          end do
       end do
       close(10)
     
      return
      end
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_diat_prod2_flux(flux_diat_arr2,energy_ini
     &                              ,nelec,netot)
      use mod_gridYpara_01y2
      implicit none
      integer :: nelec,netot
      real*8 :: flux_diat_arr2(netot,nelec)
      real*8 :: energy_ini
      character*50 :: name
      integer :: ie,ielec,iv,ivtot
      real*8 :: ephoton,factor,x2,x3,x4,eee,fff,etot
      real*8 :: fluxvib(nviniprod:nvmaxprod)

      flux_diat_arr2(:,:)=0.d0
      
      factor=(4.d0/(pi*hbr*zot2au))*CSconstant_au
     &     *(convl*convl)

      write(6,*)' nviniprod,nvmaxprod=', nviniprod,nvmaxprod
      write(6,*)' nelec,maxvibleves=',nelec,max_viblevels
      open(10,file='S2prod',status='old')
      do ie=1,netot
         read(10,*)etot,x2,x3,x4
     &             ,(fluxvib(ivtot),ivtot=nviniprod,nvmaxprod)
         ephoton=etot/au2eV-energy_ini
         do ivtot=nviniprod,nvmaxprod
            fluxvib(ivtot)=ephoton*factor*fluxvib(ivtot)
         end do

         ivtot=nviniprod-1
         do ielec=1,nelec
            if(max_viblevels(ielec).gt.0)then
               do iv=1,max_viblevels(ielec)
                  ivtot=ivtot+1
                  if(ie.eq.1)write(6,*)ivtot,iv,ielec
                  flux_diat_arr2(ie,ielec)=flux_diat_arr2(ie,ielec)
     &                 +fluxvib(ivtot)
               end do
            end if
         end do

      end do
      close(10)
      
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine  diat_prod3_flux(flux_arr2y3
     &                     ,flux_diat_arr1
     &                     ,flux_diat_arr2
     &                     ,flux_diat_arr3
     &                     ,nelec,netot)
      use mod_gridYpara_01y2
      use mod_pot_01y2, only: xm0,xm1,xm2
      implicit none
      integer :: nelec,netot
      real*8 :: flux_diat_arr1(netot,nelec)
      real*8 :: flux_diat_arr2(netot,nelec)
      real*8 :: flux_diat_arr3(netot,nelec)
      real*8 :: flux_arr2y3(netot,nelec)
      real*8 :: energy_ini
      character*50 :: name
      integer :: ie,ielec,iv,ivtot
!
      flux_diat_arr3(:,:)=0.d0

      if(inc.eq.2.or.xm0.eq.xm2.or.xm1.eq.xm2)then

         do ie=1,netot
            do ielec=1,nelec
               flux_diat_arr3(ie,ielec) = flux_diat_arr1(ie,ielec)
            end do
         end do


      else

         do ie=1,netot
            do ielec=1,nelec
               flux_diat_arr3(ie,ielec) = flux_arr2y3(ie,ielec)
     &                                  - flux_diat_arr2(ie,ielec)
            end do
         end do

      end if

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine   total_AyByC_flux(flux_arr2y3,flux_diat_arr2
     &                     ,flux_diat_arr3,flux_AyByC
     &                     ,nelec,netot)
      use mod_gridYpara_01y2
      implicit none
      integer :: nelec,netot
      real*8 :: flux_diat_arr2(netot,nelec)
      real*8 :: flux_diat_arr3(netot,nelec)
      real*8 :: flux_arr2y3(netot,nelec)
      real*8 :: flux_AyByC(netot,nelec)
      real*8 :: energy_ini
      character*50 :: name
      integer :: ie,ielec,iv,ivtot
!
      flux_AyByC(:,:)=0.d0

      do ie=1,netot
         do ielec=1,nelec
            flux_AyByC(ie,ielec) = flux_arr2y3(ie,ielec)
     &                           - flux_diat_arr3(ie,ielec)
     &                           - flux_diat_arr2(ie,ielec)
         end do
      end do

      return
      end
