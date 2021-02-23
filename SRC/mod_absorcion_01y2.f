      module mod_absorcion_01y2

!
!     Initialize absorption functions for r=R1 and R=R2
!

      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      
      implicit none
      real*8 :: absr1,absalp1
     &         ,absr2,absalp2     
      integer :: n1expo,ninterv1
     &          ,n2expo,ninterv2 
      real*8,allocatable :: absfr1(:),absfr2(:)
      integer :: ir1abs,ir2abs    

      save

      contains
********************************************
*   functions of   mod_absorcion_01y2.f    *
********************************************
!=======================================================

!--------------------------------------------------
      subroutine ini_absorcion
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2

      implicit none

      include "mpif.h"

      integer :: ir1,ir2     
      real*8 :: r1,potabs1,delta,deltan,expo,f1
      real*8 :: r2,potabs2,f2
      integer :: ierror
*********************************************************
      namelist /inputabsorcion/absr1,absalp1,n1expo,ninterv1
     &                        ,absr2,absalp2,n2expo,ninterv2     
      

      ninterv1=1
      ninterv2=1
      n1expo=2
      n2expo=2

      allocate(absfr1(npun1),absfr2(npun2)
     &       , stat=ierror)
      norealproc_mem=norealproc_mem
     &    +npun1+npun2
      write(6,*)'norealproc_mem= ',norealproc_mem
      write(6,*)'nointegerproc_mem= ',nointegerproc_mem
      call flush(6)
     
      write(6,'(40("="),/,10x,"absorption_mod",/,40("="))')
      write(6,*)
      write(6,*)'  absorption data'
      write(6,*)'  -------------------'
         open(10,file='input.dat',status='old')
         read(10,nml = inputabsorcion)
         write(6,nml = inputabsorcion)
         call flush(6)
         close(10)

* Initialization of absorption parameters
      ir1abs=npun1
      ir2abs=npun2
      do ir1=1,npun1
         f1=1.d0
         potabs1=0.d0
         r1=dble(ir1-1)*ah1+rmis1
         if(r1.ge.absr1)then
            if(ir1.lt.ir1abs)ir1abs=ir1
            delta=(r1-absr1)/dble(ninterv1)
            deltan=(delta)**n1expo
            expo = absalp1*deltan
            f1 = dexp(-expo)
c            potabs1=-hbr*expo/tstep
         endif
         absfr1(ir1)=f1
      enddo

      do ir2=1,npun2
         f2=1.d0
         potabs2=0.d0
         r2=dble(ir2-1)*ah2+rmis2
         if(r2.ge.absr2)then
            if(ir2.lt.ir2abs)ir2abs=ir2
            delta=(r2-absr2)/dble(ninterv2)
            deltan=(delta)**n2expo
            expo = absalp2*deltan
            f2 = dexp(-expo)
c            potabs2=-hbr*expo*0.1d0/tstep
         endif
         absfr2(ir2)=f2
      enddo

      write(6,*)'   --> Absorption starts at ir1,ir2 = '
     &                                  ,ir1abs,ir2abs


      return
      end subroutine ini_absorcion
!--------------------------------------------------
!=======================================================
      end module mod_absorcion_01y2
