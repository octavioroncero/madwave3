      program CIPtot_From_OmgRef
!!      use ieee_arithmetic
      implicit real*8(a-h,o-z)

      character*50 name,viejo,posicion,system,frpini,fRgini,fcoefini
*
*     Sum the  Cumulative Inelastic Probabilities (CIP)
*     for a given inoitial Omgref, summing in Jtot and parity
*     obtained in the program CIP-OmgRef.f
*     input in eV of translational energy: has to be changed..
*     output energies are in eV traslational energy
*
*     second column in output CIPxxx files is: k^2 in Angstroms^2
*        CIP is normalized, i.e. it is divided by the (2j_ref+1)
*     sigma= CIP* pi/(k^2 (2j_ref+1) )
*
*

      parameter(S2max=1.d0) ! maximum value of Smat^2 allowed to avoid num. errors at low energies
      parameter(ndim=200)
*********************************************************
      namelist /inputgridbase/npun1,rmis1,rfin1,npun1min
     &                       ,npun2,rmis2,rfin2
     &                       ,nangu,nangplot
     &     ,Jtot,iparity,inc,nelecmax,iommin,iommax,j0
     &     ,jini,jmax,nvini,nvmax
     &                       ,nvref,jref,iomref,ielecref
      namelist/sigmaS2prod/Jtotmax,nCalc
     &     ,nESigma,EminSigma_eV,EmaxSigma_Ev
*********************************************************
* preparing dimensions depending on Jacobi coordinates used in propagation (iprod=1 -products- or iprod=2 reactants) 

      real*8,allocatable :: CIP(:,:,:), CIPtot(:,:)
      real*8,allocatable :: xx(:),f(:,:),Energy(:),xK2(:)
      integer, allocatable :: iomRefCalc(:)
      
!     reading data
         
         open(10,file='input.dat',status='old')
         read(10,nml = inputgridbase)
         close(10)
         open(10,file='input.dat',status='old')
         read(10,nml = sigmaS2prod)
         close(10)
         write(6,*)'  --- input data ---'
         write(6,nml = inputgridbase)
         write(6,nml = sigmaS2prod)
         write(6,*)'  --- end input data ---'

**> preparing magnitudes for iprod=1 or 2

      if(iprod.eq.1)then              ! products Jacobi coordinates
         write(6,*)'  For inelastic calculations iprod must be = 2'
         write(6,*)'      to use reactant Jacobi coordinates '
         call flush(6)
         stop

      else  !!! if(iprod.eq.2)then          ! reactants Jacobi coordinates
         iom0=iommin
         iom1=iommax
         j00=jini
         j11=jmax
         nv0=nvini
         nv1=nvmax
         write(6,'("--> Using reactant Jacobi coordinates",/)')
      endif
      write(6,*)
      write(6,*)'  iom0,iom1= ',iom0,iom1
      write(6,*)'  j00,j11= ',j00,j11
      write(6,*)'  nv0,nv1= ',nv0,nv1
      write(6,*)'  jref= ',jref
      write(6,*)

 !     Reading the Omgref for which  CIP have been Calculated

       open(10,file='CalculatedOmgRef.dat',status='old')
       read(10,*)nOmgref
       write(6,*)'  Reading CIP for ',nOmgref,'  Omega projections'
       allocate(iomRefCalc(nOmgref))
       do i=1,nOmgref
          read(10,*)iomRefCalc(i)
          write(6,*)'  Omg= ',iomRefCalc(i)
       enddo
       close(10)
       if(iomRefCalc(1).ne.0)then
          write(6,*)' First Calculated Omegref must be zero '
          call flush(6)
          stop
       endif
        
!     allocating matrices
       allocate(
     &     CIP(nESigma,j00:j11,nOmgref)
     &     ,CIPtot(nESigma,j00:j11)
     &     ,Energy(neSigma),xK2(neSigma)
     &     ,f(nEsigma,2),xx(nEsigma)
     &     ,stat=ierror)

       CIP(:,:,:)=0.d0
       CIPtot(:,:)=0.d0
       Energy(:)=0.d0
       xK2(:)=0.d0       

       do ielec=1,nelecmax
       do iv=nv0,nv1
!     Reading properties for OmgRef=0

          i=1
          iomref=0
          factor=0.d0
          write(name,'("../CIP-Omg",i2.2,"/CIP.Omgref",i2.2
     &       ,".vf",i2.2,".ef",i1.1)')iomref,iomref,iv,ielec
          write(6,*)' Opening file ',name
          open(10,file=name,status='old')
          do iE=1,neSigma
             read(10,*)energy(ie),xK2(ie),(CIP(ie,j,1),j=j00,j11)
           enddo                    ! ie
          close(10)

          do i=2,nOmgref
             iomref=iomRefCalc(i)
             factor=factor+2.d0
             write(name,'("../CIP-Omg",i2.2,"/CIP.Omgref",i2.2
     &           ,".vf",i2.2,".ef",i1.1)')iomref,iomref,iv,ielec
             write(6,*)' Opening file ',name
             open(10,file=name,status='old')
             do iE=1,neSigma
                read(10,*)ee,xk,(CIP(ie,j,i),j=j00,j11)
                if(dabs(ee-energy(ie)).gt.1.d-5)then
                   write(6,*)' Energy grid is not the same '
                   write(6,*) 'ee,energy= ',ee,energy(ie),ie,j00,j11
                   call flush(6)
                   stop
                endif
             enddo
! ie
             close(10)

          enddo                 ! i reading iomref

!     Interpolating for all Omega_ref values

          ! contribuions for Omega_ref
          if(jref.gt.0)then
          do iom=0,jref
             ifail=1
             i1=-1
             i2=-1
             
             do  i=1,nOmgref
                if(iom.eq.iomRefCalc(i))then
                   ifail=0
                   i1=i
                   i2=i
                endif
             enddo
             if(ifail.eq.1)then
                do  i=2,nOmgref
                if(iom.gt.iomRefCalc(i-1).and.iom.lt.iomRefCalc(i))then
                   ifail=0
                   i1=i-1
                   i2=i
                endif
                enddo
             endif

             if(i1.eq.i2.and.ifail.eq.0)then
                do iE=1,neSigma
                   do j=j00,j11
                      CIPtot(ie,j)=CIPtot(ie,j)+CIP(ie,j,i1)
                   enddo
                enddo
             elseif(i1.ne.i2.and.ifail.eq.0)then
                delta=iomRefCalc(i2)-iomRefCalc(i1)
                del2=iom-iomRefCalc(i1)
                del1=iomRefCalc(i2)-iom
                do iE=1,neSigma
                   do j=j00,j11
                      xxx=CIP(ie,j,i1)*del1+CIP(ie,j,i2)*del2
                      CIPtot(ie,j)=CIPtot(ie,j)+xxx/delta
                   enddo
                enddo
             else
                write(6,*)' Problem finding interval for iom=',iom
                call flush(6)
                stop
             endif
             write(6,*)' iom ',iom,' i1,i2= ',i1,i2
          enddo !iom=1,jrec
          endif                    ! jref>0
       
          write(name,'("CIP.vf",i2.2,".ef",i1.1)')iv,ielec
          open(11,file=name,status='new')

          factor=dble(2*jref+1)
          write(6,*)' normalization factor= ',factor
          do iE=1,neSigma
             write(11,'(1000(1x,e15.7))')energy(ie),xK2(ie)
     &           ,(CIPtot(ie,j)/factor,j=j00,j11)
          enddo                    ! ie
          close(11)
          
!  ending
       enddo !iv
       enddo !ielec
       

      stop
      end
