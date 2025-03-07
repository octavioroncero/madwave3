      program CRPprogram
      use ieee_arithmetic
      implicit real*8(a-h,o-z)

      character*50 name,viejo,posicion,system,frpini,fRgini,fcoefini
*
* calculates the Cumulative Reaction Probabilities (CRP) to calculate rates
*     input in eV of translational energy: has to be changed..
*     output energies are in eV traslational energy
*
*  second column in output CRPxxx files is: k^2 in Angstroms^2
*     sigma= CRP* pi/(k^2 (2j_ref+1) )
*

      parameter(S2max=1.d0) ! maximum value of Smat^2 allowed to avoid num. errors at low energies
      parameter(ndim=200)
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
      namelist/sigmaS2prod/Jtotmax,nCalc
     &     ,nESigma,EminSigma_eV,EmaxSigma_Ev
*********************************************************
* preparing dimensions depending on Jacobi coordinates used in propagation (iprod=1 -products- or iprod=2 reactants) 

      real*8,allocatable,dimension(:,:,:,:,:) :: CRP
      real*8,allocatable,dimension(:,:,:) :: S2J
      real*8,allocatable,dimension(:,:) :: ES2,SJ1,SJ2,PJ
      real*8,allocatable,dimension(:) :: S2mat,CRPv,CRPomg
      real*8,allocatable :: xx(:),f(:,:)
      integer,allocatable :: Jcalc(:)
      real*8, allocatable :: pm(:),BeJ(:)
      integer :: ifail,iom1real
      
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
      au2eV=27.2113957d0

!     reading data
         
         open(10,file='input.dat',status='old')
         read(10,nml = inputgridbase)
         close(10)
         iommaxprod=Jtot
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
         open(10,file='input.dat',status='old')
         read(10,nml = sigmaS2prod)
         close(10)

         write(6,*)'  --- input data ---'
         write(6,nml = inputgridbase)
         write(6,nml = inputcol)
         write(6,nml = inputflux)
         write(6,nml = inputprod)         
         write(6,nml = inputpotmass)
         write(6,nml = sigmaS2prod)
         write(6,*)'  --- end input data ---'

         emindistri=ekinmin_eV
         emaxdistri=ekinmax_eV
         nener=netot
         enermineV=EminSigma_eV        
         enermaxeV=EmaxSigma_eV
         nenerdif=nEsigma

**> preparing magnitudes for iprod=1 or 2

      xmtot=xm0+xm1+xm2
      if(iprod.eq.1)then              ! products Jacobi coordinates
         iom0=iommin
         iom1=iommax
         j00=jini
         j11=jmax
         nv0=nvini
         nv1=nvmax
         nel=nelecmax
         xmasa=(xm1*(xm0+xm2))/xmtot
         xmasaprod=(xm2*(xm0+xm1))/xmtot
         write(6,'("--> Using product Jacobi coordinates",/
     &         ,10x,"with BC mases= ",2(F12.8,1x))')xm0,xm2

      elseif(iprod.eq.2)then          ! reactants Jacobi coordinates
         iom0=iomminprod
         iom1=iommaxprod
         j00=jiniprod
         j11=jmaxprod
         nv0=nviniprod
         nv1=nvmaxprod
         nel=1
         xmasa=(xm2*(xm0+xm1))/xmtot
         xmasaprod=(xm1*(xm0+xm2))/xmtot
         write(6,'("--> Using reactant Jacobi coordinates",/
     &         ,10x,"with BC mases= ",2(F12.8,1x))')xm0,xm1
      endif
      write(6,*)
      write(6,*)'  iom0,iom1= ',iom0,iom1
      write(6,*)'  j00,j11= ',j00,j11
      write(6,*)'  nv0,nv1= ',nv0,nv1
      write(6,*)'  jref= ',jref
      write(6,*)'  nener,nenerdif= ',nener,nenerdif
      write(6,*)

**> allocating matrices

      allocate(S2J(nener,ncalc,j00:j11)
     &        ,ES2(nener,ncalc)
     &        ,SJ1(nenerdif,j00:j11)
     &        ,SJ2(nenerdif,j00:j11)
     &        ,PJ(nenerdif,0:Jtotmax)
     &        ,S2mat(j00:j11)
     &        ,CRPv(nv0:nv1)
     &        ,CRPomg(-jref:jref)
     &  ,CRP(nenerdif,j00:j11,nv0:nv1,nelecmax,-jref:jref)
     &     ,f(nener,2),xx(nener)
     &     ,Jcalc(ncalc),pm(0:ndim),BeJ(0:Jtotmax)
     &     ,stat=ierror)

      PJ(:,:)=0.d0
**>> J at which WP calculations have been performed
!**   >>  Reading calculated J
      
      open(10,file="CalculatedJ.dat",status="old",err=999)
      nJcalc=nCalc
      BeJ(:)=0.d0
      do i=1,nCalc
         read(10,*)Jcalc(i),BB
         BeJ(Jcalc(i))=BB
      enddo
      close(10)
      estep=(enermaxeV-enermineV)/dble(nenerdif-1)
      write(6,*)'  --> Jtot read, and J-shift. rot. cte (in eV)  '
      do j=1,nJcalc
         write(6,*)j,Jcalc(j),BeJ(Jcalc(j))
      enddo
      write(6,*)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   big Loop in iomref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(6,*)'enermineV,estep=',enermineV,estep

      CRP(:,:,:,:,:)=0.d0
      do iomref0=-jref,jref
      do ielec=1,nel
      do iv=nv0,nv1

**Reading reaction probabilities for calculated  J's
         S2J(:,:,:)=0.d0 !  S2J(ie,iJ,j,iom)=0.d0     
         do iJ=1,nJcalc
         Jtot0=Jcalc(iJ)
         if(Jtot0.ge.iabs(iomref0))then

               if(iomref0.eq.0)then
                  write(name,'("../Omg",i1,"/J",i3.3 
     &                 ,"/distriS2prod.v",i2.2
     &                 ,".e",i1.1)')
     &                           iabs(iomref0),Jtot0,iv,ielec
               elseif(iomref0.lt.0)then
                  if(mod(Jtot0,2).eq.0)then
                     write(name,'("../Omg",i1,"/m/J",i3.3
     &                 ,"/distriS2prod.v",i2.2
     &                 ,".e",i1.1)')
     &                    iabs(iomref0),Jtot0,iv,ielec
                  else
                     write(name,'("../Omg",i1,"/p/J",i3.3
     &                 ,"/distriS2prod.v",i2.2
     &                 ,".e",i1.1)')
     &                    iabs(iomref0),Jtot0,iv,ielec
                  endif
               elseif(iomref0.gt.0)then
                  if(mod(Jtot0,2).eq.0)then
                     write(name,'("../Omg",i1,"/p/J",i3.3
     &                 ,"/distriS2prod.v",i2.2
     &                 ,".e",i1.1)')
     &                    iabs(iomref0),Jtot0,iv,ielec
                  else
                     write(name,'("../Omg",i1,"/m/J",i3.3
     &                 ,"/distriS2prod.v",i2.2
     &                 ,".e",i1.1)')
     &                    iabs(iomref0),Jtot0,iv,ielec
                  endif
               endif

               write(6,*)iomref0,Jtot0,name
               open(5,file=name,status='old')

               if(iomref0.eq.0)then
                     write(name,'("../Omg",i1,"/J",i3.3
     &                 ,"/distriS2prod.elec")')iomref0,Jtot0
               elseif(iomref0.lt.0)then
                     write(name,'("../Omg",i1,"/m/J",i3.3
     &                 ,"/distriS2prod.elec")')
     &                    iabs(iomref0),Jtot0
               elseif(iomref0.gt.0)then
                     write(name,'("../Omg",i1,"/p/J",i3.3
     &                 ,"/distriS2prod.elec")')iomref0,Jtot0
               endif
!               open(4,file=name,status='old')


               do ie=1,nener
                  read(5,*)es2(ie,iJ),(S2J(ie,iJ,j),j=j00,j11)
                                    
!                  do j=j00,j11
!                     if(ieee_is_nan(S2J(ie,iJ,j)))then
!                        write(55,*)' NaN for ',ie,iJ,j
!                        S2J(ie,iJ,j)=0.d0
!                     endif
!                  enddo
                  
!                  read(4,*)xyz,S2tot
!!                  if(ieee_is_nan(S2tot))S2tot=0.d0
!!correcting too high reaction probabilities at low energies
!                  if(S2tot.gt.S2max.and.dabs(S2tot).gt.1.d-20)then
!                     fac=S2max/S2tot
!                     do j=j00,j11
!                        S2J(ie,iJ,j)=S2J(ie,iJ,j)*fac
!!                     if(ieee_is_nan(S2J(ie,iJ,j)))then
!!                        write(55,*)' NaN in renorm. for ',ie,iJ,j
!!                        S2J(ie,iJ,j)=0.d0
!!                     endif
!                     enddo
!                  endif
!!end-correction

               enddo  ! ie

               close(5)
               close(4)
               
         endif               ! Jtot0>iomref0   
         enddo  ! iJ index of calculated Jtot0

********Loop in energy  

**>> Interpolating the S^2-matrix for intermediate J's
**** state resolved opacity and integral cross section

         do j=j00,j11
            minJtot=iabs(iomref0)
           
            do Jtot0=minJtot,Jtotmax

                  S2mat(j)=0.d0
                  J1=0
                  J2=0
                  delta=0.d0
                  if(Jtot0.lt.Jcalc(nJcalc))then
                  do iJ=1,nJcalc-1
                  if(Jtot0.ge.Jcalc(iJ).and.Jtot0.lt.Jcalc(iJ+1))then
                     J1=Jcalc(iJ)
                     J2=Jcalc(iJ+1)
*Interpolation according to J1
                     do ie=1,nener
                        xx(ie)=es2(ie,iJ) 
                        f(ie,1)=S2J(ie,iJ,j)
                        f(ie,2)=0.d0
                     enddo
        
                     call splset(f,xx,nener,nener)
                  shift=BeJ(J1)*(dble(Jtot0*(Jtot0+1))-dble(J1*(J1+1)))
                     iold=2
                     do  iedif=1,nenerdif
                       enerdifeV=enermineV+dble(iedif-1)*estep
                       eshift=enerdifeV-shift
                       SJ1(iedif,j)=0.d0
                       if(eshift.gt.es2(nener,iJ))then
                         SJ1(iedif,j)=S2J(nener,iJ,j)
                       elseif(eshift.lt.es2(1,iJ))then
                         SJ1(iedif,j)=0.d0
                       else
                      call splinqq(f,xx,iold,nener,eshift,nener,spl)
                      SJ1(iedif,j)=spl
                       endif
                     enddo
*Interpolation according to J2
                     do ie=1,nener
                        xx(ie)=es2(ie,iJ+1)
                        f(ie,1)=S2J(ie,iJ+1,j)
                        f(ie,2)=0.d0
                     enddo
         
                     call splset(f,xx,nener,nener)
                  shift=BeJ(J2)*(dble(J2*(J2+1))-dble(Jtot0*(Jtot0+1)))
  
                     iold=2
                     do  iedif=1,nenerdif
                       enerdifeV=enermineV+dble(iedif-1)*estep
                       eshift=enerdifeV+shift 
                       SJ2(iedif,j)=0.d0
                       if(eshift.gt.es2(nener,iJ+1))then
                          SJ2(iedif,j)=S2J(nener,iJ+1,j)
                       elseif(eshift.lt.es2(1,iJ+1))then
                          SJ2(iedif,j)=0.d0
                       else
                      call splinqq(f,xx,iold,nener,eshift,nener,spl)
                      SJ2(iedif,j)=spl
                       endif
                     enddo
                     delta=dble(J2-J1)
*     building S2J
                  endif !
                  enddo ! do iJ calculated
                  elseif(Jtot0.ge.Jcalc(nJcalc))then
*Extrapolation according to J1
                     write(6,*)'  Jcalc < Jtot0= ',Jtot0
                     J1=Jcalc(nJcalc)
                     do ie=1,nener
                        xx(ie)=es2(ie,nJcalc) 
                        f(ie,1)=S2J(ie,nJcalc,j)
                        f(ie,2)=0.d0
                     enddo
        
                     call splset(f,xx,nener,nener)
                  shift=BeJ(J1)*(dble(Jtot0*(Jtot0+1))-dble(J1*(J1+1)))

                     iold=2
                     do  iedif=1,nenerdif
                       enerdifeV=enermineV+dble(iedif-1)*estep
                       eshift=enerdifeV-shift
                       SJ1(iedif,j)=0.d0
                       if(eshift.gt.es2(nener,nJcalc))then
                          SJ1(iedif,j)=S2J(nener,nJcalc,j)
                       elseif(eshift.lt.es2(1,nJcalc))then
                          SJ1(iedif,j)=0.d0
                       else
                       call splinqq(f,xx,iold,nener,eshift,nener,spl)
                         SJ1(iedif,j)=spl
                       endif
                     enddo
                  endif
!--> end intra-extrapolation

                  do iedif=1,nenerdif
                     enerdifeV=enermineV+dble(iedif-1)*estep
                     S2mat(j)=0.d0
                     if(Jtot0.eq.J1)then
                        S2mat(j)= SJ1(iedif,j)  
                     elseif(Jtot0.eq.J2)then
                        S2mat(j)= SJ2(iedif,j)
                     elseif(Jtot0.gt.Jcalc(nJcalc))then
                        S2mat(j)= SJ1(iedif,j)
                     else       !if(J1.ge.iomref0)then
                        
                        if(ieee_is_nan(SJ1(iedif,j)))then
                 write(55,*)' NaN in sJ1 for ',iedif,j,Jtot0,iomref0,iv
                           SJ1(iedif,j)=0.d0
                        endif
                        if(ieee_is_nan(SJ2(iedif,j)))then
                 write(55,*)' NaN in sJ2 for ',iedif,j,Jtot0,iomref0,iv
                           SJ2(iedif,j)=0.d0
                        endif
                     
                        S2mat(j) = dble(J2-Jtot0)*SJ1(iedif,j)
     &                           + dble(Jtot0-J1)*SJ2(iedif,j)
                        S2mat(j)=S2mat(j)/delta
                        

!                     else
!                        S2mat(j)=SJ2(iedif,j)       
                     endif          
                     if(S2mat(j).lt.0.d0)S2mat(j)=0.d0

                     CRP(iedif,j,iv,ielec,iomref0)=
     &                      CRP(iedif,j,iv,ielec,iomref0)
     &                    +S2mat(j)*dble(2*Jtot0+1)


                     PJ(iedif,Jtot0)=PJ(iedif,Jtot0)+S2mat(j)
                    
                  enddo ! iedif
    
            enddo  ! Jtot0
            
         enddo ! j

      enddo  ! iv
      enddo  ! ielec
      enddo  ! iomref0


* writting files

      open(7,file='CRPomg.res',status='new')
      open(9,file='CRPvib.res',status='new')
      open(8,file='CRPtot.res',status='new')
      open(10,file='PJ.res',status='new')
      ifile=10
      do iv=nv0,nv1
      do ielec=1,nel
         ifile=ifile+1
         write(name,'("CRP.vf",i2.2,".ef",i1.1)')iv,ielec
         open(ifile,file=name,status='unknown')
      enddo
      enddo

**>> Writting results   

      write(8,*)'# Energy(eV)  k^2(Angstrom^2)  CRP '
      do iedif=1,nenerdif
         enerdifeV=enermineV+dble(iedif-1)*estep 
         ener=enerdifeV/ceVcm*conve1
         xkini2=2.d0*xmasa*(ener)/hbr/hbr
         CRPtot=0.d0
         ifile=10
         do iomref0=-jref,jref
            CRPomg(iomref0)=0.d0
         enddo
         do Jt=0,Jtotmax
            if(PJ(iedif,Jt).lt.1.d-50)PJ(iedif,Jt)=0.d0
         enddo
      write(10,'(1000(1x,e16.7))')enerdifeV,(PJ(iedif,Jt),Jt=0,Jtotmax)
         do iv=nv0,nv1
            CRPv(iv)=0.d0

            do ielec=1,nel
               ifile=ifile+1
               do j=j00,j11
                  S2mat(j)=0.d0
                  do iomref0=-jref,jref
                
                     S2mat(j)=S2mat(j)+CRP(iedif,j,iv,ielec,iomref0)
                     CRPv(iv)=CRPv(iv)+CRP(iedif,j,iv,ielec,iomref0)
                     CRPomg(iomref0)=CRPomg(iomref0)
     &                       +CRP(iedif,j,iv,ielec,iomref0)
                     CRPtot=CRPtot+CRP(iedif,j,iv,ielec,iomref0)
                  
                  enddo
                  if(S2mat(j).lt.1.d-60)S2mat(j)=0.d0
               enddo  ! j
               write(ifile,'(50(1x,e15.7))')enerdifeV,xkini2
     &                                   ,(S2mat(j),j=j00,j11)
            enddo  ! ielec
         enddo  ! iv

         do iv=nv0,nv1
            if(CRPv(iv).lt.1.d-30)CRPv(iv)=0.d0
         enddo

         write(9,'(1000(1x,e15.7))')enerdifeV,xkini2
     &                    ,(CRPv(iv),iv=nv0,nv1)
         write(8,'(50(1x,e15.7))')enerdifeV,xkini2,CRPtot
         write(7,'(1000(1x,e15.7))')enerdifeV,xkini2
     &                    ,(CRPomg(iomref0),iomref0=-jref,jref)



      enddo ! ienerdif 
     
      stop
 999  write(6,*)' file "CalculatedJ.dat" not found '
      call flush(6)


      stop
 9999 format(80('*'),/,20x,' ',/
     &,30x,'for the ',a8,'  system',/,30x
     &     ,/,80('*'),//)

      end

