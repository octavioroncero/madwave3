       module mod_colini_01y2
!
!     generation of photo initiated process
!     iphoto=1 ----> rpaq0=   <Jtot | d . d | Psi^Jtotini_nvbound >
!     iphoto=2 ---->  rpaq0=  Psi^Jtotini_nvbound
!
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      
      implicit none

      save
      real*8 :: Rcolini,Rcolini_ang,ecol,ecol_ev,deltaE,deltaE_ev
      real*8 :: xk0,xmasa0,alpha0,factor0
      integer :: il0

      contains
********************************************
*   functions of   mod_colini_01y2.f     *
********************************************
!=======================================================
!--------------------------------------------------
      subroutine set_colini
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      implicit none
      include "mpif.h"
       
      double precision :: det,xxx,r2,yyy,xr2
      integer :: i,ierr,ir2
*********************************************************
      namelist /inputcol/Rcolini_ang,ecol_ev,deltaE_ev
      
!  reading data in namelist
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
      
      xmasa0=xm2reac
      if(iprod.eq.1)then
          xmasa0=xm2prod
      endif
      xk0=dsqrt(2.d0*xmasa0*ecol)/hbr
      alpha0=dsqrt(2.d0*ecol/xmasa0)*hbr/deltae
      factor0=(2.d0/pi/alpha0/alpha0)**(0.25d0)
      DET=DBLE(Jtot*(Jtot+1)+jref*(jref+1)-2*iomref*iomref)
      il0=0.5D0*(-1.D0+DSQRT(1.D0+4.d0*DET))+0.5D0
      
      write(6,*)
      write(6,*)'      ** Initial wavepacket for a collision **'
      write(6,*)
      write(6,*)'    iprod= ',iprod
      write(6,*)
      write(6,*)'  ** R1: initial (e,v,j,Omg): '
     &     ,ielecref,nvref,jref,iomref
      write(6,*)'        Jtot,il= ',Jtot,il0
      write(6,*)'  **  R2 initial conditions **'
      write(6,*)'      colision energy(eV)= ',ecol_ev
      write(6,*)'      energy width(eV)= ',deltaE_ev
      write(6,*)'  r20 =',rcolini,'   k= ',xk0,'   alpha= ',alpha0
      write(6,*)
      if(npun1.eq.1.and.iprod.eq.1)then
          write(6,*)' for rigid rotor (npun1=1) iprod must be 0'
          call flush(6)
          stop
      endif
      xxx=0.d0
      do ir2=1,npun2
         r2=rmis2+dble(ir2-1)*ah2
         call  rgauscolini(xr2,r2,rcolini,alpha0,xk0,factor0,il0)
         xxx=xxx+xr2*xr2
      enddo
      xxx=xxx*ah2
      write(6,*)'  norm of incident gaussian in R2= ',xxx

      if(iprod.eq.1)then
         call colini_prodJacobi
      else
         call colini_reacJacobi
      endif

**> checking norm

      xxx=0.0d0
      do i=1,ntotproc(idproc)
          xxx=xxx+rpaqproc(i)*rpaqproc(i)
       enddo
       write(6,*)' initial norm proc=',idproc,' is ',xxx
      call MPI_REDUCE(xxx,yyy,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
      if(idproc.eq.0)then
         write(6,*)
         write(6,*)'  ** Norm of initial wavepacket= ',yyy
         call flush(6)
      endif
      
      return
      end subroutine set_colini
!--------------------------------------------------
      subroutine colini_reacJacobi
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      implicit none
      integer :: i,ican,icanp,ielec,iom,iang,iangp,ir,ir1,ir2
      double precision :: r1,r2,rp,Rg,angpart,fdiatom
      double precision :: xr2,xxx
      include "mpif.h"
 
      rpaqproc(:)=0.d0
      write(6,*)'   Initial wave packet in Reactants Jacobi coordinates'

!      write(name,"('funini.id',i2.2)")idproc
!      open(68,file=name,status='unknown')
      do i=1,ntotproc(idproc)
         call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
         ican=ibasproc(icanp,idproc)
         iang=indangreal(iangp,idproc)
         if(iom.eq.iomref.and.ielec.eq.ielecref)then
            r2=rmis2+dble(ir2-1)*ah2
            r1=rmis1+dble(ir1-1)*ah1
            Rg=R2

            angpart=Djmmp(iang,jref,ielec,iom)
            if(npun1.eq.1)then
               fdiatom=1.d0
            else
               fdiatom=fd(ir1,nvref,jref,ielecref)
            endif

            rp=r1
            call  rgauscolini(xr2,rg,rcolini,alpha0,xk0,factor0,il0)
            xxx=-xr2*angpart
            rpaqproc(i)=xxx*fdiatom*dsqrt(ah2)
!            write(68,*)r1,r2,iang,fdiatom,xr2,angpart,iom,ielec
!     &                  ,i,rpaqproc(i)
         endif
      enddo
!      close(68)
      
      return
      end subroutine colini_reacJacobi 
!--------------------------------------------------
      subroutine colini_prodJacobi
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      
      implicit none
      include "mpif.h"
      double precision , allocatable :: ff(:,:),xx(:)
      double precision , allocatable :: popmin(:),poptot(:),pop(:)
      double precision , allocatable :: fr1(:),fr2(:),fang(:)
      double precision , allocatable :: distri(:,:),distritot(:,:)
      double precision :: A,B,A1,r1,r2,cos,xxx,cbeta,ssign,x1
      double precision :: beta,dd,Rg,rp,angjac,angpart,fdiatom
      double precision :: rpaqref,ss,xn,xr2
      integer :: ierr,ir1,ir2,jtot2,itot,icanp,iangp,ir,ican,ielec
      integer :: iom,m2,mp2,iang,mmm2,jref2,ii1,ii2,nsend
      integer :: ierror,iold,mmmp2
      
!***>  building initial wavepacket in product Jacobi coordinates
!*     *> transforming from reactant to product Jacobi coordinates
      
      write(6,*)' Starting collision in products channel:'
      write(6,*)'                    tranforming initial wvp!!!'
      write(6,*)'       iprod',iprod

      rpaqproc(:)=0.d0
! allocating variables 

      allocate(ff(npun2,2),xx(npun2),pop(iommin:iommax)
     &     ,poptot(iommin:iommax),fr1(npun1),fr2(npun2),fang(nangu)
     &     ,distri(npun1,npun2),distritot(npun1,npun2)
     &     , stat=ierror)

      fr1(:)=0.d0
      fr2(:)=0.d0
      fang(:)=0.d0
      distri(:,:)=0.d0
      distritot(:,:)=0.d0
      
! calculating initial wavepacket in products Jacobi coordinates
          do ir2=1,npun2
            xx(ir2)=rmis2+dble(ir2-1)*ah2
            ff(ir2,1)=fdprod(ir2,nvref,jref,ielecref)*dsqrt(ah1/ah2)
            ff(ir2,2)=0.d0
            write(67,*)xx(ir2),ff(ir2,1)
         enddo
         call splset(ff,xx,npun2,npun2)

         jtot2=jtot*2
         B=(xm0+xm1)/xm1-xm2/(xm0+xm2)
         A=(xm0+xm2)/xm2-xm1/(xm0+xm1)
         A1=xm1/(xm0+xm1)
         if(idproc.eq.0)write(6,*)'  A,A1,B= ',A,A1,B

         do itot=1,ntotproc(idproc)
            icanp=indcanproc(itot)
            iangp= indangproc(itot)
            ir=indrgridproc(itot)

            ican=ibasproc(icanp,idproc)
            ielec=nelebas(ican)
            iom=iombas(ican)
            m2=2*iom
            mp2=2*iomref
            if(ielec.eq.ielecref)then
               iang=indangreal(iangp,idproc)
               ir1=indr1(ir,iang)
               ir2=indr2(ir,iang)
               r1=rmis1+dble(ir1-1)*ah1
               r2=rmis2+dble(ir2-1)*ah2
               cos=cgamma(iang)          
               xxx=A*r1*cos
               cbeta=(xxx-r2)/dsqrt(r1*r1*a*a+r2*r2-2.d0*xxx*r2)
               ssign=(-1.d0)**(iabs(iom))
               x1=0.d0
               if(dabs(cbeta).le.1.d0)then
                  beta=dacos(cbeta)
                  call dwigner(pm,jtot2,mp2,m2,beta,ndim)
                  x1=pm(Jtot2)
                  if(iomref.ne.0)then
                    call dwigner(pm,jtot2,-mp2,m2,beta,ndim)
                    ss=dble(iparity)*((-1.d0)**(iabs(Jtot+iomref)))
                    x1=x1+ss*pm(Jtot2)
                    x1=x1/dsqrt(2.d0)
                  endif

                  if(iom.ne.0)x1=x1*dsqrt(2.d0)

**building wvp in internal variables of Product Jacobi coordinates

                  Rg=A*A*r1*r1+R2*R2-2.d0*A*r1*r2*cgamma(iang)
                  Rg=dsqrt(Rg)*xm2/(xm0+xm2)

                  rp=A1*a1*r1*r1+R2*r2+2.d0*a1*r1*r2*cgamma(iang)
                  rp=dsqrt(rp)
                  angjac=a*a1*r1*r1+(a-a1)*r1*r2*cgamma(iang)-r2*r2
                  angjac=angjac*xm2/(rp*Rg*(xm0+xm2))
                  angjac=dacos(angjac)

                  if(angjac.ge.0.d0.and.angjac.le.pi)then
                     mmm2=2*(iomref-iomatom(ielecref))
                     mmmp2=2*iomdiat(ielecref)
                     jref2=2*jref
                     call dwigner(pm,jref2,mmm2,mmmp2,angjac,ndim)
                     xn=dsqrt((2.d0*dble(jref)+1.d0)*0.5d0)
                     angpart=pm(jref2)*xn*dsqrt(wreal(iang))
                  else
                     angpart=0.d0
                  endif

                  if(npun1.eq.1)then
                     fdiatom=1.d0
                  else
                     iold=2
                     if(rp.gt.xx(1).and.rp.le.xx(npun2))then
                        call splinqq(ff,xx,iold,npun2,rp,npun2,fdiatom)
!                        fdiatom=splinq(ff,xx,iold,npun2,rp,npun2)
                     else
                        fdiatom=0.d0
                     endif
                  endif

               call  rgauscolini(xr2,rg,rcolini,alpha0,xk0,factor0,il0)

                  xxx=-xr2*angpart*r2*r1/(Rg*rp)
  
                  rpaqref=xxx*fdiatom*dsqrt(ah2)
***end paqref
                  rpaqproc(itot)=rpaqproc(itot)+ssign*x1*rpaqref
                  distri(ir1,ir2)=distri(ir1,ir2)+ssign*x1*rpaqref
                  ii1=int( (rp-rmis1)/ah1 +1)
                  if(ii1.ge.1.and.ii1.le.npun1)then
                     fr1(ii1)=fr1(ii1)+fdiatom*fdiatom
                  endif
                  ii2=int( (Rg-rmis2)/ah2 +1)
                  if(ii2.ge.1.and.ii2.le.npun2)then
                     fr2(ii2)=fr2(ii2)+xr2*xr2
                  endif
                  fang(iang)=fang(iang)+ssign*x1*rpaqref
c                     write(6,'(11(1x,e15.7))')r1,r2,rp,Rg
c     &                   ,dacos(cgamma(iang))*180.d0/pi,angjac*180/pi
c     &                   ,fdiatom,xr2,angpart,xxx,paqref
               endif  ! cbeta.le.1
            endif ! elec=ielecref 
         enddo  !  itotp

         nsend=npun1*npun2
         call MPI_REDUCE(distri,distritot,nsend,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
       
            write(name,'("distri.",i3.3)')idproc
            open(33,file=name,status='unknown')
            do ir1=1,npun1
               do ir2=1,npun2
                  r1=rmis1+dble(ir1-1)*ah1
                  r2=rmis2+dble(ir2-1)*ah2
                  write(33,*)r1,r2,distritot(ir1,ir2)
               enddo
               write(33,'()')
            enddo
            close(33)

         do iom=iommin,iommax
            pop(iom)=0.d0
         enddo
         do itot=1,ntotproc(idproc)
            call indiproc(itot,icanp,ielec,iom,iangp,ir,ir1,ir2)
            pop(iom)=pop(iom)+rpaqproc(itot)*rpaqproc(itot)
         enddo
!         do iom=iommin,iommax
!            write(6,*)idproc,iom,pop(iom)
!         enddo

         nsend=iommax-iommin+1
         call MPI_REDUCE(pop,poptot,nsend,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)


         if(idproc.eq.0)then
            write(6,*)'  Helicity distribution in product Jacobi coor.'
            write(6,*)
            do iom=iommin,iommax
               write(6,*)iom,poptot(iom)
            enddo
            write(6,*)
         endif
*     deallocating auxiliary variables
         
      deallocate(ff,xx,pop,poptot,fr1,fr2,fang
     &     ,distri,distritot)
         
      return
      end subroutine colini_prodJacobi 
!--------------------------------------------------
!=======================================================
      end module mod_colini_01y2
