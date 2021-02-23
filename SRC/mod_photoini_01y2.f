       module mod_photoini_01y2
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

      integer :: Jtotini,iparini,nvbound,nprocbnd,maxbnddim
      integer :: ifileauto,igridbnd
      real*8, allocatable :: dipol(:,:,:,:)
      real*8 :: energy_ini_eV,xnorm_ini
      real*8,allocatable :: factresj(:,:,:)
      
      contains
********************************************
*   functions of   mod_photoini_01y2.f     *
********************************************
!=======================================================

!--------------------------------------------------
      subroutine set_trans_dipole
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      implicit none
      include "mpif.h"
      
      integer :: iang_proc,iang,ir,ir1,ir2,ielec,iq,ierror
      real*8 :: r1,r2,cgam
      
      real*8 :: dx(nelecmax),dy(nelecmax),dz(nelecmax)
      real*8 :: dip(-1:1,nelecmax)

      write(6,*)
      write(6,*)'   --------------------------------'
      write(6,*)'   Setting transition dipole moment  '
      write(6,*)'   --------------------------------'
      write(6,*)
      call flush(6)
      
      call setdipini
      
      do iang_proc=1,nanguprocdim
         iang=indangreal(iang_proc,idproc)
         write(name,'("dip/dip.",i3.3,".dat")')iang
         write(6,*)' in set_trans_dipole ',idproc,iang_proc,iang,name
         open(10,file=name,status='new')
         do ir=1,npunreal(iang)
            ir1=indr1(ir,iang)
            ir2=indr2(ir,iang)
            r1=rmis1+dble(ir1-1)*ah1
            r2=rmis2+dble(ir2-1)*ah2
            cgam=cgamma(iang)
            call dipele(r1,r2,cgam,dx,dy,dz,nelecmax,nelecmax)

            do ielec=1,nelecmax
               if(dabs(dx(ielec)).gt.1.d-10
     &           .or.dabs(dz(ielec)).gt.1.d-10)then
                  dip(0,ielec)=dz(ielec)
                  dip(1,ielec)=-dx(ielec)/dsqrt(2.d0)
                  dip(-1,ielec)=dx(ielec)/dsqrt(2.d0)
               elseif(dabs(dy(ielec)).gt.1.d-10)then
                  dip(0,ielec)=0.d0
                  dip(1,ielec)=dy(ielec)/dsqrt(2.d0)
                  dip(-1,ielec)=dy(ielec)/dsqrt(2.d0)
               else
                  dip(0,ielec)=0.d0
                  dip(1,ielec)=0.d0
                  dip(-1,ielec)=0.d0
               endif
            enddo

            write(10,*)ir1,ir2
     &              ,((dip(iq,ielec),ielec=1,nelecmax),iq=-1,1)
         enddo                  ! ir
         close(10)
      enddo  ! iang_proc

      write(6,*)'  transition electric dipole files written'
      call flush(6)
    
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      write(6,*)ierror
      
      return
      end subroutine set_trans_dipole
!--------------------------------------------------
      subroutine read_trans_dipole
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      implicit none
      include "mpif.h"
     
      integer :: ierror
      integer :: iang_proc,iang,ir,ir1,ir2,ielec,iq,nsend
      real*8 :: r1,r2,cgam,xxx
      
      allocate(dipol(npuntot,nangu,nelecmax,-1:1)
     &       , stat=ierror)
      norealproc_mem=norealproc_mem
     &    +npuntot*nangu*nelecmax*3
      write(6,*)'norealproc_mem= ',norealproc_mem
      write(6,*)'nointegerproc_mem= ',nointegerproc_mem
      call flush(6)
      
      write(6,*)
      write(6,*)'   --------------------------------'
      write(6,*)'   Reading transition dipole moment  '

      do iang=1,nangu
         write(name,'("../dip/dip.",i3.3,".dat")')iang
         open(10,file=name,status='old')
         
         do ir=1,npunreal(iang)
            ir1=indr1(ir,iang)
            ir2=indr2(ir,iang)
            read(10,*)ir1,ir2
     &   ,((dipol(ir,iang,ielec,iq),ielec=1,nelecmax),iq=-1,1)
         enddo                  ! ir
         close(10)
      enddo                     ! iang
      
      write(6,*)
      write(6,*)'    transition dipole moment READ  '
      write(6,*)'   --------------------------------'
      write(6,*)

      return
      end subroutine read_trans_dipole
!--------------------------------------------------
      subroutine dip_bnd
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      implicit none
      include "mpif.h"
      character*50 :: filegrid,filebnd

      real*8, allocatable :: xnormOmg(:),xnormOmgProc(:)
      integer :: ierror,j2,ifail,iom,ijk,iomini,iq,i,it
      integer :: icanp,ielec,iangp,ir,ir1,ir2,ierr,nnn,it0,indt,iloop
      real *8 :: ssign,yyy,coef,coef2,xxx,xnorm,xnormtot,xx,yy

*********************************************************
      namelist /inputbnd/Jtotini,iparini,nvbound,nprocbnd,maxbnddim
     &                  ,igridbnd

      igridbnd=1
      write(6,'(40("="),/,10x,"dip_bndgrid",/,40("="))')
      write(6,*)
      write(6,*)'  grid and basis data'
      write(6,*)'  -------------------'
      open(10,file='input.dat',status='old')
      read(10,nml = inputbnd)
      write(6,nml = inputbnd)
      call flush(6)
      close(10)

      allocate(factresj(0:Jtotini,-1:1,iommin:iommax)
     &     ,xnormOmg(iommin:iommax),xnormOmgProc(iommin:iommax)
     &     , stat=ierror)

**>> tresj factor

       write(6,*)
       write(6,*)'     ** Initial wavepacket **'
       if(iphoto.eq.1)then
         write(6,*)'        Phi^J_v(t=0)=< JM | d . e | Psi^Jini_v >'
         write(6,*)'     considering linearly polarized light p=0 '
         write(6,*)'     F_{MM_i} included later in treating spectrum'
         write(6,*)'          in program cheby spectra'
         write(6,'(/,10x,"3j coeff. for transition operator",/)')
         write(6,'(10x,"iom_i",7x,"q",7x,"iom",/)')
         factresj(:,:,:)=0.d0
         j2=1
         ifail=0
         do iom=iommin,iommax
            ijk=0
            call factorial
         if(iparini*((-1.d0)**(Jtotini)).lt.0.d0)ijk=1
         do iomini=ijk,Jtotini
         do iq=-1,1
            ssign=(-1.d0)**dble(iom+iq)
            yyy=1.d0
            if(iom.eq.0)yyy=yyy/dsqrt(2.d0)
            if(iomini.eq.0)yyy=yyy/dsqrt(2.d0)
            if(iparity*iparini.ge.0)yyy=0.d0

            call tresj(Jtotini,j2,Jtot,iomini,iq,-iom,coef)
            call tresj(Jtotini,j2,Jtot,-iomini,iq,-iom,coef2)
            xxx=coef+coef2*dble(iparini*((-1)**Jtotini))
            factresj(iomini,iq,iom)=ssign*yyy*xxx
            if(dabs(factresj(iomini,iq,iom)).gt.1.d-10)ifail=ifail+1
         write(6,*)iomini,iq,iom,factresj(iomini,iq,iom),ssign,yyy,xxx
         enddo
         enddo
         enddo
         call flush(6)

         if(ifail.eq.0)then
           write(6,*)' Forbidden transition for iphoto=1 '
           write(6,*)' Check rotational quantum numbers'
           call flush(6)
           call MPI_BARRIER(MPI_COMM_WORLD, ierror)
           stop
         endif
      
      elseif(iphoto.eq.2)then
         write(6,*)'        Phi^J_v(t=0)= Psi^Jini_v >'

         if(Jtot.ne.Jtotini.or.iparity.ne.iparini)then
           write(6,*)' Forbidden transition for iphoto=2'
           write(6,*)' jtot, iparity should be Jtotini,iparini'
           call flush(6)
           call MPI_BARRIER(MPI_COMM_WORLD, ierror)
           stop
         endif
      elseif(iphoto.gt.2)then
         write(6,*)' iphoto =',iphoto,' is out of range!'
         write(6,*)'       Change it to 0,1 or 2'
         call flush(6)
         call MPI_BARRIER(MPI_COMM_WORLD, ierror)
         stop
      endif                     ! iphoto=1

      write(6,*)'              J,   ipar= ',Jtot,iparity
      write(6,*)'           Jini,iparini= ',Jtotini,iparini
      write(6,*)'              v= ',nvbound
      write(6,*)' '

      if(igridbnd.eq.1)then
           call dip_bndgrid
      else
           call dip_bndbase
      endif
  
**> checking norm

      xnorm=0.d0
      xnormOmgProc(:)=0.d0
      do i=1,ntotproc(idproc)
         call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
         xnorm=xnorm+rpaq0(i)*rpaq0(i)
         xnormOmgProc(iom)=xnormOmgProc(iom)+rpaq0(i)*rpaq0(i)
      enddo
      write(6,*)' xnorm partial = ',xnorm

      call MPI_REDUCE(xnorm,xnormtot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xnormtot,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      if(xnormtot.lt.1.d-15)then
         write(6,*)' norm of initial wvp too low= ',xnormtot
         call flush(6)
         call MPI_BARRIER(MPI_COMM_WORLD, ierror)
         stop
      endif

      xnorm_ini=xnormtot
      do i=1,ntotproc(idproc)
         rpaq0(i)=rpaq0(i)/dsqrt(xnormtot)
      enddo
      
      write(6,*)
      if(iphoto.eq.1)then
      write(6,*)'   iv= ',nvbound,'  |d.e Psi|^2= ',xnormtot
      write(6,*)'         renormalizing...'
      else
      write(6,*)'   iv= ',nvbound,'  | Psi|^2= ',xnormtot
      endif
      nnn=iommax-iommin+1
      call MPI_REDUCE(xnormOmgproc,xnormOmg,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xnormOmg,nnn,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      do iom=iommin,iommax
         write(6,*)' Omg= ',iom,'  Prob= ',xnormOmg(iom)/dsqrt(xnormtot)
      enddo
      write(6,*)      

      call flush(6)
      
!     preparing autocorrelation file

      ifileauto=11
      write(name,'("auto.v",i2.2,".j",i2.2,".J",i3.3)')nvref,jref,Jtot
      it0=0
      indt=0
      iloop=0

      open(3,file='cont.data',status='old',err=1)
      read(3,*)it0,indt,iloop
      close(3)
 1    continue
      
      if(idproc.eq.0)then
         if(it0.eq.0)then
            open(ifileauto,file=name,status='new')
            write(ifileauto,'(4(1x,d20.12))')xnorm_ini,energy_ini_eV
     &           ,emindlt,delta2
c            write(ifileauto,*)0,1.d0
            call flush(ifileauto)
         else
            open(ifileauto,file=name,status='old')
            read(ifileauto,*)xx
            do it=1,it0
               read(ifileauto,*)xx,yy
            enddo
         endif
      endif
      
!     deallocate auxiliar matrices
      
      deallocate(factresj,xnormOmg,xnormOmgProc)

  
      return
      end subroutine dip_bnd
!--------------------------------------------------
!--------------------------------------------------
      subroutine dip_bndgrid
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      implicit none
      include "mpif.h"
      character*50 :: filegrid,filebnd
      integer, allocatable:: ielecbnd(:),iombnd(:)
     &                     ,ir1bnd(:),ir2bnd(:),iangbnd(:)
      real*8,allocatable :: bndvec(:)
 
      integer :: ielecbndmin,iombndmin,ibnd,ijk,iiicanp,iii
      integer :: iiir,iangp,ierror,i,iang,icanp,ielec,ierr
      integer :: iprocbnd,nbnddim,npun1bnd,npun2bnd,incbnd
      integer :: iom,iomini,iq,ir,ir1,ir2,nangubnd,j2,ifail
      integer :: ivX,it,it0,indt,iloop,nnn
      real*8 :: rmis1bnd,rfin1bnd,rmis2bnd,rfin2bnd
      real*8 :: ah1bnd,ah2bnd,x123,xnorm,xnormtot,y123,ssign
      real*8 :: yyy,xxx,coef,coef2,xx,yy

      allocate(ielecbnd(maxbnddim),iombnd(maxbnddim)
     &     ,ir1bnd(maxbnddim),ir2bnd(maxbnddim),iangbnd(maxbnddim)
     &     ,bndvec(maxbnddim)
     &     , stat=ierror)

* reading bnd 
         write(6,*)'   reading bound nvbound= ',nvbound
     &                             ,' state from ',nprocbnd
     &                                    ,' files in bnd dir'
          write(6,*)'      calculated with boundlanz-madwave.f '
          write(6,*)'           in the same grid!!!!!'
          write(6,*)' '

**>> reading grid of bnd calculation for Jtotini
*    to adapt it to the dissociation dynamics with Jtot and nelecmax

         do iprocbnd=0,nprocbnd-1
            if(iparini.eq.1)then
               write(filebnd
     &         ,'("../../bndJ",i1,"p/bnd.iv",i3.3,".ip",i3.3)')
     &                             Jtotini,nvbound,iprocbnd

               write(filegrid,'("../../bndJ",i1,"p/grid.ip",i3.3)')
     &                           Jtotini,iprocbnd
            elseif(iparini.eq.-1)then
               write(filebnd
     &         ,'("../../bndJ",i1,"m/bnd.iv",i3.3,".ip",i3.3)')
     &                             Jtotini,nvbound,iprocbnd

               write(filegrid,'("../../bndJ",i1,"m/grid.ip",i3.3)')
     &                           Jtotini,iprocbnd
            else
               write(6,*)' no good iparini =',iparini
               call flush(6)
               call MPI_BARRIER(MPI_COMM_WORLD, ierror)
               stop
            endif


            open(111,file=filegrid,status='old')

            open(110,file=filebnd,status='old')

            read(110,*)ivX,energy_ini_eV
         
            read(111,*)nbnddim
            write(6,*)' Reading ',nbnddim,' values in ',filebnd
     &        ,' & ',filegrid
            call flush(6)
            if(nbnddim.gt.maxbnddim)then
            write(6,*)'  nbnddim= ',nbnddim,'  > maxbnddim= ',maxbnddim
               call flush(6)
               call MPI_BARRIER(MPI_COMM_WORLD, ierror)
               stop
            endif

            read(111,*)rmis1bnd,rfin1bnd,npun1bnd
            ah1bnd=(rfin1bnd-rmis1bnd)/dble(npun1bnd-1)
            if(dabs(rmis1-rmis1bnd).gt.1.d-4
     &         .or.dabs(ah1-ah1bnd).gt.1.d-4)then
               write(6,*)'  grid in r1 for bnd non equal '
               write(6,*)'  rmis1bnd,ah1bnd= '
     &                    ,rmis1bnd,ah1bnd
               call flush(6)
               call MPI_BARRIER(MPI_COMM_WORLD, ierror)
               stop
            endif

            read(111,*)rmis2bnd,rfin2bnd,npun2bnd
            ah2bnd=(rfin2bnd-rmis2bnd)/dble(npun2bnd-1)
            if(dabs(rmis2-rmis2bnd).gt.1.d-4
     &         .or.dabs(ah2-ah2bnd).gt.1.d-4)then
               write(6,*)'  grid in r2 for bnd non equal '
               write(6,*)'  rmis2bnd,ah2bnd= '
     &                    ,rmis2bnd,ah2bnd
               call flush(6)
               call MPI_BARRIER(MPI_COMM_WORLD, ierror)
               stop
            endif

            read(111,*)nangubnd,incbnd
            if(nangu.ne.nangubnd.or.inc.ne.incbnd)then
               write(6,*)' grig in angle for bnd non equal '
               write(6,*)'  nangubnd,incbnd= ',nangubnd,incbnd
               call flush(6)
               call MPI_BARRIER(MPI_COMM_WORLD, ierror)
               stop
            endif

            if(nbnddim.gt.maxbnddim)then
            write(6,*)'  nbnddim= ',nbnddim,'  > maxbnddim= ',maxbnddim
               call flush(6)
               call MPI_BARRIER(MPI_COMM_WORLD, ierror)
               stop
            endif

c-----> end checking grids

            do ibnd=1,nbnddim
               read(111,*)iii,iiicanp,ielecbnd(ibnd),iombnd(ibnd),iangp
     &            ,iiir,ir1bnd(ibnd),ir2bnd(ibnd),iangbnd(ibnd)
            enddo
            do ibnd=1,nbnddim
               read(110,*)bndvec(ibnd)
            enddo

            close(110)
            close(111)

           write(6,*)' building initial wave packet for iphoto= ',iphoto
           call flush(6)

            do i=1,ntotproc(idproc)
               call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
               iang=indangreal(iangp,idproc)
               do ibnd=1,nbnddim

                  iomini=iombnd(ibnd)
                  iq=iom-iombnd(ibnd)
              
                  if(iphoto.eq.2)then                  
                     if(iang.eq.iangbnd(ibnd).and.iom.eq.iomini
     &                  .and.ielec.eq.ielecbnd(ibnd)
     &                  .and.ir1.eq.ir1bnd(ibnd)
     &                  .and.ir2.eq.ir2bnd(ibnd))then
  
                         rpaq0(i)=bndvec(ibnd)
                     endif                  
                  elseif(iphoto.eq.1)then
                     if(iang.eq.iangbnd(ibnd).and.iabs(iq).le.1
     &                .and.ielec.eq.ielecbnd(ibnd)
     &                .and.ir1.eq.ir1bnd(ibnd)
     &                .and.ir2.eq.ir2bnd(ibnd))then
c                  write(6,*)iang,ielec,ir1,ir2
c     &                 ,iangbnd(ibnd),ielecbnd(ibnd),ir1bnd(ibnd)
c     &                 ,ir2bnd(ibnd),iabs(iq)

                         y123=dipol(ir,iang,ielec,iq)
                         x123=factresj(iomini,iq,iom)
                         rpaq0(i)=rpaq0(i)+bndvec(ibnd)*y123*x123
                     endif
                  endif
               enddo ! ibnd
            enddo  ! i
            
         enddo                     ! iprocbnd

      call MPI_BARRIER(MPI_COMM_WORLD, ierror)

      deallocate(ielecbnd,iombnd
     &     ,ir1bnd,ir2bnd,iangbnd
     &     ,bndvec)
      return
      end subroutine dip_bndgrid
!--------------------------------------------------
!--------------------------------------------------
      subroutine dip_bndbase
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      implicit none
      include "mpif.h"

      character*40 :: frpini,fRgini,fcoef
      real*8, allocatable :: cdis(:,:),x(:),f(:,:)
      real*8, allocatable :: fdis(:,:),fr(:,:),frpdis(:,:),frp(:,:)
      real*8, allocatable :: faux(:,:),gaux(:,:),haux(:,:)
      real*8, allocatable :: hauxtot(:,:),Yjomgini(:,:,:),Psini(:,:,:)
      integer,allocatable :: jd(:),ld(:),nv1(:),nv2(:)
      integer :: ifile,nuno,npunlie1,nv1ini,nv1max
      integer :: npunlie2,nv2ini,nv2max,nbaslie,nnn,nntot,jjd
      integer :: j0ini,lmax,iang,iom,i,ir,iv2,iv1,ir2,ibas,nbnd,ir1
      real*8 :: energy_ini_cm,sum,popomg(0:Jtotini),y,xnorm,xang,xxx
      real*8 :: pop,fctr,y123,xx,yy,xnormtot,r1,r2
      integer :: iq,ierror,j,icanp,ielec,iangp,ican,iomini

      ifile=111
      nbnd=nvbound

      write(6,*)
      write(6,*)'     ** Initial bound state **'
      write(6,*)'        Phi(t=0)= d . e  Psi^Jini'
      write(6,*)'           Jini= ',Jtotini,' parity=',iparini
      write(6,*)'           v= ',nvbound
      write(6,*)'   bound state in basis with bndele program'
      write(6,*)
 
*             Numerical radial function in rp
      write(6,*)'   reading r functions for the wvp grid'
      call flush(6)
      frpini='../../bndele/rpwf1'
      write(6,"(/,5x,'Openining ifile= ',i3,' file: ',a40)")ifile
     &                                                      ,frpini
       open(ifile,file=frpini,status='old')
       read(ifile,*)npunlie1,nuno,nv1ini,nv1max
       close(ifile)
       open(ifile,file=frpini,status='old')
      allocate(x(npunlie1),f(npunlie1,2)
     &     ,frpdis(npunlie1,nv1ini:nv1max),frp(npun1,nv1ini:nv1max) 
     &     , stat=ierror)

       call rfunread(ifile,npunlie1,nv1ini,nv1max,frpdis,nuno,
     &              frp,npun1,rmis1,ah1,f,x)
       close(ifile)
       deallocate(x,f)

*             Numerical radial function in Rg

       write(6,*)'    reading R functions for the wvp grid'
       call flush(6)
       fRgini='../../bndele/Rgwf1'
       write(6,"(/,5x,'Openining ifile= ',i3,' file: ',a40)")ifile
     &                                                      ,fRgini
       open(ifile,file=fRgini,status='old')
       read(ifile,*)npunlie2,nuno,nv2ini,nv2max
       close(ifile)
       open(ifile,file=fRgini,status='old')
      allocate(x(npunlie2),f(npunlie2,2)
     &     ,fdis(npunlie2,nv2ini:nv2max),fr(npun2,nv2ini:nv2max) 
     &     , stat=ierror)

       call rfunread(ifile,npunlie2,nv2ini,nv2max,fdis,nuno,
     &              fr,npun2,rmis2,ah2,f,x)
      close(ifile)

*        Coefficient and basis set of the initial state

      if(iparini.gt.0)then
          write(fcoef,'("../../bndele/Xlie.J",i2.2,".p.j0")')Jtotini
      else
          write(fcoef,'("../../bndele/Xlie.J",i2.2,".m.j0")')Jtotini
      endif
      write(6,*)' Reading coefficients in file= ',fcoef
      call flush(6)
      open(ifile,file=fcoef,status='old')
      read(ifile,*)nnn
      read(ifile,*)nnn,nbaslie
      close(ifile)
      allocate(cdis(nbaslie,nvbound:nvbound)
     &     ,jd(nbaslie),ld(nbaslie),nv1(nbaslie),nv2(nbaslie)
     &     ,faux(nv1ini:nv1max,nv2ini:nv2max),gaux(nv1ini:nv1max,npun2)
     &     ,haux(npun1,npun2), hauxtot(npun1,npun2)
     &     , stat=ierror)
    
       nntot=nbaslie
       open(ifile,file=fcoef,status='old')
       write(6,*)ifile,nvbound,nvbound,nbaslie,nbaslie
       call flush(6)

       call coefread(ifile,nvbound,nvbound,nbaslie,nbaslie
     &              ,cdis,jd,ld,nv1,nv2,energy_ini_cm)

       close(ifile)
       energy_ini_eV=energy_ini_cm/8065.5d0

       if(nv1(nntot).gt.nv1max.or.nv2(nntot).gt.nv2max)then   
          write(6,*)' ** Attention in BNDCPL **'
          write(6,*)'   v1= ',nv1(nntot),'  while nv1max= ',nv1max
          write(6,*)'   v2= ',nv2(nntot),'  while nv2max= ',nv2max
          call flush(6)
          call MPI_BARRIER(MPI_COMM_WORLD, ierror)
          stop
       endif
       if(jd(nntot).gt.nangu2-1.and.nangu.gt.1)then
          write(6,*)' ** Attention in BNDCPL **'
          write(6,*)'   jd= ',jd(nntot),'  while nangu2-1= ',nangu2-1
          call flush(6)
c         stop
       endif
      
       write(6,*)
       write(6,*)'  Omega distribution'
       write(6,*)'-----------------------'
       do iom=0,Jtotini
          popomg(iom)=0.d0
       enddo
       sum=0.d0
       
       do i=1,nntot
          sum=sum+cdis(i,nbnd)*cdis(i,nbnd)
          popomg(ld(i))=popomg(ld(i))+cdis(i,nbnd)*cdis(i,nbnd)
        
       enddo
       do iom=0,Jtotini
          if(idproc.eq.0)write(6,*)iom,popomg(iom)
       enddo
       write(6,*)'     sum coef**2 = ',sum
       j0ini=jd(1)

***> building Spherical harmonics for Jtotini       
    
      allocate(Yjomgini(nangu,0:nangu*inc,0:jtotini)
     &        ,Psini(npuntot,nangu,0:Jtotini)
     &     , stat=ierror)

      write(6,*)' in bndbase jmax= ',jmax
      Psini(:,:,:)=0.d0
      do iom=0,Jtotini
         lmax=max0(jd(nntot),jmax)
         lmax=min0(lmax,ndim)
         do iang=1,nangu
            y=cgamma(iang)
            call nplegm(pm,lmax,iom,y,ndim)
            do j=0,lmax
               Yjomgini(iang,j,iom)=pm(j)*dsqrt(wreal(iang))
            enddo
         enddo
      enddo

***> building the initial state wavefunction

      write(6,*)' building Psi ini '
      call flush(6)

* initialization

      xnorm=0.d0
      steptot=ah1*ah2
      Psini(:,:,:)=0.d0
      do iom=0,Jtotini
      do iang=1,nangu
         faux(:,:)=0.d0
         gaux(:,:)=0.d0
 
* step 1: angular part

         do ibas=1,nntot
         if(jd(ibas).le.lmax.and.iom.eq.ld(ibas))then
            iv1=nv1(ibas)
            iv2=nv2(ibas)
            jjd=jd(ibas)
            xang=Yjomgini(iang,jjd,iom)*cdis(ibas,nbnd)
            faux(iv1,iv2)=faux(iv1,iv2)+xang
         endif
         enddo

* step 2: R2 part

         do iv2=nv2ini,nv2max
         do ir2=1,npun2
             xxx=fR(ir2,iv2)

             do iv1=nv1ini,nv1max
                gaux(iv1,ir2)= gaux(iv1,ir2)
     &                   +  faux(iv1,iv2)*xxx
             enddo
         enddo
         enddo

* step 3: R1 part and norm

         do ir=1,npunreal(iang)
            ir1=indr1(ir,iang)
            ir2=indr2(ir,iang)
            do iv1=nv1ini,nv1max
                xxx=frp(ir1,iv1)
                Psini(ir,iang,iom)=Psini(ir,iang,iom)
     &              + gaux(iv1,ir2)*xxx
             enddo
             xnorm=xnorm+Psini(ir,iang,iom)**2
         enddo

      enddo
      enddo

      write(6,*)' en building Psi ini'
      call flush(6)

      write(6,*)'  Omega dist. in initial state recons.'
      do iom=0,Jtotini
         pop=0.d0
         do iang=1,nangu
         do ir=1,npunreal(iang)
            pop=pop+  Psini(ir,iang,iom)**2
         enddo
         enddo
         write(6,*)iom,pop*steptot
      enddo

      write(6,*)'  j dist. in  Psini'
c      xnorm=0.d0
      do jjd=jini,jd(nntot)
         pop=0.d0
         haux(:,:)=0.d0

         do iom=0,Jtotini
         do iang=1,nangu
         do ir=1,npunreal(iang)
            ir1=indr1(ir,iang)
            ir2=indr2(ir,iang)
            xx=Psini(ir,iang,iom)*Yjomgini(iang,jjd,iom)
            haux(ir1,ir2)=haux(ir1,ir2)+xx
         enddo
         enddo
         enddo
  
         pop=0.d0
         do ir1=1,npun1
         do ir2=1,npun2
            yy=haux(ir1,ir2)
            pop=pop+yy*yy
         enddo
         enddo
         write(6,*)jjd,pop*steptot
c         xnorm=xnorm+pop
      enddo
 
      xnorm=xnorm*steptot
      xnormtot=xnorm

      write(6,*)
      write(6,*)'       norm of bound state= ',xnormtot
      write(6,*)

      if(dabs(xnormtot).lt.1.d-10)then
         write(6,*)' --> Norm too low!! <--'
         call flush(6)
         call MPI_BARRIER(MPI_COMM_WORLD, ierror)
         stop
      endif

***>  building initial wavepacket

      rpaq0(:)=0.d0
      do i=1,ntotproc(idproc)
         call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
         ican=ibasproc(icanp,idproc)
         iang=indangreal(iangp,idproc)

         if(iphoto.eq.1)then
            do iomini=0,Jtotini
               iq=iom-iomini
               fctr=factresj(iomini,iq,iom)*dsqrt(ah1*ah2)
               y123=fctr*dipol(ir,iang,ielec,iq)
               rpaq0(i)=rpaq0(i)+Psini(ir,iang,iomini)*y123
            enddo
            
         elseif(iphoto.eq.2)then

            if(Jtotini.ne.Jtot)then
               write(6,*)'  STOP!!: Jtotini= ',Jtotini
     &                  ,' must be equal to Jtot= ',Jtot
     &              ,' in the bright apprx. of iphoto= ',iphoto
               call flush(6)
               call MPI_BARRIER(MPI_COMM_WORLD, ierror)
               stop
            endif
                          
            rpaq0(i)=0.d0
            if(ielec.eq.ielecref)then
               do iomini=0,Jtotini
                  if(iomini.eq.iom)then
                     xx=Psini(ir,iang,iomini)*dsqrt(steptot)
                     rpaq0(i)=xx
                  endif
               enddo
            endif
            
         endif
      enddo

      deallocate (cdis,x,f,fdis,fr,frpdis,frp
     &     , faux,gaux,haux,hauxtot,Yjomgini,Psini
     &     , jd,ld,nv1,nv2 )

      return
      end subroutine dip_bndbase
!--------------------------------------------------
 !=======================================================
      end module mod_photoini_01y2
