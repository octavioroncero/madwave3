      program Aeinstein

! *********************************************************************
!*     *  Calculate Einstein coefficients between bound states
!     *     *     of a triatomic molecule
!  from state A  --> stateX
!*     *       Using  Jacobi coordinates                            **
!*               Program works internally in a.u.                    **
! *********************************************************************

      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      use mod_photoini_01y2

      implicit none
      include "mpif.h"     
      integer ::  nv1_fin,nv2_fin,Jtot_fin,ipar_fin,iv_fin,iv_ini
      double precision ::  EShift_fin_eV, Energy_fin_eV
      double precision ::  EShift_fin, Energy_fin, Energy_ini
      double precision ::  Ephoton_eV,Ephoton,Atot,autocor,xmat,xmat_tot
      integer :: ierror,i
      integer :: icanp,ielec,iom,iangp,ir,ir1,ir2,iang
      character*50 :: filename

      double precision, allocatable :: Efin(:),Apartial(:),funA(:,:,:,:)

!     initial state is read in dip_bnd routine of module mod_photoini_01y2.f
!     in directory ../bnd
!     using the namelist inputbnd with data
!      namelist /inputbnd/Jtotini,iparini,nvbound,nprocbnd,maxbnddim
!     &                  ,igridbnd
!     whose initial energy is energy_ini_eV
!
      
*********************************************************
      namelist /einstein/nv1_fin,nv2_fin,EShift_fin_eV
*********************************************************

c! Initialize MPI environment and get proc's ID and number of proc in
c! the partition.

      call MPI_INIT(ierror)
!      call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)

      write(filename,'("sal."i3.3)')idproc
      open(6,file=filename,status='unknown')

!     determining grids and costants
      
      call input_grid
      call pot0
      call  paralelizacion
      call basis
      
!     reading potential

      call pot2
      Jtot_fin=Jtot
      ipar_fin=iparity
!     setting rpaq and rpaq0 dimensions

      allocate(rpaq0(nbastotproc)
     &    , funA(npun1,npun2,nangu2,0:Jtot_fin)
     &     )

      rpaq0(:)=0.d0
      funA(:,:,:,:)=0.d0

!     reading data for rovibroelectronic transitions
      
      write(6,*)
      write(6,*)'  AEinstein data'
      write(6,*)'  -------------------'
         open(10,file='input.dat',status='old')
         read(10,nml =einstein)
         write(6,nml =einstein)
         call flush(6)
         close(10)
         EShift_fin=Eshift_fin_eV/au2eV
       
      allocate(Efin(nv1_fin:nv2_fin),Apartial(nv1_fin:nv2_fin)
     &       , stat=ierror)
      

!     reading transition dipole moments
      
      call read_trans_dipole
      
!     reading data for initial state: d.e | Psi^Jini > of energy= energy_ini_eV
      
      call dip_bnd
      energy_ini=energy_ini_eV/au2eV
      iv_ini=nvbound
      
!     loop over final vibrational states
      
      Atot=0.d0
      xmat_tot=0.d0
      do iv_fin=nv1_fin,nv2_fin
         call read_bndgridA(funA
     &                 ,iv_fin,Jtot_fin,ipar_fin
     &        ,energy_fin_eV)
         
         energy_fin=energy_fin_eV/au2eV
         Ephoton_eV=energy_fin_eV- EShift_fin_eV-Energy_ini_eV
         Ephoton=Ephoton_ev/au2eV


**>> electric dipole matrix elements

         autocor=0.d0
         do i=1,ntotproc(idproc)
            call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
            iang=indangreal(iangp,idproc)
            autocor=autocor+rpaq0(i)*funA(ir1,ir2,iang,iom)
         enddo
         call MPI_REDUCE(autocor,xmat,1,MPI_REAL8,MPI_SUM
     &        ,0,MPI_COMM_WORLD,ierror)

         write(6,*)xmat*xmat,' overlap**2 with renormalized d.e Psi'
         xmat_tot=xmat_tot+xmat*xmat
         write(6,*)energy_ini_eV,energy_fin_eV,Eshift_fin_eV,Ephoton_eV
     &    ,'  E_i   E_f Shift Photon (eV)'
         
         Apartial(iv_fin)=dabs(xmat*xmat*xnorm_ini
     &                     *Ephoton*Ephoton*Ephoton *Aconstant_au)
         if(idproc.eq.0)then
         write(6,'(2(1x,i3),2(1x,d15.7)
     &              ," v_i  v_f  A_{vi vf} E_photon")')
     &           iv_ini,iv_fin,Apartial(iv_fin),Ephoton_eV
         endif
         Atot=Atot+Apartial(iv_fin)

       enddo  ! ivX
         
       if(idproc.eq.0)then
          write(6,*)xmat_tot,' total overlap^2 with renorm. d.e Psi'
      
          write(6,*)' v= ',iv_ini,' total rate(a.u.)= ',Atot
     &          ,' lifetime(s)=',2.418884396d-17/Atot
       endif                                        

      stop
      end


!----------------------------------------

      subroutine read_bndgridA(funA,nvbound_fin,Jtot_fin,ipar_fin
     &                       ,energy_fin_eV)
      
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      use mod_photoini_01y2
      implicit none
      include "mpif.h"
      character*50 filegrid,filebnd

      integer :: ielecbnd,iombnd
     &                     ,ir1bnd,ir2bnd,iangbnd
      real*8 :: bndvec
            integer :: ielecbndmin,iombndmin,ibnd,ijk,iiicanp,iii
      integer :: iiir,iangp,ierror,i,iang,icanp,ielec,ierr
      integer :: iprocbnd,nbnddim,npun1bnd,npun2bnd,incbnd
      integer :: iom,iomini,iq,ir,ir1,ir2,nangubnd,j2,ifail
      integer :: ivX,it,it0,indt,iloop,nnn
      real*8 :: rmis1bnd,rfin1bnd,rmis2bnd,rfin2bnd
      real*8 :: ah1bnd,ah2bnd,x123,xnorm,xnormtot,y123,ssign
      real*8 :: yyy,xxx,coef,coef2,xx,yy,energy_fin_eV
      real*8 :: xmat,autocor
      integer :: nvbound_fin,Jtot_fin,ipar_fin
      double precision :: funA(npun1,npun2,nangu2,0:Jtot_fin)

!      allocate(ielecbnd(maxbnddim),iombnd(maxbnddim)
!     &     ,ir1bnd(maxbnddim),ir2bnd(maxbnddim),iangbnd(maxbnddim)
!     &     ,bndvec(maxbnddim)
!     &     , stat=ierror)

*     reading bnd
      write(6,*)'--------------------------------'
      write(6,*)'   State A kept in ../bndAJXYp '
      call flush(6)
         write(6,*)'   reading bound nvbound= ',nvbound_fin
     &                             ,' state from ',nprocbnd
     &                                    ,' files in bnd dir'
          write(6,*)'      calculated with boundlanz-madwave.f '
          write(6,*)'           in the same grid!!!!!'
          write(6,*)' '

**>> reading grid of bnd calculation for Jtotini
*    to adapt it to the dissociation dynamics with Jtot and nelecmax

          funA(:,:,:,:)=0.d0
          xnorm=0.d0
         do iprocbnd=0,nprocbnd-1
            if(ipar_fin.eq.1)then
               write(filebnd
     &         ,'("../bndAJ",i2.2,"p/bnd.iv",i3.3,".ip",i3.3)')
     &                             Jtot_fin,nvbound_fin,iprocbnd

               write(filegrid,'("../bndAJ",i2.2,"p/grid.ip",i3.3)')
     &                           Jtot_fin,iprocbnd
            elseif(ipar_fin.eq.-1)then
               write(filebnd
     &         ,'("../bndAJ",i2.2,"m/bnd.iv",i3.3,".ip",i3.3)')
     &                             Jtot_fin,nvbound_fin,iprocbnd

               write(filegrid,'("../bndAJ",i2.2,"m/grid.ip",i3.3)')
     &                           Jtot_fin,iprocbnd
            else
               write(6,*)' no good ipar_fin =',ipar_fin
               call flush(6)
               call MPI_BARRIER(MPI_COMM_WORLD, ierror)
               stop
            endif

            open(111,file=filegrid,status='old')

            open(110,file=filebnd,status='old')

            read(110,*)ivX,energy_fin_eV
         
            read(111,*)nbnddim
            write(6,*)' Reading ',nbnddim,' values in ',filebnd
     &        ,' & ',filegrid
            call flush(6)

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


c-----> end checking grids

            do ibnd=1,nbnddim
               read(111,*)iii,iiicanp,ielecbnd,iombnd,iangp
     &            ,iiir,ir1bnd,ir2bnd,iangbnd
               read(110,*)xxx
               ir1=ir1bnd
               ir2=ir2bnd
               iang=iangbnd
               iom=iombnd
               funA(ir1,ir2,iang,iom)=xxx
!               write(69,*)ir1,ir2,iang,iom,xxx
               xnorm=xnorm+xxx*xxx
            enddo

            close(110)
            close(111)
         
         enddo                     ! iprocbnd

         write(6,*)'     norm= ',xnorm
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)

 !     deallocate(ielecbnd,iombnd
 !    &     ,ir1bnd,ir2bnd,iangbnd
 !    &     ,bndvec)
      return

      end subroutine read_bndgridA
