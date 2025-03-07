      module mod_gridYpara_01y2
      implicit none 

!---------------------------------------------------------------------------------------!
! Determination of  grid for the AB + C  in a grid on R=R2, r=R1,gam                    !
!                                01 + 2
!
!  input_grid: determines grid 
!  paralelizacion: determines index matrices and distributes quantities
!_______________________________________________________________________________________!
      save
*     constants
      complex*16, parameter :: zero = dcmplx(0.d0,0.d0)
      complex*16, parameter :: zeye = dcmplx(0.d0,1.d0)
      real*8,parameter :: pi = dacos(-1.d0)
      real*8,parameter :: conve1 = 1.197d-4
      real*8,parameter :: convl = 0.52917726D0
      real*8,parameter :: convm = 0.182288853D4
      real*8,parameter :: conve = 4.55633538D-6
      real*8,parameter :: eV2cm = 8065.5d0
      real*8,parameter :: hbr = convl/dsqrt(convm*conve/conve1)
      real*8,parameter :: au2eV = 27.21113957d0
      real*8,parameter ::  zot2au=(1.d0/conve1)*conve

      real*8,parameter ::  cluz_au=137.d0
      real*8,parameter :: epsilon0_au=0.07957747d0 
      real*8,parameter :: Aconstant_au
     &              = 1.d0/(3.d0*pi*(cluz_au**3)*epsilon0_au) ! = 1/(3 pi hbar^4 Epsilon_0 c^3) 
       real*8,parameter :: CSconstant_au= 1.d0/(cluz_au*epsilon0_au)    ! = 1/(hbar^2 Epsilon_0 c) 
      
*     public ! for input data in namelist, grid, basis

*     grid and basis data in namelist
      integer :: npun1,npun1min,npun2,nangu,nangu2
      integer :: Jtot,iparity,inc,nelecmax,iommin,iommax,j0
      integer :: nvini,nvmax,jini,jmax
      integer :: nvref,jref,iomref,ielecref,ncan
      integer :: nproc,idproc

*     write options
      integer :: iwrt_pot,iwrt_wvp,iwrt_reac_distri
     &           ,n1plot,n2plot,nangplot

*     process

      integer :: iphoto
      real*8 :: photonorm

*     products states analysis

      integer :: iprod,nviniprod,nvmaxprod,jiniprod,jmaxprod
      integer :: iomminprod,iommaxprod,n2prod1,n2prod0
      integer :: nangproj0,nangproj1
      real*8 :: rbalinprod
      real*8 :: Emincut_prod,Emincut_prod_eV
      
*     for dimension of matrices, total quantites and per processor
      integer :: nomgproc,nanguproc,nanguprocdim,ncanmax
      integer :: nomgprocdim,nprocdim,nomgdim,ncanprocdim
      integer :: ncouprocmax
      
      real*8 :: rmis1,rfin1,rmis2,rfin2,ah1,ah2,steptot
      real*8, allocatable :: wreal(:),weight(:),cgamma(:)
     &                        ,cgamprod(:),weiprod(:)
      
      integer, allocatable ::  ncanproc(:),ibasproc(:,:)
     &              ,indiomreal(:,:),indangreal(:,:),indproc(:,:)
     &              ,ntotproc(:), ncouproc(:),ipcou(:,:)
      integer, allocatable :: iombas(:),nelebas(:)

*     to calculate total memory
      integer*8 :: nointeger_mem,noreal_mem
      integer*8 :: nointegerproc_mem,norealproc_mem
      contains
***************************************
*   functions of  mod_gridYpara_01y2  *
***************************************
!--------------------------------------------------
      subroutine input_grid
!--------------------------------------------------
      implicit none
      include "mpif.h"
      integer :: ierror,ir1,ir2,iang
      real*8 :: div
*********************************************************
      namelist /inputgridbase/npun1,rmis1,rfin1,npun1min
     &                       ,npun2,rmis2,rfin2
     &                       ,nangu
     &     ,Jtot,iparity,inc,nelecmax,iommin,iommax,j0
     &     ,jini,jmax,nvini,nvmax
     &                       ,nvref,jref,iomref,ielecref
*********************************************************
      namelist /inputprod/iprod
     &                   ,nviniprod,nvmaxprod
     &                   ,jiniprod,jmaxprod
     &     ,iomminprod,iommaxprod
     &     ,Rbalinprod,n2prod0,n2prod1,nangproj0,nangproj1
     &     ,Emincut_prod_eV
*********************************************************
      namelist /inputprocess/iphoto
*********************************************************
      namelist /inputwrite/iwrt_pot,iwrt_wvp,iwrt_reac_distri
     &                     ,n1plot,n2plot,nangplot              ! grid step on r1, r2, ang to print
      
      nointegerproc_mem=0
      norealproc_mem=0
      nointeger_mem=0
      noreal_mem=0
      
      npun1min=32
      iomref=0
      ielecref=1
      iommin=0
      
      write(6,'(40("="),/,10x,"GridYpara_mod",/,40("="))')
      write(6,*)
      write(6,*)'  grid and basis data'
      write(6,*)'  -------------------'
         open(10,file='input.dat',status='old')
         read(10,nml = inputgridbase)
         write(6,nml = inputgridbase)
         call flush(6)
         close(10)
         nangu2=nangu*inc
         iomminprod=0
         iommaxprod=0
         n2prod1=npun2
         n2prod0=1
         nangproj0=1
         nangproj1=nangu
         iprod=0
      write(6,*)
      write(6,*)'  products data'
      write(6,*)'  -------------'
         Emincut_prod_eV=0.d0
         open(10,file='input.dat',status='old')
         read(10,nml = inputprod)
         if(iommaxprod.eq.0)iommaxprod=min0(Jtot,jmaxprod)
         Emincut_prod=Emincut_prod_eV/conve1/eV2cm

         write(6,nml = inputprod)
         call flush(6)
         close(10)
         if(n2prod1.gt.npun2)then
            write(6,*)'  n2prod1= ',n2prod1,'  > npun2= ',npun2
            write(6,*)'   setting n2prod1=npun2'
            n2prod1=npun2
         endif

      write(6,*)
      write(6,*)'  process data'
      write(6,*)'  ------------'
         open(10,file='input.dat',status='old')
         read(10,nml = inputprocess)
         write(6,nml = inputprocess)
         call flush(6)
         close(10)
         
      iwrt_pot=0
      iwrt_wvp=0
      iwrt_reac_distri=0
      n1plot=1
      n2plot=1
      nangplot=1
      write(6,*)
      write(6,*)'  write data'
      write(6,*)'  ------------'
         open(10,file='input.dat',status='old')
         read(10,nml = inputwrite,err=2)
 2       continue
         write(6,nml = inputwrite)
         call flush(6)
         close(10)
         if(iwrt_reac_distri.eq.1)then
            write(6,*)' iwrt_reac_distri=1 --> writing distriREAC files'
         elseif(iwrt_reac_distri.eq.2)then
            write(6,*)' iwrt_reac_distri=2 --> writing Cvj files'
         endif
         
!     radial grid integration steps (in angstroms)
         
      if(npun2.gt.1)then
          div=dble(npun2-1)
          ah2=(rfin2-rmis2)/div
      else
          ah2=0.d0
      endif

      if(npun1.gt.1)then
        div=dble(npun1-1)
        ah1 = (rfin1-rmis1)/div
      else
        ah1=0.d0
      endif
      if(npun1.gt.1.and.npun2.gt.1)then
           steptot=ah1*ah2
      elseif(npun1.gt.1.and.npun2.le.1)then
           steptot=ah1
      elseif(npun1.le.1.and.npun2.gt.1)then
           steptot=ah2
      else
           steptot=1.d0
      endif

!     angular grid: Gauss legendre

      allocate(wreal(nangu2),weight(nangu2),cgamma(nangu2)
     &     ,cgamprod(nangu2),weiprod(nangu2)
     &       , stat=ierror)

      norealproc_mem=norealproc_mem
     &     + 5*nangu2
      write(6,*)'norealproc_mem= ',norealproc_mem
      call flush(6)

      if(nangu.eq.1)then
         wreal(1)=1.d0
         weight(1)=1.d0
         cgamma(1)=1.d0
      else
         if(inc.eq.1)then
            if(j0.ne.0)then
               if(idproc.eq.0)write(6,*)' inc=1 --> j0=0, here j0=',j0
               call flush(6)
               call MPI_BARRIER(MPI_COMM_WORLD, ierror)
              stop
            endif
         elseif(inc.eq.2)then
            if(j0.gt.1)then
            if(idproc.eq.0)write(6,*)' inc=2 --> j0=0 or 1, here j0=',j0
            call flush(6)
            call MPI_BARRIER(MPI_COMM_WORLD, ierror)
            stop
            endif
         else    
            if(idproc.eq.0)write(6,*)' inc=1 or 2, and it is inc= ',inc
            call flush(6)
            call MPI_BARRIER(MPI_COMM_WORLD, ierror)
            stop
         endif
         call gauleg(weight,cgamma,nangu2)
         
         do iang=1,nangu2
            cgamprod(iang)=cgamma(iang)
            weiprod(iang)=weight(iang)
         enddo
         do iang=1,nangu
            cgamma(iang)=cgamprod(iang)
            wreal(iang)=weiprod(iang)*dble(inc)
         enddo
         
         if(idproc.eq.0)then
            write(6,"(//,10x
     &         ,'Angular grid to plot pes and wvp every ',i3
     &         ,//)")nangplot
            if(nangu2.eq.2*nangu)then
               write(6,*)'  >>>>>is nangu even????, nangu= ',nangu
               write(6,*)'   otherwise there could be some problems...'
            endif
            do iang=1,nangu,nangplot
               write(6,*)iang,dacos(cgamma(iang))*180.d0/pi
            enddo
         endif
      endif
      call flush(6)
      
!     basis set conditions  
      
      if(iabs(iparity).ne.1)then
          write(6,*)'  ipar= ',iparity
          write(6,*)' | iparity| must be 1 '
          call flush(6)
          call MPI_BARRIER(MPI_COMM_WORLD, ierror)
          stop
      endif
!      if(iparity*(-1)**(Jtot).lt.0)then
!          if(iommin.eq.0)then
!             write(6,*)' iommin=0 not for this parity'
!             call flush(6)
!             call MPI_BARRIER(MPI_COMM_WORLD, ierror)
!             stop
!          endif
!      endif
      if(jmax.gt.(nangu2-1)*inc)then
        write(6,*)'  jmax= ',jmax
     &        ,' > (nangu2-1)*inc= ',(nangu2-1)*inc
        call flush(6)
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        stop
      endif  
     
      return
      end subroutine input_grid

!--------------------------------------------------
!--------------------------------------------------
      
      end module mod_gridYpara_01y2
