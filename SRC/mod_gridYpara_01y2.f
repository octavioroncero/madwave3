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

*     products states analysis

      integer :: iprod,nviniprod,nvmaxprod,jiniprod,jmaxprod
      integer :: iomminprod,iommaxprod,n2prod1,n2prod0
      integer :: nangproj0,nangproj1
      real*8 :: rbalinprod
      
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
     &                       ,nangu,nangplot
     &     ,Jtot,iparity,inc,nelecmax,iommin,iommax,j0
     &     ,jini,jmax,nvini,nvmax
     &                       ,nvref,jref,iomref,ielecref
*********************************************************
      namelist /inputprod/iprod
     &                   ,nviniprod,nvmaxprod
     &                   ,jiniprod,jmaxprod
     &     ,iomminprod,iommaxprod
     &     ,Rbalinprod,n2prod0,n2prod1,nangproj0,nangproj1
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
         iommaxprod=Jtot
         n2prod1=npun2
         n2prod0=1
         nangproj0=1
         nangproj1=nangu
         iprod=0
      write(6,*)
      write(6,*)'  products data'
      write(6,*)'  -------------'
         open(10,file='input.dat',status='old')
         read(10,nml = inputprod)
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
         
         
!     radial grid integration steps (in angstroms)
         
      ah2 = (rfin2-rmis2)/dble(npun2-1)
      if(npun1.gt.1)then
        div=dble(npun1-1)
        ah1 = (rfin1-rmis1)/div
        steptot=ah1*ah2
      else
        ah1=0.d0
        steptot=ah2
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
      if(iparity*(-1)**(Jtot).lt.0)then
          if(iommin.eq.0)then
             write(6,*)' iommin=0 not for this parity'
             call flush(6)
             call MPI_BARRIER(MPI_COMM_WORLD, ierror)
             stop
          endif
      endif
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
      subroutine paralelizacion
!--------------------------------------------------
      implicit none
      include "mpif.h"
      real*8 :: div
      integer :: ierror,lmax,ip,iparbc,ielec,iom,iomdi,iomat
      integer :: isignexp,isign,ipar,iomtot,ifail,jmin,isi
      integer :: i,iang,iangproc,ican,io,iomc,iomind,iomindc,ncan
      integer :: ipang,ipc,ipomg,nnproc,nomg,nprocang

      write(6,*)
      write(6,*)' Memory allocating among processors, id= ',idproc
      write(6,*)
      call flush(6)      
!--------------------------------
!  dimensions for parallelization
!--------------------------------
      nprocdim=nproc
      nomgdim=iommax-iommin+1
      nomg=nomgdim
      ncanmax=nomg*nelecmax
      write(6,*)'nprocdim= ',nprocdim,'  nomgdim= ',nomgdim
      if(nprocdim.le.nomgdim)then
        if(mod(nomgdim,nprocdim).ne.0)then
            write(6,*)' no. of Omegas no divisible by nproc'
            call MPI_BARRIER(MPI_COMM_WORLD, ierror)
            call flush(6)
            stop
        endif
        nomgproc=nomg/nprocdim
        nomgprocdim=nomgdim/nprocdim
        nprocang=1
        nanguproc=nangu
        nanguprocdim=nangu
      else
        if(mod(nprocdim,nomgdim).ne.0)then
            write(6,*)' no. of nproc no divisible by Omegas'
            call flush(6)
            call MPI_BARRIER(MPI_COMM_WORLD, ierror)
            stop
        endif
        if(mod(nangu*nomgdim,nprocdim).ne.0)then
            write(6,*)' no. of angles,Omegas no divisible by nproc'
            call flush(6)
            call MPI_BARRIER(MPI_COMM_WORLD, ierror)
            stop
        endif
        nomgproc=1
        nomgprocdim=1
        nprocang=nprocdim/nomg
        nanguproc=nangu/nprocang
        nanguprocdim=nangu/nprocang
      endif
      ncanprocdim=nelecmax*nomgprocdim

      if(nomgproc.gt.nomgprocdim)then
         write(6,*)'  nomgproc= ',nomgproc,' > nomgprocdim '
         call flush(6)
         call MPI_BARRIER(MPI_COMM_WORLD, ierror)
         stop
      endif

* index parallel and mpi initialization

        write(6,*)'  allocating indexes '
        allocate( ncanproc(0:nprocdim-1)
     &    ,ibasproc(ncanmax,0:nprocdim-1)
     &    ,indiomreal(nomgprocdim,0:nprocdim-1)
     &    ,indangreal(nanguprocdim,0:nprocdim-1)
     &    ,indproc(iommin:iommax,nangu)
     &    ,ntotproc(0:nprocdim-1)
     &    , ncouproc(0:nprocdim-1)
     &    ,ipcou(nprocdim,0:nprocdim-1)
     &       , stat=ierror)

        nointegerproc_mem=nointegerproc_mem
     &    +nprocdim*ncanmax*nprocdim
     &    +nomgprocdim*nprocdim
     &    +nanguprocdim*nprocdim
     &    +(iommax-iommin+1)*nangu
     &    +nprocdim*(nprocdim+2)
      write(6,*)'nointegerproc_mem= ',nointeger_mem
      call flush(6)

      if(ierror.ne.0)then
         write(*,*)" error in initmem for indexes for parallelitation "
         call flush(6)
         stop
      endif
      do ip=0,nprocdim-1
         ncanproc(ip)=0
         ntotproc(ip)=0
         ncouproc(ip)=0
         do ican=1,ncanmax
            ibasproc(ican,ip)=0
         enddo
         do io=1,nomgprocdim
            indiomreal(io,ip)=0
         enddo
         do iang=1,nanguprocdim
            indangreal(iang,ip)=0
         enddo
         do nnproc=1,nprocdim
            ipcou(nnproc,ip)=0
         enddo
      enddo
      do iang=1,nangu
      do iom=iommin,iommax
         indproc(iom,iang)=0
      enddo
      enddo

* assignment of omega's and iang's to different processors

      write(6,*)'  Omegas/Angles per procesor '
      call flush(6)

      ip=-1
      do ipomg=1,nomgdim/nomgproc
         do ipang=1,nangu/nanguproc
            ip=ip+1

            do iomind=1,nomgproc
               iom=iommin+(ipomg-1)*nomgproc+iomind-1
               indiomreal(iomind,ip)=iom
            enddo

            iangproc=0
            do iang=nangu-ipang+1,1,-nprocang
               iangproc=iangproc+1
               indangreal(iangproc,ip)=iang
            enddo
            write(6,*)' procesor ',ip
            write(6,'(10x,"Omega= ",50(1x,i3))')(indiomreal(iomind,ip)
     &                              ,iomind=1,nomgproc)

            write(6,"(10x,'Angles= ',50(1x,i3))")
     &           (indangreal(iangproc,ip),iangproc=1,nanguproc)
            call flush(6)

         enddo
      enddo

      ncan=0
      do ielec=1,nelecmax
      do iom=iommin,iommax
         ncan=ncan+1
         do ip=0,nproc-1
         do iomind=1,nomgproc
            if(iom.eq.indiomreal(iomind,ip))then
                ncanproc(ip)=ncanproc(ip)+1
                ibasproc(ncanproc(ip),ip)=ncan
                write(6,*)'    in processor= ',ip
                call flush(6)
            endif
         enddo
         enddo
      enddo
      enddo
**>> coupled processors

      write(6,*)'  coupled Processors '
      write(6,*)
      call flush(6)

      ncouprocmax=0
      do ip=0,nprocdim-1
         ncouproc(ip)=0
         do ipc=0,nprocdim-1
            ifail=0
            do iomind=1,nomgproc
            do iomindc=1,nomgproc
               iom=indiomreal(iomind,ip)
               iomc=indiomreal(iomindc,ipc)
               if(iabs(iom-iomc).le.1)ifail=1
            enddo
            enddo
            if(ifail.eq.1)then
               ncouproc(ip)=ncouproc(ip)+1
               ipcou(ncouproc(ip),ip)=ipc
            endif
         enddo
         if(ncouproc(ip).gt.ncouprocmax)ncouprocmax=ncouproc(ip)

         write(6,*)' Processor ',ip,' coupled to '
     &        ,(ipcou(ipc,ip),ipc=1,ncouproc(ip))
         call flush(6)

      enddo
      write(6,*)' Max. no. of coupled processors= ',ncouprocmax
      call flush(6)

      return
      end subroutine paralelizacion

!--------------------------------------------------
!--------------------------------------------------
      
      end module mod_gridYpara_01y2
