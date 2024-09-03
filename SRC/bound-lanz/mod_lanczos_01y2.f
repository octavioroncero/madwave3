      module mod_lanczos_01y2
      
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      implicit none

      save
      double precision, allocatable :: rflanzproc(:,:),etrue(:)
      double precision ::  Emax_lanzini_eV,Emax_lanzini
     &                    ,criter_eV,criter

      integer :: nkryini,nloop,ntimes,nkrymax,nvecmin,nvecmax,nomxiter
      integer :: ntrue,plotbnd

      contains
!---------------------------------
      subroutine set_lanczos
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      implicit none
      integer :: imem

      namelist /lanczos_parameters/ Emax_lanzini_eV
     &                             ,ntimes
     &                             ,nkryini,nkrymax
     &                             ,nvecmin,nvecmax
     &                             ,criter_eV,nomxiter
     &                             ,plotbnd

      ntimes=1000
      nkrymax=10000
      nkryini=9000
      criter_eV=1.d-4
      nomxiter=nkrymax
      nvecmin=1
      nvecmax=10
      plotbnd=0
      write(6,'(40("="),/,10x,"Lanczos_01y2_mod",/,40("="))')
      write(6,*)
      write(6,*)'  lanczos_parameters'
      write(6,*)'  ------------------'
      call flush(6)
         open(10,file='input.dat',status='old')
         read(10,nml = lanczos_parameters)
         write(6,nml = lanczos_parameters)
         call flush(6)
         close(10)

      nloop=nkrymax/ntimes
      Emax_lanzini=Emax_lanzini_eV*ev2cm*conve1
      criter=criter_eV*ev2cm*conve1

      imem=3*nbastotproc
      write(6,*)' flanzproc in set_lanczos '
     &    ,'  in proc= ',idproc,' is= '
     &     ,dble(imem)*1.d-9,' Gb'
      norealproc_mem=norealproc_mem+imem+nkrymax
      
      allocate(rflanzproc(nbastotproc,0:2),etrue(nkrymax))

       write(6,*)' ending set_lanczos, proc= ',idproc
     &     ,' norealproc_mem=',norealproc_mem
     &     ,' nointegerproc_mem=',nointegerproc_mem
     &     ,' memory(Gb)= '
     &   ,dble(nointegerproc_mem*4+norealproc_mem*8)*1.d-9
      call flush(6)

     
      
      call flush(6)
      
         
      end subroutine set_lanczos
!---------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end  module mod_lanczos_01y2
