      module mod_PsiE_01y2
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      save
*     public !

      integer :: nofe,ifile_PsiE
      character*50 :: file_resonancesE
      double precision,allocatable :: resonancesE_ev(:),efe(:)
      complex*16, allocatable :: zfE(:,:)
      contains

**********************
*     input_PsiE     *
**********************
!--------------------------------------------------
       subroutine input_PsiE
!--------------------------------------------------
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2

      implicit none
      integer :: imem,iE,ierror
      logical :: exists
*********************************************************
      namelist /inputPsiE/nofe,file_resonancesE
*********************************************************

      nofe=0
      if( iphoto > 0 )then
         write(6,'(40("="),/,10x,"mod_PsiE",/,40("="))')
         write(6,*)
         write(6,*)'  grid and basis data'
         write(6,*)'  -------------------'
         open(10,file='input.dat',status='old')
         read(10,nml = inputPsiE,err=999) 
         write(6,nml = inputPsiE)
         call flush(6)
         close(10)
      
         if( nofe > 0 )then

            allocate(resonancesE_eV(nofe),EfE(nofe)
     &         ,zfe(nbastotproc,nofe)
     &        ,stat=ierror)
         
            if(ierror.ne.0)then
               write(*,*)" error in input_PsiE for zfe "
               call flush(6)
               stop
            else
               ifile_PsiE=10
               imem= nbastotproc*2
              norealproc_mem=norealproc_mem+imem
              write(6,*)' zfe in input_PsiE, proc= '
     &                    ,idproc
     &          ,' imem=',imem
     &          ,' memory(Gb)= '
     &           ,dble(imem*8)*1.d-9
               call flush(6)
            end if

            write(6,*)' --- Reading resonance energies in file= '
     &           ,file_resonancesE
            call flush(6)
            open(ifile_PsiE,file=file_resonancesE,status='old')
            do ie=1,nofe
               read(ifile_PsiE,*) resonancesE_eV(ie)
               write(6,*)ie,resonancesE_eV(ie)
               EfE(ie) = resonancesE_eV(ie)*ev2cm*conve1
            end do
            close(ifile_PsiE)

            zfe(:,:) = dcmplx(0.d0,0.d0)
         
         end if  ! nnofe > 0

      end if                    ! iphoto > 0
      
      return
 999  write(6,*)'  no PsiE information, skipping !!'
      call flush(6)
      close(10)
      end
                 
***************************funE  ******************************

      subroutine funE(it)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      implicit none
      integer :: ikcheb,it,ie,i
      double precision :: d,expo,deno,E,Es
      complex*16 :: zfactor,zexpo,zCkcheby
      include "mpif.h"

      ikcheb=it
      d=2.d0
      if(ikcheb.eq.0)d=1.d0
      zfactor=dcmplx(d*hbr/delta2,0.d0)

      do ie=1,nofe

         E=efe(ie)
         Es=(E-emindlt)/delta2
         expo=-dble(ikcheb)*dacos(Es)
         zexpo=cdexp(dcmplx(0.d0,expo))
         deno=1.d0-Es*Es
         deno=dsqrt(deno)
         zCkcheby=zfactor*zexpo/dcmplx(deno,0.d0)
         do i=1,ntotproc(idproc)
            zfe(i,ie)=zfe(i,ie)+ zCkcheby*dcmplx(rpaqproc(i),0.d0)
         enddo
      enddo ! ie=1,nofe

      return
      end
**********************************************************************

      subroutine plot_zfe(iloop)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      implicit none
      double precision :: fun(npun1,nangu), funtot(npun1,nangu)
      double precision :: funrp(npun2,nangu), funtotrp(npun2,nangu)
      double precision :: funr1r2(npun1,npun2,nangu)
     &                   , funr1r2tot(npun1,npun2,nangu)
      integer :: ifE,jr2,jelec,iang,nnn,iloop,ierr
      integer :: i,icanp,ielec,iom,iangp,ir,ir1,ir2,jom
      double precision :: r2,r1,sum,cang
      
      include "mpif.h"


* writting the wvp at each iloop
* in the same angular grid as the pes
* (i.e. each nangplot points of the angular grid)

      do ifE=1,nofe

      do jelec=1,nelec
         if(idproc.eq.0)then
            write(name,'("fe",i2.2,".elec",i2.2,"r1ang")')
     &                            ife,jelec
            open(10,file=name,status='unknown')
         endif

         fun(:,:)=0.d0
         funtot(:,:)=0.d0

         do i=1,ntotproc(idproc)
            call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
            iang=indangreal(iangp,idproc)
            if(ielec.eq.jelec)then
                fun(ir1,iang)=fun(ir1,iang)
     &               +dreal(zfe(i,ife)*dconjg(zfe(i,ife)))
            endif
         enddo

         nnn=npun1*nangu
         call MPI_REDUCE(fun,funtot,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)

         if(idproc.eq.0)then
            do ir1=1,npun1,n1plot
               r1=rmis1+dble(ir1-1)*ah1
               sum=0.d0
               do iang=1,nangu
                  if(dabs(funtot(ir1,iang)).lt.1.d-90)then
                     funtot(ir1,iang)=0.d0
                  endif
                  sum=sum+funtot(ir1,iang)
                  cang=cgamma(iang)
                  write(10,'(500(1x,e15.7))')r1,cang
     &                    ,(funtot(ir1,iang))
               enddo
               write(10,'()')
            enddo


         endif ! idproc = 0

         if(idproc.eq.0)then
            close(10)
         endif

* plotting wf_E for fixed r_peq

         if(idproc.eq.0)then
            write(name,'("fe",i2.2,".elec",i2.2,".r2gam")')
     &                      ife,jelec
            open(10,file=name,status='unknown')
         endif

         funrp(:,:)=0.d0
         funtotrp(:,:)=0.d0

         do i=1,ntotproc(idproc)
            call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
            iang=indangreal(iangp,idproc)
            if(ielec.eq.jelec)then
                funrp(ir2,iang)=funrp(ir2,iang)
     &               +dreal(zfe(i,ife)*dconjg(zfe(i,ife)))
            endif
         enddo

         nnn=npun2*nangu
         call MPI_REDUCE(funrp,funtotrp,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)

         if(idproc.eq.0)then
            do iang=1,nangu,nangplot
               do ir2=1,npun2,n2plot
                  r2=rmis2+dble(ir2-1)*ah2
                  if(dabs(funtotrp(ir2,iang)).lt.1.d-90)then
                          funtotrp(ir2,iang)=0.d0
                  endif
                       write(10,'(500(1x,e15.7))')r2,cgamma(iang)
     &                 ,funtotrp(ir2,iang)
               enddo
               write(10,'()')
            enddo

            close(10)
         endif


         do jom=iommin,iommax

            funr1r2(:,:,:)=0.d0
            funr1r2tot(:,:,:)=0.d0
            do i=1,ntotproc(idproc)
               call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
               iang=indangreal(iangp,idproc)
               if(ielec.eq.jelec.and.iom.eq.jom)then
                funr1r2(ir1,ir2,iang)=funr1r2(ir1,ir2,iang)
     &               +dreal(zfe(i,ife)*dconjg(zfe(i,ife)))
            endif
         enddo

         nnn=npun1*npun2*nangu
         call MPI_REDUCE(funr1r2,funr1r2tot,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
            
            
            if(idproc.eq.0)then
               write(name,'("fer1r2",i2.2,".elec",i2.2,".iom",i2.2)')
     &                            ife,jelec,jom
               open(10,file=name,status='unknown')
               do ir1=1,npun1,n1plot
                  r1=rmis1+dble(ir1-1)*ah1
                  do ir2=1,npun2,n2plot
                     r2=rmis2+dble(ir2-1)*ah2
                     write(10,'(2000(1x,e15.7))')r1,r2
     &                 ,(funr1r2tot(ir1,ir2,iang),iang=1,nangu)
                  enddo
                  write(10,'()')
               enddo
               close(10)
            endif

         enddo

      enddo ! ielec
      enddo ! ife

      return
      end
!--------------------------------------------------
!--------------------------------------------------

      end module mod_PsiE_01y2
