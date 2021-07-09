      program boundlanz
*********************************************************************
**                      program bondlanz                           **
**     version 2    ,  February   2021                             **
**                 parallelized for MPI                            **
**  references:                                                    **
**              J. Chem. Phys.  107 (1997) 10085                   **
**              J. Chem. Phys.  109 (1998) 9391                    **
**              J. Chem. Phys.  123 (2005) 194309                  **
**              J. Chem. Phys.  125 (2006) 054102                  **
**              J. Phys. Chem. A 113 (2009) 14488                  ** 
**                                                                 **
**  Calculates A+BC bound states in grids                          **
**               b) in reactant Jacobi coordinates (iprod=2)       **
**                  with 
**                                                                 **
**       including several coupled diabatic electronic states      **
**                                                                 **
**       Using  Jacobi coordinates                                 **
**                                                                 **
**                                        1                        **  
**                                        ^                        **
**                      R=R_2             |                        **
**                  2<--------------------| r=R_1                  **
**                                        |                        **
**                                        |                        **
**                                        0                        **
**                                                                 **
**                                                                 **
**       For a general value of Jtot  in a body-fixed frame        **
**      with the z-axis along the R=R_2 reactant Jacobi vector.    **
**                                                                 **
**   Input
**  ------
**        The progagation is done in a grid of gamma, r and R      **
**                      and a basis for Omega                      **
**             The angular kinetic term is evaluated               **
**            using a DVR  method for each Omega value             ** 
**                                                                 **
**                           units: zots                           **
**       this program uses the files:                              **
**             parameter.h  :  definition of dimensions            **
**                                  and input data                 **
**       the potential matrix (for several states) is provided     **
**            externally as in the input of the example            **
*********************************************************************
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      use mod_lanczos_01y2
      implicit none
      integer :: ierr,nntotproc,nnbastotproc,i,iv,n1,n2,ival
      double precision :: xn,xnorm,ener,enertot,ediag,uno

      include "mpif.h"

c! Initialize MPI environment and get proc's ID and number of proc in
c! the partition.

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      write(name,'("sal.id",i3.3)')idproc
      open(6,file=name,status='unknown')
      write(6,'(40("="),//)')
      write(6,'(10x,"BoundLanz3 version 6 ",//)')
      write(6,*)' output of proc. idproc= ',idproc,' of nproc= ',nproc
      write(6,'(/,40("="),//)')

!     initialization of data

      call input_grid
      call ini_absorcion
      call  paralelizacion
      call pot0

!     determining basis

      call basis

!     reactants and products functions calculation

      call angular_functions
      
      call radial_functions01_read

      if(iprod.ne.2)then
         write(6,*)' iprod = ',iprod
         write(6,*)' it must be iprod=2 to assign bnd states'
         call flush(6)
         stop
      endif
      
      call product_radialf_read
      
!     reading potential

      call pot2
      
!    dimensioning vectors for propagation,
!     calculating kinetic energy terms, etc
!    to evaluate H phi in a radial/angular grid & basis for Omega/electronic

      call set_vectors_Hphi

!     Lanczos dimensions initialization

      call set_lanczos
      
**>>  initial guess for Lanczos iteration

      nntotproc=ntotproc(idproc)
      nnbastotproc=nbastotproc
      if(idproc.eq.0)then
         write(6,*)'  Lanczos initialization'
         call flush(6)
      endif

      call initia_lanzvec(rflanzproc,nntotproc,nnbastotproc)

* initial guess for Lanczos iterations

      do i=1,ntotproc(idproc)
         rpaqproc(i)=rflanzproc(i,2)
      enddo

* fftw3 initialization

         write(6,*)'  FFW3 initialization'
         call flush(6)
         call Hmatrix(rpaqproc,rHpaqrec,nntotproc,nnbastotproc)

**>> checking norm and energy of initial guess

      xn=0.d0
      do i=1,ntotproc(idproc)
         xn=xn+rpaqproc(i)*rpaqproc(i)
      enddo
      call MPI_REDUCE(xn,xnorm,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xnorm,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      if(idproc.eq.0)then
          write(6,*)idproc,' norm= ',xnorm
       endif
       
      ener=0.d0
      do i=1,ntotproc(idproc)
         ener=ener+rpaqproc(i)*rHpaqrec(i)
      enddo
      call MPI_REDUCE(ener,enertot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
      if(idproc.eq.0)then
         write(6,*)idproc, '  initial energy= '
     &                   ,enertot/conve1,enertot
         write(6,*)
         write(6,*)' entering lanzvalpar '
         call flush(6)
      endif
****>> DIAGONALIZATION using lanczos

      ival=0

      open(3,file='cont.data',status='old',err= 1234)
         ival=1
         read(3,*)n1,n2
         read(3,*)ntrue
         do i=1,ntrue
            read(3,*)n1,etrue(i)
         enddo
      close(3)
 1234 continue

      nkryini=1
      uno=1.d0
      if(ival.eq.0)then
         call lanzvalpariter(rflanzproc,rHpaqrec,rpaqproc
     &                    ,nntotproc,nbastotproc
     &                    ,etrue,ntrue,
     &                    nloop,ntimes,nkryini,criter,uno,
     &                    idproc,nproc)
      endif
      if(idproc.eq.0)then
         write(6,*)
         write(6,*)' --> Calculated eigenvalues with nkryini ',nkryini
         write(6,*)
         write(6,*)'    i          zots            cm-1             eV'
         write(6,*)' -------------------------------------------------'
         write(6,*)

         do i=1,ntrue
            if(etrue(i).lt.vcutmax)then
               write(6,'(1x,i5,3(2x,e15.7))')i,etrue(i)
     &                 ,etrue(i)/conve1
     &                 ,etrue(i)/conve1/8065.5d0
            endif
         enddo
         call flush(6)
      endif

      if(plotbnd.ne.1)then
***> calculating bound states
         write(6,*)' ------------------------------------------- '
         write(6,*)' --> starting calculation of eigenstates <-- '
         write(6,*)'        ',nvecmax,'  of ntrue= ',ntrue
         write(6,*)' ------------------------------------------- '
         write(6,*)
         call flush(6)
         do iv=nvecmin,min0(nvecmax,ntrue)

            ediag=etrue(iv)

            if(idproc.eq.0)then
               write(6,*)
               write(6,*)' ->  E(',iv,')= ',etrue(iv)/conve1/8065.5
               write(6,*)
               call flush(6)
            endif 

            call cgradvecpar(ediag,rflanzproc,rHpaqrec,rpaqproc
     &           ,nntotproc,nbastotproc
     &           ,nomxiter,criter,uno,idproc,nproc)

* output

            call bndassign(iv)
            call write_bnd(iv)
              
         enddo  ! iv
      else
*** > printing bound states to plot
         do iv=nvecmin,min0(nvecmax,ntrue)
            ediag=etrue(iv)

            if(idproc.eq.0)then
               write(6,*)
           write(6,*)' -> printing E(',iv,')= ',etrue(iv)/conve1/8065.5
               write(6,*)
               call flush(6)
            endif

            call read_bnd(iv)
            call plot_bnd(iv)

         enddo
      endif
*****************************
      call MPI_FINALIZE(ierr)
      
      stop
      end
      
*********************** initialanzvec  **********************************

      subroutine initia_lanzvec(rflanzproc,nntot,nnbastotproc)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      use mod_lanczos_01y2, only:Emax_lanzini
      implicit none

      integer :: nntot,nnbastotproc,i
      double precision :: xx,xx_tot
      integer :: itotp,icanp,ielec,iom,iangp,ir,ir1,ir2,ierr
      double precision :: rflanzproc(nnbastotproc,0:2)

      include "mpif.h"

*     ************************************
*     * Initialization of lanczos vector *
*     * for eigenvalues and eigenvectors *
*     ************************************

      do itotp=1,nnbastotproc
         rflanzproc(itotp,0)=0.d0
         rflanzproc(itotp,1)=0.d0
         rflanzproc(itotp,2)=0.d0
      enddo

      do itotp=1,ntotproc(idproc)

         call indiproc(itotp,icanp,ielec,iom,iangp,ir,ir1,ir2)    
         if(vvv(ir,iangp,ielec,ielec).lt.Emax_lanzini)then
            rflanzproc(itotp,2)=1.d0
         endif
      enddo

* renormalizing

      xx=0.d0
      do i=1,ntotproc(idproc)
         xx=xx+rflanzproc(i,2)*rflanzproc(i,2)
      enddo

      call MPI_REDUCE(xx,xx_tot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(xx_tot,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      if(dabs(xx_tot).lt.1.d-10)then
         write(6,*)' --> Norm too low initia_lanzvec <--'
         write(6,*)'  Emax_lanzini= ',Emax_lanzini
         call flush(6)
        stop
      endif

       do i=1,ntotproc(idproc)
          rflanzproc(i,2)=rflanzproc(i,2)/dsqrt(xx_tot)
       enddo

      return
      end


*************************************************************

      subroutine Hmatrix(eigen,Hu,nntotproc,nnbastotproc)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      implicit none
      integer :: nntotproc,nnbastotproc
      double precision :: Hu(nnbastotproc),eigen(nnbastotproc)
      integer :: i

      include "mpif.h"

      do i=1,nnbastotproc
         rpaqproc(i)=eigen(i)
         rHpaqrec(i)=0.d0
      enddo

      call difs

      do i=1,nnbastotproc
         Hu(i)=rHpaqrec(i)
      enddo
     
      return
      end

******************************************************************
      subroutine write_bnd(iv)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      use mod_lanczos_01y2
      implicit real*8(a-h,o-z)
      include "mpif.h"

      dimension proj(nangu,nvini:nvmax,ncanmax)
      dimension fun(nangu,nvini:nvmax), funtot(nangu,nvini:nvmax)
*
*  writting bound state into files for each proccessor 
*

       write(name,'("bnd.iv",i3.3,".ip",i3.3)')iv,idproc
       open(10,file=name,status='unknown')
       write(name,'("grid.ip",i3.3)')idproc
       open(11,file=name,status='unknown')

         write(10,*)iv,etrue(iv)/conve1/8065.5d0

         write(11,*)ntotproc(idproc)
         write(11,*)rmis1,rfin1,npun1
         write(11,*)rmis2,rfin2,npun2
         write(11,*)nangu,inc
         do i=1,ntotproc(idproc)
            call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
            iang=indangreal(iangp,idproc)

         write(10,*)rpaqproc(i)
         write(11,'(9(1x,i7))')i,icanp,ielec,iom,iangp,ir,ir1,ir2,iang

          enddo

       close(10)
       close(11)
     
      return
      end subroutine write_bnd
******************************************************************
      subroutine read_bnd(iv)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      use mod_lanczos_01y2
      implicit real*8(a-h,o-z)
      include "mpif.h"

      integer :: i,j,k,l
      double precision :: x1,x2,x3
*
*  writting bound state into files for each proccessor 
*

       write(name,'("bnd.iv",i3.3,".ip",i3.3)')iv,idproc
       open(10,file=name,status='old')

         read(10,*)

         do i=1,ntotproc(idproc)
            call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
            iang=indangreal(iangp,idproc)

            read(10,*)rpaqproc(i)
          enddo

       close(10)
     
      return
      end subroutine read_bnd


******************************************************************
      subroutine bndassign(iv)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      use mod_lanczos_01y2
      implicit real*8(a-h,o-z)
      include "mpif.h"
*
*  assigning bound states 
* 
      dimension disv1(nvini:nvmax),funv1(npun2,nangu)
      dimension disv2(nviniprod:nvmaxprod),funv2(npun1,nangu)
      dimension disomg(iommin:iommax)

      iv1max=-1
      iv2max=-1
      costet=0.d0
      j=0
*
*  v1 distribution in r1 
*

      write(6,*)'   dis in v1 '
      pop=0.d0
      do iv1=nvini,nvmax
         disv1(iv1)=0.d0
         yyy=0.d0
         yyytot=0.d0
         do jcanproc=1,ncanproc(idproc)
            do iang=1,nangu
            do ir2=1,npun2
               funv1(ir2,iang)=0.d0
            enddo
            enddo
          
            do i=1,ntotproc(idproc)
               call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
               iang=indangreal(iangp,idproc)
               if(icanp.eq.jcanproc)then
                  xxx=fd(ir1,iv1,j,ielec)*rpaqproc(i)
                  funv1(ir2,iang)=funv1(ir2,iang)+xxx
               endif
            enddo

            xxx=0.d0
            do iang=1,nangu
            do ir2=1,npun2
               xxx=xxx+funv1(ir2,iang)*funv1(ir2,iang)
            enddo
            enddo

            yyy=yyy+xxx
         enddo  ! jcanproc

         nnn=1
         call MPI_REDUCE(yyy,yyytot,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
         disv1(iv1)=yyytot
         write(6,*)iv1,disv1(iv1)
         pop=pop+disv1(iv1)
      enddo  ! iv1
      write(6,*)' end   dis in v1, totpop= ',pop
     
*
*  v2 distribution in r2 
*

      write(6,*)
      write(6,*)'   dis in v2 '
      pop=0.d0
      do iv2=nviniprod,nvmaxprod
         disv2(iv2)=0.d0
         yyy=0.d0
         yyytot=0.d0
         do jcanproc=1,ncanproc(idproc)
            do iang=1,nangu
            do ir1=1,npun1
               funv2(ir1,iang)=0.d0
            enddo
            enddo
          
            do i=1,ntotproc(idproc)
               call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
               iang=indangreal(iangp,idproc)
               if(icanp.eq.jcanproc)then
                  xxx=fdprod(ir2,iv2,j,ielec)*rpaqproc(i)
                  funv2(ir1,iang)=funv2(ir1,iang)+xxx
               endif
            enddo

            xxx=0.d0
            do iang=1,nangu
            do ir1=1,npun1
               xxx=xxx+funv2(ir1,iang)*funv2(ir1,iang)
            enddo
            enddo

            yyy=yyy+xxx
         enddo  ! jcanproc

         nnn=1
         call MPI_REDUCE(yyy,yyytot,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
         disv2(iv2)=yyytot
         write(6,*)iv2,disv2(iv2)
         pop=pop+disv2(iv2)
      enddo  ! iv2
      write(6,*)' end   dis in v2, totpop= ',pop

*
*   angular mean value
*

      xmean=0.d0
      do i=1,ntotproc(idproc)
         call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
         iang=indangreal(iangp,idproc)
         xxx=rpaqproc(i)*rpaqproc(i)
         xcos=cgamma(iang)
         xmean=xmean+xxx*xcos
      enddo
      nnn=1
      call MPI_REDUCE(xmean,xmeantot,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)

      write(6,*)
      write(6,*)'   <cos(theta)> = ',xmeantot

*
*  Omega distribution
*
      write(6,*)'   dis in Omega'
      pop=0.d0
      do jom=iommin,iommax
         disomg(jom)=0
            yyy=0.d0
            do i=1,ntotproc(idproc)
               call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
               iang=indangreal(iangp,idproc)
               if(iom.eq.jom)then
                  xxx=rpaqproc(i)
                  yyy=yyy+xxx*xxx
               endif
            enddo

            nnn=1
         call MPI_REDUCE(yyy,yyytot,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
         disomg(jom)=yyytot
         write(6,*)jom,disomg(jom)
         pop=pop+disomg(jom)
      enddo ! iom
      write(6,*)' end   dis in Omg, totpop= ',pop
*
*  electronic distribution
*
      write(6,*)'   dis in e-state'
      pop=0.d0
      do jelec=1,nelec
            yyy=0.d0
            do i=1,ntotproc(idproc)
               call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
               iang=indangreal(iangp,idproc)
               if(ielec.eq.jelec)then
                  xxx=rpaqproc(i)
                  yyy=yyy+xxx*xxx
               endif
            enddo

            nnn=1
         call MPI_REDUCE(yyy,yyytot,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
         write(6,*)jelec,yyytot
         pop=pop+yyytot
      enddo ! jelec
      write(6,*)' end   dis in e-state, totpop= ',pop

      return
      end subroutine bndassign

******************************************************************
      subroutine plot_bnd(iv)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      use mod_lanczos_01y2
      implicit real*8(a-h,o-z)
      include "mpif.h"

      dimension proj(nangu,ncanmax)
      dimension fun(nangu), funtot(nangu)

      dimension projr(npun1,ncanmax)
      dimension funr(npun1), funrtot(npun1)
      dimension funr1r2bond(npun1)
      dimension funr1r2bond_tot(npun1)
      dimension dis_opt(npun1),iang_opt(npun1),ir2_opt(npun1)
*        
*   plots bound state averaging over r_small
*
      if(idproc.eq.0)then   
         write(name,'("bndplot.i",i3.3)')iv
         open(10,file=name,status='unknown')
      endif

      do jr2=1,npun2
         r2=rmis2+dble(jr2-1)*ah2
         proj(:,:)=0.d0

         do i=1,ntotproc(idproc)
            call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
            iang=indangreal(iangp,idproc)
            if(ir2.eq.jr2)then
               proj(iang,icanp)=proj(iang,icanp)+rpaqproc(i)*rpaqproc(i)
            endif
         enddo

         do iang=1,nangu
            fun(iang)=0.d0
            do ican=1,ncanproc(idproc)
               fun(iang)=fun(iang)+proj(iang,ican)
            enddo
         enddo
         
         nnn=nangu
         call MPI_REDUCE(fun,funtot,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)

         if(idproc.eq.0)then
         do iang=1,nangu
            ang=dacos(cgamma(iang))
            write(10,'(100(1x,e15.7))')r2,ang,funtot(iang)            
         enddo
         write(10,'()')
         endif
      enddo ! ir2

      if(idproc.eq.0) close(10)

*        
*   plots bound state averaging over gamma
*
      if(idproc.eq.0)then   
         write(name,'("bndr1r2.i",i3.3)')iv
         open(10,file=name,status='unknown')
      endif

      do jr2=1,npun2
         r2=rmis2+dble(jr2-1)*ah2
         projr(:,:)=0.d0

         do i=1,ntotproc(idproc)
            call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
            iang=indangreal(iangp,idproc)
            if(ir2.eq.jr2)then
               projr(ir1,icanp)=projr(ir1,icanp)+rpaqproc(i)*rpaqproc(i)
            endif
         enddo

         do ir1=1,npun1
            funr(ir1)=0.d0
            do ican=1,ncanproc(idproc)
               funr(ir1)=funr(ir1)+projr(ir1,ican)
            enddo
         enddo
         
         nnn=npun1
         call MPI_REDUCE(funr,funrtot,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)

         if(idproc.eq.0)then
         do ir1=1,npun1
            r1=rmis1+dble(ir1-1)*ah1
            write(10,'(100(1x,e15.7))')r1,r2,funrtot(ir1)            
         enddo
         write(10,'()')
         endif
      enddo ! ir2

      if(idproc.eq.0) close(10)

*        
*   plots bound state averaging over gamma
*
      if(idproc.eq.0)then   
         write(name,'("bndbondx1x2.i",i3.3)')iv
         open(10,file=name,status='unknown')
      endif

      xbondang=dcos(135.5d0*pi/180.d0)

      do jr1=1,npun1
         x1=rmis1+dble(jr1-1)*ah1
         ir1=jr1
         r1=x1

         dis_opt(:)=1.d10
         iang_opt(:)=0
         ir2_opt(:)=0
         do jr2=1,npun1
            x2=rmis1+dble(jr2-1)*ah1

            xm01=xm0+xm1
            gam=xm1/xm01

            Rg=x2*x2+x1*x1*gam*gam-2.d0*xbondang*gam*x1*x2
            Rg=dsqrt(Rg)
            ir2_opt(jr2)=(Rg-rmis2)/ah2
            if(ir2_opt(jr2).lt.1)then
                ir2_opt(jr2)=1
            elseif(ir2_opt(jr2).gt.npun2)then
               ir2_opt(jr2)=npun2
            endif
            costet=(x2*xbondang-gam*x1)/Rg
            if(dabs(costet).gt.1.d0)costet=costet/dabs(costet)
            dis=1.d10
            do iang=1,nangu
               xxx=dabs(cgamma(iang)-costet)
               xxx=xxx*xxx
               if(dabs(xxx).lt.dis)then
                  dis=xxx
                  iang_opt(jr2)=iang
               endif     
            enddo
           
!            r2=rmis2+dble(ir2_opt(jr2)-1)*ah2
!            r2bond=(xm1*xm1/((xm0+xm1)*(xm0+xm1)))*r1*r1+r2*r2
!     &           +2.d0*(xm1/(xm0+xm1))*r1*r2*cgamma(iang_opt(jr2))
!            r2bond=dsqrt(r2bond)
!            ccc=(xm1/(xm0+xm1))*x1+R2*cgamma(iang_opt(jr2))
!            ccc=ccc/r2bond        
!            write(69,*)x1,x2,xbondang
!     &            ,r2bond,ccc
!     &           ,rmis2+dble(ir2_opt(jr2)-1)*ah2
!     &           ,cgamma(iang_opt(jr2))
!     &         ,dis_opt(jr2),iang_opt(jr2),ir2_opt(jr2)
         enddo                  ! jr2

         funr1r2bond(:)=0.d0
         funr1r2bond_tot(:)=0.d0
         do i=1,ntotproc(idproc)
            call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
            iang=indangreal(iangp,idproc)
            if(jr1.eq.ir1)then
            do jr2=1,npun1
            if(ir2.eq.ir2_opt(jr2).and.iang.eq.iang_opt(jr2))then
               funr1r2bond(jr2)=rpaqproc(i)
            endif
            enddo
            endif
         enddo
         
         nnn=npun1
         call MPI_REDUCE(funr1r2bond,funr1r2bond_tot,nnn
     &                     ,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 
         if(idproc.eq.0)then
         do ir2=1,npun1
            x2=rmis1+dble(ir2-1)*ah1
            write(10,'(100(1x,e15.7))')x1,x2,funr1r2bond_tot(ir2)
         enddo
         write(10,'()')
         endif
          
      enddo ! jr1

     
      if(idproc.eq.0) close(10)
      
      return
      end subroutine plot_bnd
