      program madwave3
*********************************************************************
**                      program MadWave3                           **
**     version 6    ,  February   2018                             **
**                 parallelized for MPI and openMP                 **
**                 last modification, January 2020                 **
**  references:                                                    **
**              J. Chem. Phys.  107 (1997) 10085                   **
**              J. Chem. Phys.  109 (1998) 9391                    **
**              J. Chem. Phys.  123 (2005) 194309                  **
**              J. Chem. Phys.  125 (2006) 054102                  **
**              J. Phys. Chem. A 113 (2009) 14488                  ** 
**                                                                 **
**  For dynamical studies of A+BC reaction dynamics                **
**       for:    a) collisions (iphoto=0)                          ** 
**               b)  photoinitiated (iphoto=1)                     **
**               c)  "bright" approximation (iphoto=2)             **
**       and:                                                      ** 
**   state-2-state probabilities (and amplitudes):                 **
**               a) in product jacobi coordinates (iprod=1)        **
**               b) in reactant Jacobi coordinates (iprod=2)       **
**               c) only evaluate total reaction prob (iprod=0)    **
**                                                                 **
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
**       Transforming either to reactant or product Jacobi         **
**          to solve state-to-state probabilitites                 **
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
**               Chebychev propagator in K iterations              **
**           using only real part of the wavepacket                **
**            equivalent to Gray, Balint-Kurti JCP (1998)          **
**                                                                 **
**       this program uses the files:                              **
**             input.dat  :  all input data distributed among      **
**                           different namelist:                   **
**                                  /inputgridbase/
**                                  /inputpotmass/
**                                  /inputprod/
**                                  /inputprocess/                 **
**                                  /inputbnd/ --> photo initiated **
**                                  /inputcol/ --> for collisions
**                                  /inputtime/                    **
**                                                                 **
**             cont.data    :  continuation file                   **
**                            it must contain    0 0 0 1           **
**                             at the beginning                    **
**             files in directories pot,dip,func                   **
**                   generated by the preparation progeam          **
**                           main_potini                           **
**              
**                                                                 **
**       the potential matrix (for several states) is provided     **
**            externally as in the input of the example            **
**             corresponding to H3 case.                           **
**
**   Output
**   ------
**        S2prod.XXXXX  provides energy (cm^{-1}) 
**                               total reaction probability
**                               total probability (equal to 1 if converged)
**
*********************************************************************

      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_photoini_01y2
      use mod_colini_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      use mod_flux_01y2
      use mod_coortrans01_02
      
      implicit none
      include "mpif.h"     
      character*40 filename
      integer ierr
      integer :: ntimes, nloop,kminelastic
      integer :: it,iit0,itinter,iloop,iloop0,indt,iloop00
      integer :: i,icanp,ielec,iom,iangp,ir,ir1,ir2
      real*8 :: fabs,xxx,rHnorm,autocor,autotot,time,t1,t0
      real*8 :: emean,emeantot
      real*8 :: xnorm1,xnorm2,xnorm1tot,xnorm2tot

      
*********************************************************
      namelist /inputtime/ntimes, nloop,kminelastic
*********************************************************
c! Initialize MPI environment and get proc's ID and number of proc in
c! the partition.

      call MPI_INIT(ierr)
!      call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

      write(filename,'("sal."i3.3)')idproc
      open(6,file=filename,status='unknown')
      write(6,'(40("="),//)')
      write(6,'(10x,"MadWave3 version 6 ",//)')
      write(6,*)' output of proc. idproc= ',idproc,' of nproc= ',nproc
      write(6,'(/,40("="),//)')

!     initialization of data

      call input_grid
      call ini_absorcion
!      call  paralelizacion
      call pot0

!     determining basis

      call basis
      
      call  paralelizacion
      
!     reactants and products functions calculation

      call angular_functions

      call radial_functions01_read

      if(iprod.gt.0.and.npun1.ne.1)then

         call product_radialf_read

      endif

!     reading potential

      call pot2
      
!    dimensioning vectors for propagation,
!    calculating kinetic energy terms, etc

      call set_vectors_Hphi

!     preparing initial wave packet reading electric dipole transition

      if(iphoto.eq.0)then
                              ! initial collision wave packet
         call set_colini
      
      elseif(iphoto.eq.1.or.iphoto.eq.2.or.iphoto.eq.3)then
         if(iphoto.eq.1)call read_trans_dipole
         call dip_bnd
         do i=1,ntotproc(idproc)
            rpaqproc(i)=rpaq0(i)
         enddo
!      elseif(iphoto.eq.3)then
!         call read_coupling
!         call coup_bndgrid
!         do i=1,ntotproc(idproc)
!            rpaqproc(i)=rpaq0(i)
!         enddo
       else
         write(6,*)' iphoto =',iphoto,' out of range'
         call flush(6)
         call MPI_BARRIER(MPI_COMM_WORLD, IERR)
         stop
      endif

!     initialize total flux quantities

      call ini_flux
      
! for flux on 02 diatomics in collisions (iphoto=0) of photodissociation (iphoto>0)
      if(iprod.eq.2.and.npun1.gt.1)then
         call ini_transcoor
      endif

* fftw3 initialization

      call difs

**>> initialating propagation
      
      write(6,'(40("="),/,10x,"Chebysev integration",/,40("="))')
      write(6,*)
      write(6,*)'  inputtime'
      write(6,*)'  ---------'
      open(10,file='input.dat',status='old')
      read(10,nml = inputtime)
      write(6,nml = inputtime)
      call flush(6)
      close(10)

      iit0=0
      indt=0
      iloop0=0
      
      open(3,file='cont.data',status='old',err=1)
      read(3,*)iit0,indt,iloop0
      close(3)      
      write(6,*)' reading initial iit0= ',iit0,indt,iloop0
     &            ,' in file cont.data'
      call flush(6)
 1    continue
      
      if(iit0.eq.0)then
         do i=1,ntotproc(idproc)
            rflanz0proc(i)=rpaqproc(i)
         enddo
         call difs

         do i=1,ntotproc(idproc)
            rpaqproc(i)=(rHpaqrec(i)
     &                   -emindlt*rflanz0proc(i))/delta2
         enddo

         iit0=1

      else
         write(name,'("cont",i2.2,".paq")')idproc
         write(6,*)' Reading wvp for continuation in file= ',name
         write(6,*)'         at it= ',iit0
         write(6,*)
         open(4,file=name,status='unknown',form='unformatted')
         read(4)rflanz0proc
         read(4)rpaqproc
         close(4)
      endif

!     *>> checking average energy

      emean=0.d0
      xnorm1=0.d0
      xnorm2=0.d0
      do i=1,ntotproc(idproc)
         emean=emean+rflanz0proc(i)*rHpaqrec(i)
         xnorm1=xnorm1+rflanz0proc(i)**2
         xnorm2=xnorm2+rHpaqrec(i)**2
      enddo

      call MPI_REDUCE(emean,emeantot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(xnorm1,xnorm1tot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(xnorm2,xnorm2tot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)

      if(idproc.eq.0)then
         write(6,*)
         write(6,*)'    Mean energy (eV) and norm of initial wvp= '
     &        ,emeantot/conve1/8065.5d0

         
         write(6,*)
         call flush(6)
         emean=emeantot
      endif

***********************
** Main loop in time **
***********************

      write(6,1111)

      iloop0=iit0/ntimes+1
      time=dble(iit0)
      if( iwrt_wvp == 1)then
         iloop00=iloop0-1
         call plot_wvp(iloop00)
      endif
      iloop=iloop0
      it=iit0
      call check(time,xnorm1tot,it,iloop)
     
      do iloop=iloop0,nloop+iloop0

         t0=MPI_Wtime()
         if(idproc == 0)then
           if(iphoto.ge.1.or.iprod.eq.1.or.iwrt_reac_distri.eq.2)then
               write(name,'("Cvj.",i4.4)')iloop
               open(15,file=name,status='new',form='unformatted')
           endif
           if(iprod == 2)then
              write(name,'("Cvjprod.",i4.4)')iloop
              open(16,file=name,status='new',form='unformatted')
           endif
         endif

         do it = iit0+1,iit0+ntimes

            time=it
            call difs
            do i=1,ntotproc(idproc)
               call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
               fabs=absfr1(ir1)*absfr2(ir2)
               rHnorm=(rHpaqrec(i)-emindlt*rpaqproc(i))/delta2
               xxx=(2.d0*rHnorm-fabs*rflanz0proc(i))*fabs
               rflanz0proc(i)=rpaqproc(i)
               rpaqproc(i)=xxx
            enddo

            call CoefEcheby(it)
            
            if(npun1.gt.1)call totalflux_k(it,kminelastic)
            
            call Cvjflux_k(it,kminelastic)
            if(iprod.eq.2.and.npun1.gt.1)then
               call prodpaq
               call prodcvj
            endif
            call check(time,xnorm1tot,it,iloop)

            if(idproc == 0)then
              if(iphoto.ge.1.or.iprod.eq.1.or.iwrt_reac_distri.eq.2)then
                  write(15)it,Cvj
                  call flush(15)
              endif
              if(iprod == 2)then
                  write(16)it,Cvjprod
                  call flush(16)
               endif
            endif
            
         enddo                  ! do  it

         t1=MPI_Wtime()
         if(idproc.eq.0)write(6,*)' time per loop= ',t1-t0
         iit0=iit0+ntimes
         call flush(6)

         if(npun1.gt.1)call totalflux
         call wvpcont(iit0,iloop,indt)
         
         call MPI_BARRIER(MPI_COMM_WORLD, IERR)
         
         if(idproc.eq.0)call printprob(iit0,iloop,indt)
         if( iwrt_wvp == 1)then
           iloop00=iloop
           call plot_wvp(iloop00)
         endif
         if(idproc == 0)then
           if(iphoto.ge.1.or.iprod.eq.1)then
               close(15)
           endif
           if(iprod == 2)then
              close(16)
           endif
         endif
        
         call MPI_BARRIER(MPI_COMM_WORLD, IERR)

         if(xnorm1tot.lt.1.d-8)then
            write(6,*)' norm = ',xnorm1tot,' below threshold: stop '
            stop
         endif

      enddo  ! iloop,xnorm1tot

      call MPI_FINALIZE(ierr)

      stop
 1111 format(/,11x,'k',12x,'norm',/,30('*'),/)
      end program

!----------------------------------------------------------

      subroutine check(time,xnorm1tot,it,iloop)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_photoini_01y2
      use mod_colini_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2

      implicit none
      include "mpif.h"     
      double precision :: time,autocor,autotot
      real*8 :: emean,emeantot
      real*8 :: xnorm1,xnorm2,xnorm1tot
      integer :: it,iloop,ierr,i

!--norm
      xnorm1=0.d0
      do i=1,ntotproc(idproc)
         xnorm1=xnorm1+rpaqproc(i)**2
      enddo

      call MPI_REDUCE(xnorm1,xnorm1tot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xnorm1tot,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      if(idproc.eq.0)then
         write(6,*)it,xnorm1tot
      endif

      if(iphoto.ge.1)then
            autocor=0.d0
            do i=1,ntotproc(idproc)
               autocor=autocor+rpaq0(i)*rpaqproc(i)
            enddo
            call MPI_REDUCE(autocor,autotot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
            if(idproc.eq.0)then
                write(ifileauto,*)it,autotot
                call flush(ifileauto)
            endif
      endif ! iphoto=1

      return
      end subroutine check

!---------------------------------------------------------------------

      subroutine CoefEcheby(ikcheb)
      use mod_flux_01y2
      use mod_pot_01y2, only : emindlt,delta2
      use mod_gridYpara_01y2, only : hbr 
      implicit none

      complex*16 :: zexpo,zfactor
      real*8 :: d,E,Es,expo,deno
      integer :: ie,ikcheb

!     for the representation of Energy resolved quantitites from a Chevishev expansion
! zCkcheby is the factor for Energy resolved reaction probabilities and state2state prob.
!    summing on the Chebyshev iterations.

      d=2.d0
      if(ikcheb.eq.0)d=1.d0
      zfactor=dcmplx(d*hbr/delta2,0.d0)
      do ie=1,nEtot
         E=etotS2(ie)
         Es=(E-emindlt)/delta2
         expo=-dble(ikcheb)*dacos(Es)
         zexpo=cdexp(dcmplx(0.d0,expo))
         deno=1.d0-Es*Es
         deno=dsqrt(deno)
         zCkcheby(iE)=zfactor*zexpo/dcmplx(deno,0.d0)
      enddo

      return
      end subroutine CoefEcheby

!---------------------------------------------------------------------

      subroutine wvpcont(iit0,iloop,indt)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_photoini_01y2
      use mod_colini_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      use mod_flux_01y2
      use mod_coortrans01_02
      implicit none
      real*8 :: abstot
      integer :: iit0,iloop,indt
* printing only if continuation files are needed

      if(ncontfile.eq.1)then

**>> continuation files

         write(name,'("cont",i2.2,".paq")')idproc
         open(4,file=name,status='unknown',form='unformatted')
         write(4)rflanz0proc
         write(4)rpaqproc
         close(4)

         write(name,'("cont",i2.2,".S2")')idproc
         open(4,file=name,status='unknown',form='unformatted')
         write(4)zS2
         write(4)zCR
         write(4)zCdR

         close(4)

         if(iprod.gt.1)then
            write(name,'("cont",i2.2,".S2prod")')idproc
            open(4,file=name,status='unknown',form='unformatted')
            write(4)zS2prod
            close(4)
         endif
         
!         if(nofe.gt.0)then
!               write(name,'("cont",i2.2,".zfE")')idproc
!               open(4,file=name,status='unknown',form='unformatted')
!               write(4)zfE
!               close(4)
!         endif

         if(idproc.eq.0)then
            open(3,file='cont.data',status='unknown')
            write(3,*)iit0,indt+1,iloop
            write(3,*)abstot,nflux,ncan
            close(3)

         endif  ! idproc.eq.0

      endif !cont

      return
      end
!---------------------------------------------------------------------

      subroutine printprob(iit0,iloop,indt)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_photoini_01y2
      use mod_colini_01y2
      use mod_absorcion_01y2
      use mod_Hphi_01y2
      use mod_flux_01y2
      use mod_coortrans01_02
      implicit none
      include "mpif.h"
      integer :: ifile,iii,ielec,ivprod,iomprod,ican,iele,j,iv
      integer :: ie,iom,ifilelec,ivfile,iiielec,iiiv
      integer :: iit0,iloop,indt
      real*8 :: S2reac,Av,S2no,S2nobis
      real*8 :: vibprod(nviniprod:nvmaxprod)
      real*8 :: rotdistri(jini:jmax)

      real*8 :: S2pro(nviniprod:nvmaxprod,jiniprod:jmaxprod
     &                      ,iomminprod:iommaxprod)
      complex*16 :: zzz
 
      ifilelec=50
      ivfile=100

*     writes information at each loop in time

      if(npun1.gt.1)then
         write(name,'("S2prod.v",i2.2,".J",i3.3,".k",i5.5)')
     &              nvref,Jtot,iloop
         open(20,file=name,status='unknown')
      endif
      if(iwrt_reac_distri.eq.1)then
!         write(name,'("S2mat.J",i3.3,".k",i5.5)')
!     &              Jtot,iloop
!     open(20,file=name,status='unknown')

         do ielec=1,nelec
            iiielec=ifilelec+ielec
            write(name,"('distriS2reac.elec',i2.2)")ielec
            open(iiielec,file=name,status='unknown')
         enddo
         iiiv=ivfile
         do ielec=1,nelec
            do iv=nvini,noBCstates(jini,ielec)
               iiiv=ivfile+(iv+1-nvini)+(nvmax-nvini+1)*(ielec-1)
               write(name,"('distriS2reac.v',i2.2,'.e',i1.1)")
     &                 iv,ielec
               open(iiiv,file=name,status='unknown')
            enddo
         enddo

      endif

      do ie=1,netot

****> state-2-state for reactants

         rotdistri(:)=0.d0
         S2reac=0.d0
         do ican=1,ncan
            iele=nelebas(ican)
            do j=j00(ican),jmax,inc
               do iv=nvini,noBCstates(j,iele)
                  zzz=zS2(ie,iv,j,ican)/2.d0/pi
                  Av=dreal(zzz*dconjg(zzz))
                  S2(iv,j,ican)=Av*S2factor(ie,iv,j,iele)
                  S2reac=S2reac+S2(iv,j,ican)
                  rotdistri(j)=rotdistri(j)+S2(iv,j,ican)
               enddo
            enddo
         enddo

!****> state-2-state for products

         if(npun1.gt.1)then
            if(iprod.eq.1)then
               S2no=0.d0
               vibprod(:)=0.d0
               do ican=1,ncan
               iele=nelebas(ican)
               do j=j00(ican),jmax,inc
               do iv=nvini,min0(nvmaxprod,noBCstates(j,iele))
                  vibprod(iv)=vibprod(iv)+S2(iv,j,ican)
               enddo
               enddo
               enddo
            elseif(iprod.gt.1)then
               S2no=0.d0
               vibprod(:)=0.d0
               S2pro(:,:,:)=0.d0
               do iom=iomminprod,iommaxprod
               do j=jiniprod,jmaxprod
                  do iv=nviniprod,nvmaxprod
                     zzz=zS2prod(ie,iv,j,iom)
                     Av=dreal(zzz*dconjg(zzz))*0.25d0/(pi*pi)
                     S2pro(iv,j,iom)=Av*S2prodfac(ie,iv,j)
                     S2no=S2no+S2pro(iv,j,iom)
                     vibprod(iv)=vibprod(iv)+S2pro(iv,j,iom)
                  enddo
               enddo
               enddo
            endif

***> printing total flux to products  

            do iv=nviniprod,nvmaxprod
               if(vibprod(iv).lt.1.d-90)vibprod(iv)=0.d0
            enddo
            if(S2prodtot(ie).lt.1.d-90)S2prodtot(ie)=0.d0
            if(S2reac.lt.1.d-90)S2reac=0.d0

            if(iprod.eq.0)then
               write(20,"(501(1x,e15.7))")etotS2(ie)/conve1/8065.5d0
     &             ,S2prodtot(ie)*photonorm
     &             ,(S2prodtot(ie)+S2reac)*photonorm
     &             ,reacfct(ie)
            elseif(iprod.eq.1)then
               write(20,"(501(1x,e15.7))")etotS2(ie)/conve1/8065.5d0
     &             ,S2reac*photonorm
     &             ,(S2prodtot(ie)+S2reac)*photonorm
     &          ,(vibprod(iv)*photonorm,iv=nvini,min0(nvmaxprod,nvmax))
            else
               write(20,"(501(1x,e15.7))")etotS2(ie)/conve1/8065.5d0
     &             ,S2prodtot(ie)*photonorm
     &             ,(S2prodtot(ie)+S2reac)*photonorm
     &             ,S2no*photonorm
     &            ,(vibprod(iv)*photonorm,iv=nviniprod,nvmaxprod)

            endif
         endif  ! npun.gt.1
         if(iwrt_reac_distri.eq.1)then


!            write(20,"(501(1x,e15.7))")etotS2(ie)/conve1/8065.5d0
!     &             ,S2reac*photonorm
!     &            ,(rotdistri(j)*photonorm,j=jini,jmax)
            
            iiiv=ivfile
            do ielec=1,nelec
               S2nobis=0.d0
            do iv=nvini,noBCstates(jini,ielec)
               iiiv=ivfile+(iv+1-nvini)+(nvmax-nvini+1)*(ielec-1)
               rotdistri(:)=0.d0
               do ican=1,ncan
                  iele=nelebas(ican)
                  if(iele.eq.ielec)then
                     do j=j00(ican),jmax,inc
                     zzz=zS2(ie,iv,j,ican)/2.d0/pi
                     Av=dreal(zzz*dconjg(zzz))
                     rotdistri(j)=rotdistri(j)+S2(iv,j,ican)
                     S2nobis=S2nobis+S2(iv,j,ican)
                     enddo ! j
                  endif ! iele=ielec
               enddo ! ican
               write(iiiv,"(1000(1x,e15.7))")etotS2(ie)/conve1/8065.5d0
     &            ,(rotdistri(j)*photonorm,j=jini,jmax)

            enddo               ! iv
                iiielec=ifilelec+ielec
               write(iiielec,"(1000(1x,e15.7))")
     &           etotS2(ie)/conve1/8065.5d0 ,S2nobis*photonorm
 
            enddo ! ielec

         endif
      enddo ! ie=1,ne
   
!      close(20)

      do ielec=1,nelec
         iiielec=ifilelec+ielec
         close(iiielec)
      enddo
      iiiv=ivfile
      do ielec=1,nelec
         do iv=nvini,noBCstates(jini,ielec)
            iiiv=ivfile+(iv+1-nvini)+(nvmax-nvini+1)*(ielec-1)
            close(iiiv)
         enddo
      enddo

      return
      end       

!---------------------------------------------------------------------

      subroutine plot_wvp(iloop)
      use mod_flux_01y2
      use mod_pot_01y2, only : emindlt,delta2
      use mod_gridYpara_01y2, only : hbr 
      use mod_flux_01y2
      implicit none
      include "mpif.h"
      integer :: iloop,jelec,jr2,ir1,iang,ir2,i,icanp,ielec,iangp
      integer :: nnn,ierr,ir,iom
      real*8 :: r2,r1
      real*8 :: fun(npun1,nangu), funtot(npun1,nangu)

* writting the wvp at each iloop
* in the same angular grid as the pes
* (i.e. each nangplot points of the angular grid)

      do jelec=1,nelec
         if(idproc.eq.0)then   
            write(name,'("wvp.elec",i2.2,".i",i3.3)')jelec,iloop
            open(10,file=name,status='unknown')
         endif
      
         do jr2=1,npun2
         r2=rmis2+dble(jr2-1)*ah2
         if(r2.lt.absr2)then
            do ir1=1,npun1
            do iang=1,nangu
               fun(ir1,iang)=0.d0
               funtot(ir1,iang)=0.d0
            enddo
            enddo

            do i=1,ntotproc(idproc)
               call indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
               iang=indangreal(iangp,idproc)
               if(ir2.eq.jr2.and.ielec.eq.jelec)then
                fun(ir1,iang)=fun(ir1,iang)+rpaqproc(i)*rpaqproc(i)
               endif
            enddo

            nnn=npun1*nangu
            call MPI_REDUCE(fun,funtot,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)

            if(idproc.eq.0)then
               do ir1=1,npun1,n1plot
                  r1=rmis1+dble(ir1-1)*ah1
                  if(r1.le.absr1)then
                     do iang=1,nangu
                        if(dabs(funtot(ir1,iang)).lt.1.d-90)then
                           funtot(ir1,iang)=0.d0
                        endif
                     enddo
                     write(10,'(500(1x,e15.7))')r1,r2
     &                 ,(funtot(ir1,iang),iang=1,nangu,nangplot)
                  endif
               enddo
           
               write(10,'()')

            endif

         endif ! r2.le.r2abs
         enddo ! ir2

         if(idproc.eq.0) close(10)
      enddo  ! jelec

      return
      end subroutine plot_wvp
!---------------------------------------------------------------------
