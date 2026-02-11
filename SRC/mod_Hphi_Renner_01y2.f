      module mod_Hphi_01y2
!
!     Initialize magnitudes for the application of H phi
!     for a triatomic  01-2 system
!     in Jacobi coordinates in a body-fixed frame
!
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      
      implicit none

      save

      real*8, allocatable :: rpaqproc(:)
     &         ,rflanz0proc(:)
     &         ,rHpaqsend(:),rHpaqrec(:),recibo(:)
     &     ,rHpaqproc(:,:)
     &     ,rpaq0(:)

!      real*8, allocatable :: xl2p(:,:,:,:,:)
!     &     , xj2(:,:,:,:)
!---> Adding extra dimension in nelecmax for Renner coupling
      real*8, allocatable :: xl2p(:,:,:,:,:,:)
     &     , xj2(:,:,:,:,:)
      
      real*8, allocatable :: rr1m2(:,:),rr2m2(:,:)

      real*8, allocatable :: xpr1(:),xp2r1(:),xpr2(:),xp2r2(:)
     &               , pr1(:,:),p2r1(:,:)
     &                ,pr2(:,:),p2r2(:,:)
      integer, allocatable :: nopr1(:),nopr2(:)
      
      contains
********************************************
*   functions of   mod_Hphi_01y2.f     *
********************************************
!=======================================================

!--------------------------------------------------
      subroutine set_vectors_Hphi
!
!     allocate vectors required for the propagation and application of H phi
!
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      implicit none
      include "mpif.h"
      
      integer :: ierror
      integer*8 :: imem

      imem=nbastotproc*(6+ncouprocmax)
     &    +npun2*nangu*nanguprocdim*ncanprocdim*3
     &    +npun1*nangu*nanguprocdim*ncanmax
     &    +npuntot*nangu*2
     &    +npun1*2*(1+nangu)+npun2*2*(1+nangu)+2*nangu

      write(6,*)' memory to allocate in set_vectors_Hphi '
     &    ,'  in proc= ',idproc,' is= '
     &    ,dble(imem)*1.d-9,' Gb'
      call flush(6)
      
! paqs and Hpaqs

      allocate(rpaqproc(nbastotproc),recibo(nbastotproc)
     &         ,rflanz0proc(nbastotproc)
     &         ,rpaq0(nbastotproc)
     &         ,rHpaqsend(nbastotproc),rHpaqrec(nbastotproc)
     &         ,rHpaqproc(nbastotproc,ncouprocmax)
     &         ,stat=ierror)
      if(ierror.ne.0)then
        write(*,*)" error in initmem for paqs and Hpaqs "
        stop
      else
         imem= nbastotproc*(6+ncouprocmax)
         norealproc_mem=norealproc_mem+imem
         write(6,*)' rpaq,rHpaq,etc in set_vectors_Hphi, proc= '
     &       ,idproc
     &       ,' imem=',imem
     &       ,' memory(Gb)= '
     &        ,dble(imem*8)*1.d-9
         call flush(6)
      endif
      rpaq0(:)=0.d0

! centrifugal terms
!      allocate(xl2p(npun2,nangu,nanguprocdim,ncanprocdim,-1:1)
!     &       , xj2(npun1,nangu,nanguprocdim,ncanmax)
!     &       , stat=ierror)
!---> Adding extra dimension in nelecmax for Renner coupling
      allocate(xl2p(npun2,nangu,nanguprocdim,ncanprocdim,-1:1,nelecmax)
     &       , xj2(npun1,nangu,nanguprocdim,ncanmax,ncanmax)
     &       , stat=ierror)
      if(ierror.ne.0)then
        write(*,*)" error in initmem for centrifugal "
        stop
      else
         norealproc_mem=norealproc_mem+imem
         imem=
     &    npun2*nangu*nanguprocdim*ncanprocdim*3
     &        +npun1*nangu*nanguprocdim*ncanmax
         imem=imem*nelecmax
          write(6,*)' xl2p,xj2 in set_vectors_Hphi, proc= '
     &       ,idproc
     &       ,' imem=',imem
     &       ,' memory(Gb)= '
     &       ,dble(imem*8)*1.d-9
         call flush(6)
      endif

!     grid 

      allocate(rr1m2(npuntot,nangu)
     &     ,rr2m2(npuntot,nangu)
     &       , stat=ierror)
      if(ierror.ne.0)then
        write(*,*)" error in initmem for Ri^-2 "
        stop
      else
         norealproc_mem=norealproc_mem
     &    +npuntot*nangu*2
         write(6,*)'norealproc_mem= ',norealproc_mem
         write(6,*)'nointegerproc_mem= ',nointegerproc_mem
         call flush(6)
      endif

!     momenta
      
      allocate(xpr1(npun1),xp2r1(npun1)
     &               ,xpr2(npun2),xp2r2(npun2)
     &               , pr1(npun1,nangu)
     &                ,p2r1(npun1,nangu)
     &                ,pr2(npun2,nangu)
     &                ,p2r2(npun2,nangu)
     &                ,nopr1(nangu)
     &                ,nopr2(nangu)
     &       , stat=ierror)
      if(ierror.ne.0)then
        write(*,*)" error in initmem for momenta "
        stop
      else
         norealproc_mem=norealproc_mem
     &    +npun1*2*(1+nangu)+npun2*2*(1+nangu)+2*nangu
      endif

!     kinetic terms initialization
      
      write(6,'(/,"-",/,3x," Initializing Kinetic terms ",/)')

**>> Angular kinetic terms 

      call  angkinj2
      call bfgridl2mat

**>> grid dependent radial  kinetic terms

      call Tradial
      
      write(6,*)' ending set_vectors_Hphi, proc= ',idproc
     &     ,' norealproc_mem=',norealproc_mem
     &     ,' nointegerproc_mem=',nointegerproc_mem
     &     ,' memory(Gb)= '
     &   ,dble(nointegerproc_mem*4+norealproc_mem*8)*1.d-9
      call flush(6)

      return
      end subroutine set_vectors_Hphi
!--------------------------------------------------
!--------------------------------------------------
      subroutine angkinj2
!-------------------------------------------------!
!     j² terms represented in the angular grid    !
!             for evaluation of H phi             !
!                                                 !
!       1)  <theta,Omega| j^ 2 | theta',Omega>    !
!                                                 !
!                 Using the DVR  method           !
!_________________________________________________!
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      implicit none
      include "mpif.h"
      
      integer :: ierror
      integer :: jjjmax,iang,jangp,icanp,ir1,ican,ielec,iom,iomdi,iomat
      integer :: jcanp,jcan,jelec,jom,jomdi,jomat
       integer :: i,j,jp,k,lamj,kp
      real*8 :: xj,xom,r1,yyy,x1,x2,xLam,xLAMj

      real*8 :: xkinj(nangu,nangu)

      if(nangu.gt.1)then
         jjjmax=(nangu-1)*inc

         xj2(:,:,:,:,:)=0.d0


         do jcanp=1,ncanproc(idproc)
            jcan=ibasproc(jcanp,idproc)
            jelec=nelebas(jcan)
            jom=iombas(jcan)
            jomdi=iomdiat(nelebas(jcan))
            jomat=iomatom(nelebas(jcan))
      
    
            
         do icanp=1,ncanproc(idproc)
            ican=ibasproc(icanp,idproc)
            ielec=nelebas(ican)
            iom=iombas(ican)
            iomdi=iomdiat(nelebas(ican))
            iomat=iomatom(nelebas(ican))
            xom=dble(iom)  

**>> Evaluation of the matrix for j^2
            xkinj(:,:)=0.d0
            do i=1,nojbas(ican)
            do j=1,nojbas(jcan)
               xkinj(i,j)=0.d0
               if(ican.eq.jcan.and.i.eq.j)then
                  xj=dble(jbas(i,ican))
                  xLam=dble(iom-iomdi)  ! (iomdi)--> renner
                  xkinj(i,j)=xj*(xj+1.d0)+xlam*xlam !-dble((iom-iomdi)**2) ! xLam*xLam
                  xkinj(i,j)=xkinj(i,j)*hbr*hbr*0.5d0 ! change for Renner
!               else if(iom.eq.jom.and.jbas(i,ican).eq.jbas(j,jcan)
!     &                .and. ielec.ne.jelec. and. jomat.eq. -iomat)then
!                  xLam=dble(iomdi)  ! (iomdi)--> renner
!                  xom=dble(iom)  
!                  xkinj(i,j)=0.d0 !-2.d0*xLam*(xom-xLam)*hbr*hbr*0.5d0 ! 
               end if
            enddo
            enddo

            do ir1=1,npun1
               r1=rmis1+dble(ir1-1)*ah1
               do k=1,nojbas(ican)
               do kp=1,nojbas(jcan)
                  if(jbas(k,ican).le.jjjmax)then
                     yyy=xkinj(k,kp)/(xm1red*r1*r1)
                     if(yyy.gt.rotcutmax)then
                         yyy=rotcutmax
                     endif
                     do i=1,nangu
                     do jp=1,nanguproc
                      j=indangreal(jp,idproc)
                      x1=Djmmp(i,jbas(k,ican),ielec,iom)
                      x2=Djmmp(j,jbas(kp,jcan),jelec,jom)
                      xj2(ir1,i,jp,icanp,jcanp) =
     &                            xj2(ir1,i,jp,icanp,jcanp) + yyy*x1*x2
                     enddo
                     enddo
                  endif
                enddo
                 enddo
           enddo

         enddo                  ! icanp
         enddo                  !jcanp
      else
         xj2(:,:,:,:,:)=0.d0
      endif
         
      return
      end subroutine angkinj2
      
!--------------------------------------------------
!--------------------------------------------------
      subroutine bfgridl2mat
!-------------------------------------------------!
!      l² terms represented in the angular grid   !
!             for evaluation of H phi             !
!        1)  <theta,Omega| l^ 2 | theta',Omega>   !
!        2)  <theta,Omega| l^ 2 | theta',Omega'>  !
!                                                 !
!                 Using the DVR  method           !
!_________________________________________________!
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      implicit none
      include "mpif.h"
      
      integer :: ierror

      real*8,allocatable ::  bfl2mat(:,:),xl2mat(:,:)
     &                      ,hmat(:,:),T(:,:),eigval(:)
      integer,allocatable :: ind(:)
      integer :: iomdim1,iomdim0,maxomg
      integer :: ielec,jelec,j,iom0,iom1,nnn,i,jj,jom,ir2
      real*8 :: yylmin,yylmax,yylmean,r2,r2eff,emineff,emeaneff
      real*8 :: eee,emaxeff,x11,x22,yyy,facmass,erotlmax,xxx
      integer :: iom,iterm,jjcanp,jjcan,jjelec,jjom,iq,iang,jang,jangp
      integer :: sigelec,sigelecp,iomat,jomat

      
*  initialization

      iomdim1=iommax
      iomdim0=iommin
      maxomg=iomdim1-iomdim0+1

      allocate(bfl2mat(iomdim0:iomdim1,iomdim0:iomdim1)
     &       , xl2mat(iomdim0:iomdim1,iomdim0:iomdim1)
     &       , hmat(maxomg,maxomg),T(maxomg,maxomg),eigval(maxomg)
     &       , ind(maxomg)
     &       , stat=ierror)

      erotlmax=0.d0
      facmass=hbr*hbr*0.5d0/xm2red

      xl2p(:,:,:,:,:,:)=0.d0

*  starting loops

      do jelec=1,nelec
      do ielec=1,nelec
         sigelec=sigatom(ielec)
         sigelecp=sigatom(jelec)
         do j=j0,(nangu-1)*inc,inc
            iom0=iommin
            if(iparity.ne.sigelec*int((-1.d0)**Jtot)
     &         .and.iommin.eq.0)then
                 iom0=1
            endif
            if(iom0.le.j)then

            iom1=min0(j,iommax)
            nnn=iom1-iom0+1
            if(nnn.gt.0)then
               iomat=iomatom(ielec)*sigatom(ielec)
               jomat=iomatom(jelec)*sigatom(jelec)
               call l2mat_Renner(bfl2mat,facmass,iom0,iom1,j,Jtot
     &                 ,iomat,sigelec,sigelecp,iommin,iommax)
              
               if(nnn.gt.1)then
                  i=0
                  do iom=iom0,iom1
                     i=i+1
                     jj=0
                     do jom=iom0,iom1
                        jj=jj+1
                        hmat(i,jj)=bfl2mat(iom,jom)
                     enddo
                  enddo

                  call DIAGON(hmat,nnn,maxomg,T,eigval)
                  yylmin=eigval(1)
                  yylmax=eigval(nnn)
                  yylmean=0.5d0*(yylmin+yylmax)
               elseif(nnn.eq.1)then
                  yylmin=bfl2mat(iom0,iom0)
                  yylmax=yylmin
                  yylmean=yylmin
               endif
* Cut of rotational barrier to avoid numerical problems in diagonal representation

               do ir2=1,npun2
                  r2=rmis2+dble(ir2-1)*ah2

*  cutting centrifugal term
 
                  r2eff=r2
                  emineff=yylmin/(r2*r2)
                  emaxeff=yylmax/(r2*r2)
                  emeaneff=yylmean/(r2*r2)
                  if(emeaneff.gt.rotcutmax)then
                      eee=rotcutmax
                      r2eff=dsqrt(dabs(yylmean/eee))
                      emaxeff=yylmax/(r2eff*r2eff)
                  endif
                  if(emaxeff.gt.erotlmax)erotlmax=emaxeff

                  xl2mat(:,:)=0.d0
                  do iom=iom0,iom1
                  do jom=iom0,iom1
                     if(dabs(bfl2mat(iom,jom)).gt.1.d-8
     &                    .and.r2eff.gt.1.d-3)then
                         xl2mat(iom,jom)=bfl2mat(iom,jom)/(r2eff*r2eff)
                     endif
                  enddo
                  enddo
                       
* Transformation to angular grid representation

                  do iom=iom0,iom1
                  do jom=iom0,iom1
                     iterm=0
                     do jjcanp=1,ncanproc(idproc)

                        jjcan=ibasproc(jjcanp,idproc)
                        jjelec=nelebas(jjcan)
                        jjom=iombas(jjcan)

                        iq=-10
                        if(jjelec.eq.ielec.and.jjom.eq.jom
     &                     .and.iabs(jom-iom).le.1)then

                           if(iom.eq.jom)then
                              iq=0
                           elseif(jom.eq.iom+1)then
                              iq=+1
                           elseif(jom.eq.iom-1)then
                              iq=-1
                           endif
                           if(iabs(iq).le.1)then
                           do iang=1,nangu
                           do jangp=1,nanguproc
                              jang=indangreal(jangp,idproc)
                              x11=Djmmp(iang,j,ielec,iom)
                              x22=Djmmp(jang,j,ielec,jom)
                              yyy=xl2mat(iom,jom)
                              xl2p(ir2,iang,jangp,jjcanp,iq,jjelec) =
     &                             xl2p(ir2,iang,jangp,jjcanp,iq,jjelec)
     &                            + x11*x22*yyy
                           enddo
                           enddo
                           endif

!     For Renner coupling
                        
!                        else if(jjelec.ne.ielec.and.jjom.eq.jom
!     &                     .and.iabs(jom-iom).le.0)then
!                        
!                           do iang=1,nangu
!                           do jangp=1,nanguproc
!                              jang=indangreal(jangp,idproc)
!                              x11=Djmmp(iang,j,ielec,iom)
!                              x22=Djmmp(jang,j,ielec,jom)
!                              yyy=xl2mat(iom,jom)
!                              iq=0
!                              xl2p(ir2,iang,jangp,jjcanp,iq,jjelec) =
!     &                             xl2p(ir2,iang,jangp,jjcanp,iq,jjelec)
!     &                            + x11*x22*yyy
!                           enddo
!                           enddo
                        endif
                     enddo  ! ican1p

                  enddo  ! jom
                  enddo  ! iom

               enddo ! ir2
            endif ! nnn > 0

            endif  ! iom0 < j
         enddo  ! j=0,jmax
      enddo  ! ielec=1,nelec
      enddo  ! jelec=1,nelec

!     deallocating auxiliary matrices
      deallocate(bfl2mat
     &       , xl2mat
     &       , hmat,T,eigval
     &       , ind)
      
      return
      end subroutine  bfgridl2mat
      
!--------------------------------------------------
      subroutine Tradial
!-------------------------------------------------!
!                                                 !
!                                                 !
!_________________________________________________!
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      implicit none
      include "mpif.h"
      
      integer :: ierror
      real*8 ::  rr1m2mx,rr2m2mx,r1,r2,xk1max,xk2max,xk1min,xk2min
      real*8 :: box,box1,box2
      integer :: ip,iangp,iang,icanp,ican,ielec
      integer :: ir1min,ir2min,ir1max,ir2max,ir,ir1,ir2
      integer :: n0,n1,n2
      real*8 :: xpr1(npun1),xp2r1(npun1)
      real*8 :: xpr2(npun2),xp2r2(npun2)
      

      rr1m2mx=0.d0
      rr2m2mx=0.d0

      nopr1(:)=1
      nopr2(:)=1
      do ip=0,nproc-1
         do iangp=1,nanguproc
            iang=indangreal(iangp,ip)

            do icanp=1,ncanproc(ip)
               ican=ibasproc(icanp,ip)

               ielec=nelebas(ican)
               ir1min=npun1
               ir2min=npun2
               ir1max=0
               ir2max=0
               do ir=1,npunreal(iang)
                  ir1=indr1(ir,iang)
                  ir2=indr2(ir,iang)
                  r1=dble(ir1-1)*ah1+rmis1
                  r2=dble(ir2-1)*ah2+rmis2
* for kinetic terms
                  rr1m2(ir,iang)=1.d0/r1/r1
                  rr2m2(ir,iang)=1.d0/r2/r2
                  if(rr1m2(ir,iang).gt.rr1m2mx)
     &                                   rr1m2mx=rr1m2(ir,iang)
                  if(rr2m2(ir,iang).gt.rr2m2mx)
     &                                    rr2m2mx=rr2m2(ir,iang)
                  if(ir1.lt.ir1min)ir1min=ir1
                  if(ir2.lt.ir2min)ir2min=ir2
                  if(ir1.gt.ir1max)ir1max=ir1
                  if(ir2.gt.ir2max)ir2max=ir2

               enddo
               nopr1(iang)=ir1max-ir1min+1
               nopr2(iang)=ir2max-ir2min+1

            enddo  ! ican

         enddo  ! iangproc
      enddo  ! ip 

**>> Radial kinetic terms: Momenta for Fourier Transforms

      xk1max=0.d0
      xk2max=0.d0
      xk1min=1.d10
      xk2min=1.d10

      do ielec=1,nelec
         do iang=1,nangu
* for R_1 = r
            n0=nopr1(iang)
            call noptFFT(n0,n1,npun1)
            if(n1.lt.npun1min)n1=npun1min
            nopr1(iang)=n1
            if(npun1.gt.1)then
               box=dble(n1-1)*ah1
               call sinmom(box,n1,npun1,xm1red,hbr,xpr1,xp2r1)
               do ir1=1,n1
                  pr1(ir1,iang)=xpr1(ir1)
                  p2r1(ir1,iang)=xp2r1(ir1)
                  if(dabs(xp2r1(ir1)).gt.xk1max)xk1max=dabs(xp2r1(ir1))
                  if(dabs(xp2r1(ir1)).lt.xk1min)xk1min=dabs(xp2r1(ir1))
               enddo
            else
               pr1(1,iang)=0.d0
               p2r1(1,iang)=0.d0
               xk1max=0.d0
               xk1min=0.d0
            endif
* for R_2 = R
            n2=npun2
            nopr2(iang)=n2
            box=dble(n2-1)*ah2
            call sinmom(box,n2,npun2,xm2red,hbr,xpr2,xp2r2)
            do ir2=1,n2
               pr2(ir2,iang)=xpr2(ir2)
               p2r2(ir2,iang)=xp2r2(ir2)
               if(dabs(xp2r2(ir2)).gt.xk2max)xk2max=dabs(xp2r2(ir2))
               if(dabs(xp2r2(ir2)).lt.xk2min)xk2min=dabs(xp2r2(ir2))

            enddo
c            if(idproc.eq.0)then   
c               write(6,*)'  Elec. St. = ',ielec,'  no. angle= ',iang
c               write(6,*)'     n1_FFT = ',n1,' while npun1= ',npun1 
c               write(6,*)'     n2_FFT = ',n2,' while npun2= ',npun2 
c               call flush(6)
c            endif
         enddo
      enddo
      
      write(6,*)
      write(6,*)'  ** check of energies covered by radial grids **'
      write(6,*)
      write(6,*)'        Vmin(eV)= ',vmintot/conve1/eV2cm
      write(6,*)'  E_kin^1(max)= ',xk1max/conve1/8065.5d0
     &            ,' E_kin^2(max)= ',xk2max/conve1/8065.5d0
      write(6,*)'  E_kin^1(min)= ',xk1min/conve1/8065.5d0
     &            ,' E_kin^2(min)= ',xk2min/conve1/8065.5d0
 
      write(6,*)
     
***> for first derivatives: Fourier transform instead of Sinus transform 

       box1=dble(npun1-1)*ah1
       box2=dble(npun2-1)*ah2
       call fftmom(box1,npun1,npun1,xm1red,hbr,xpr1,xp2r1)
       call fftmom(box2,npun2,npun2,xm2red,hbr,xpr2,xp2r2)
     
      return
      end subroutine Tradial
      
!--------------------------------------------------
      subroutine difs
*********************************************************************
**                  Subroutine DIFS                                **
**                                                                 **
**                Evaluates  H Phi(t)                              **
**            using Reatant Jacobi coordinates                     **
**                                                                 **
*********************************************************************

      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      implicit none
      include "mpif.h"
      
      integer :: ierror

      INTEGER FFTW_R2HC
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_DFT_R2HC_ICKY
      PARAMETER (FFTW_DFT_R2HC_ICKY=512)
      INTEGER FFTW_NONTHREADED_ICKY
      PARAMETER (FFTW_NONTHREADED_ICKY=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)
      
      integer*8 plan,flags
      INTEGER :: req, status(MPI_STATUS_SIZE),ierr
      
      real*8 :: xp2r2p2r1(npun2,npun1),rphi0(npun2,npun1)
      real*8 :: auxR(npun1*npun2),raux(npun2,npun1)
      integer :: icanp,ican,ielec,iom,jangp,jang,ir1,ir2,ir
      integer :: i,nn2,nn1,iii,i1,i2,ican2p,ican2,ielec2,iom2
      integer :: itotp,jdproc,iangp,iang,jdprocc,nsend,nreceived
      integer :: semaforo
      integer*8 :: max_semaforo,icount,i0
      real*8 :: divi,ekin,r1,r2,renner

      max_semaforo=10000000
      flags=FFTW_ESTIMATE
      plan=0

      rHpaqproc(:,:)=0.d0
      rHpaqrec(:)=0.d0
                                 
******************
**  1) T Phi(t) **
******************
      do icanp=1,ncanproc(idproc)
         ican=ibasproc(icanp,idproc)
         ielec=nelebas(ican)
         iom=iombas(ican)

!$OMP parallel do schedule(dynamic)
!$OMP&     reduction(+:rHpaqrec,rHpaqproc)
!$OMP&     private(jangp,jang,ir1,ir2,rphi0,raux,xp2r2p2r1,i0,ir,i
!$OMP&            ,divi,ekin,iii,itotp,ican2p,ican2
!$OMP&            ,ielec2,iom2,jdprocc,jdproc,iang,iangp
!$OMP&            ,i1,i2,nn1,nn2,auxr,plan,icount)       
!$OMP&    shared(indangreal,rpaqproc,p2r1,p2r2,radcutmax
!$OMP&             ,vvv,ibasproc,iombas,nelebas,xj2,xl2p
!$OMP&             ,indr1,indr2,ncouproc,ncanproc,nanguproc
!$OMP&             ,ican,icanp,ielec,iom,idproc,npun1,npun2
!$OMP&             ,flags,indtotproc,semaforo   
!$OMP&             ,npunreal,nopr1,nopr2,ipcou,max_semaforo)
         do jangp=1,nanguproc
            jang=indangreal(jangp,idproc)

             
             do ir1=1,npun1
                do ir2=1,npun2
                   rphi0(ir2,ir1)=0.d0
                   raux(ir2,ir1)=0.d0
                   xp2r2p2r1(ir2,ir1)=0.d0
                enddo
             enddo

             ir=1
             i0=indtotproc(ir,icanp,jangp,idproc)-1
             do ir=1,npunreal(jang)
                i=i0+ir  !indtotproc(ir,icanp,jangp,idproc)
                ir1=indr1(ir,jang)
                ir2=indr2(ir,jang)
                raux(ir2,ir1)=rpaqproc(i)
                rphi0(ir2,ir1)=raux(ir2,ir1)
             enddo

             nn2=nopr2(jang)
             nn1=nopr1(jang)

             divi=4.d0*dble((nn1+1))*dble((nn2+1))
             divi=1.d0/(divi)

             do ir1=1,nn1
             do ir2=1,nn2
                ekin=p2r1(ir1,jang)+p2r2(ir2,jang)
                if(ekin.gt.radcutmax)ekin=radcutmax
                xp2r2p2r1(ir2,ir1)=ekin
                iii=nn2*(ir1-1)+ir2   ! iii+1
                auxR(iii)=raux(ir2,ir1)
             enddo
             enddo
!$OMP critical
             semaforo=1

             call dfftw_plan_r2r_2d(plan,nn2,nn1,auxR,auxR
     &                    ,FFTW_RODFT00
     &                    ,FFTW_RODFT00
     &            ,flags)
             
             semaforo=0
!$OMP end critical
             icount=0
 1           continue
             if(semaforo.eq.1)then
!                write(6,*)'  atrapdado en if ',jang
!                call flush(6)
                icount=icount+1
                if(mod(icount,max_semaforo).eq.0)then
                   write(6,*)'  icount in difs= ',icount,max_semaforo
                   call flush(6)
!                   stop
                endif
                go to 1
             endif
             
             call dfftw_execute(plan)


             do ir1=1,nn1
             do ir2=1,nn2
                iii=nn2*(ir1-1)+ir2  ! iii+1
                auxR(iii)=auxR(iii)*(xp2r2p2r1(ir2,ir1))*divi
             enddo
             enddo
             call dfftw_execute(plan)

!$OMP critical
             semaforo=1
             call dfftw_destroy_plan(plan)
             semaforo=0
!$OMP end critical

             do ir1=1,nn1
             do ir2=1,nn2
                iii=nn2*(ir1-1)+ir2  !iii+1
                raux(ir2,ir1)=auxR(iii)
             enddo
             enddo

             ir=1
             i0=indtotproc(ir,icanp,jangp,idproc)-1
             do ir=1,npunreal(jang)
                ir1=indr1(ir,jang)
                ir2=indr2(ir,jang)
                itotp=i0+ir !indtotproc(ir,icanp,jangp,idproc)
                rHpaqrec(itotp)=rHpaqrec(itotp)+raux(ir2,ir1)
             enddo

******************
**  2) V Phi(t) **
******************

            do ican2p=1,ncanproc(idproc)
                ican2=ibasproc(ican2p,idproc)
                ielec2=nelebas(ican2)
                iom2=iombas(ican2)
                
!                renner=0.d0
                if(iom.eq.iom2)then
                   ir=1
                   i0=indtotproc(ir,ican2p,jangp,idproc)-1

!                   if(ielec.eq.ielec2)then
!                      renner=dble(iomatom(ielec)*iomatom(ielec))
!                      renner=renner*hbr*hbr*0.5d0
!                   elseif(iomatom(ielec).eq.iomatom(ielec2)
!     &         .and.dabs(sigatom(ielec)+sigatom(ielec2)).lt.1.d-5)then
!                      renner=dble(iabs(iomatom(ielec))*iom)
!                      renner=renner*hbr*hbr*0.5d0
!                   endif
                   do ir=1,npunreal(jang)  ! assumes that grids are the same for all electronic states
                        itotp=i0+ir !indtotproc(ir,ican2p,jangp,idproc)
                        ir1=indr1(ir,jang)
                        ir2=indr2(ir,jang)
                        r1=rmis1+dble(ir1-1)*ah1
                        r2=rmis2+dble(ir2-1)*ah2
                        rHpaqrec(itotp)=rHpaqrec(itotp)
     &                       + rphi0(ir2,ir1)*vvv(ir,jangp,ielec,ielec2)
!     &                    +rphi0(ir2,ir1)*renner*( 1.d0/(xm2red*r2*r2) )
!---     &                      +renner*( 1.d0/(xm1red*r1*r1) )
                   enddo
                endif
            enddo

******************************
**>>  3) j^2 + l^2
******************************

            do jdprocc=1,ncouproc(idproc)
            jdproc=ipcou(jdprocc,idproc)
            do iangp=1,nanguproc
               iang=indangreal(iangp,jdproc)
               do ican2p=1,ncanproc(jdproc)
                  ican2=ibasproc(ican2p,jdproc)   
                  iom2=iombas(ican2)
                  ielec2=nelebas(ican2)
                  if(iom.eq.iom2.and.ielec.eq.ielec2)then
                     ir=1
                     i0=indtotproc(ir,ican2p,iangp,jdproc)-1
                     do ir=1,npunreal(iang)
                        ir1=indr1(ir,iang)
                        ir2=indr2(ir,iang)
                        itotp=i0+ir ! indtotproc(ir,ican2p,iangp,jdproc)
                      rHpaqproc(itotp,jdprocc)=rHpaqproc(itotp,jdprocc)
     &                 +rphi0(ir2,ir1)*(xj2(ir1,iang,jangp,icanp,ielec2)
     &                             +xl2p(ir2,iang,jangp,icanp,0,ielec2))
                     enddo
                  elseif(iom.eq.iom2+1.and.ielec.eq.ielec2)then
                     ir=1
                     i0=indtotproc(ir,ican2p,iangp,jdproc)-1
                     do ir=1,npunreal(iang)
                        ir1=indr1(ir,iang)
                        ir2=indr2(ir,iang)
                        itotp=i0+ir !indtotproc(ir,ican2p,iangp,jdproc)
                     rHpaqproc(itotp,jdprocc)=rHpaqproc(itotp,jdprocc)
     &               +rphi0(ir2,ir1)*xl2p(ir2,iang,jangp,icanp,1,ielec2)
                     enddo
                  elseif(iom.eq.iom2-1.and.ielec.eq.ielec2)then
                     ir=1
                     i0=indtotproc(ir,ican2p,iangp,jdproc)-1
                     do ir=1,npunreal(iang)
                        ir1=indr1(ir,iang)
                        ir2=indr2(ir,iang)
                        itotp=i0+ir !indtotproc(ir,ican2p,iangp,jdproc)
                     rHpaqproc(itotp,jdprocc)=rHpaqproc(itotp,jdprocc)
     &             +rphi0(ir2,ir1)*xl2p(ir2,iang,jangp,icanp,-1,ielec2)
                     enddo
                  elseif(iom.eq.iom2.and.ielec.ne.ielec2)then
                     ir=1
                     i0=indtotproc(ir,ican2p,iangp,jdproc)-1
                     do ir=1,npunreal(iang)
                        ir1=indr1(ir,iang)
                        ir2=indr2(ir,iang)
                        itotp=i0+ir !indtotproc(ir,ican2p,iangp,jdproc)
                     rHpaqproc(itotp,jdprocc)=rHpaqproc(itotp,jdprocc)
     &               +rphi0(ir2,ir1)*xl2p(ir2,iang,jangp,icanp,0,ielec2)
     &               +rphi0(ir2,ir1)*xj2(ir1,iang,jangp,icanp,ican2p)
                     enddo
                  endif
               enddo ! ican2p

            enddo  ! iangp
            enddo  ! jdprocc

         enddo   ! jangp

      enddo ! icanp

** sending all parts to the corresponding processor

      do jdprocc=1,ncouproc(idproc)
         jdproc=ipcou(jdprocc,idproc)

         if(jdproc.eq.idproc)then

            do itotp=1,ntotproc(idproc)
              rHpaqrec(itotp)=rHpaqrec(itotp)+rHpaqproc(itotp,jdprocc)
            enddo

         else

             nsend=ntotproc(jdproc)
             nreceived=ntotproc(idproc)
             do itotp=1,ntotproc(jdproc)
                rHpaqsend(itotp)=rHpaqproc(itotp,jdprocc)
             enddo
          
             call MPI_IRECV(recibo,nreceived,MPI_REAL8,jdproc,jdproc
     &                  ,MPI_COMM_WORLD,req,ierr)
             call MPI_SEND(rHpaqsend,nsend,MPI_REAL8,jdproc,idproc
     &                  ,MPI_COMM_WORLD,ierr)
             call MPI_WAIT(req,status, ierr)
         
             do itotp=1,ntotproc(idproc)
                rHpaqrec(itotp)=rHpaqrec(itotp)+ recibo(itotp)
             enddo

          endif
      enddo  ! jdprocc
     
**end mpi communication
     
      return
      end subroutine difs
      
!--------------------------------------------------
!=======================================================
      end module mod_Hphi_01y2
