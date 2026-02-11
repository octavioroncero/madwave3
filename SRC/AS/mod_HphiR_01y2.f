      module mod_HphiR_01y2
!
!     Initialize magnitudes for the application of H phi
!                 at a FIXED R value
!     for a triatomic  01-2 system
!     in Jacobi coordinates in a body-fixed frame
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      
      implicit none

      save
      
      real*8, allocatable :: Easymp(:)
      integer :: nchan,nchan_dim
      integer, allocatable :: ielec_chan(:),nv_chan(:),j_chan(:)
      integer, allocatable :: iom_chan(:)
      real*8, allocatable :: rpaqproc(:)
     &         ,rflanz0proc(:)
     &         ,rHpaqsend(:),rHpaqrec(:),recibo(:)
     &     ,rHpaqproc(:,:)
     &     ,rpaq0(:)

      real*8, allocatable :: xl2p(:,:,:,:,:)
     &     , xj2(:,:,:,:)
      
      real*8, allocatable :: rr1m2(:,:),rr2m2(:,:)

      real*8, allocatable :: xpr1(:),xp2r1(:),xpr2(:),xp2r2(:)
     &               , pr1(:,:),p2r1(:,:)
     &     ,pr2(:,:),p2r2(:,:)
     &     ,r1d2(:,:),r1d1(:,:)  ,r1m2(:)
      integer, allocatable :: nopr1(:),nopr2(:)

      integer :: ir2_adiabatic
      real*8 :: r2_adiabatic
      integer,allocatable :: indtotRproc(:,:,:,:),ntotRproc(:)
      integer :: nbastotRproc
      integer, allocatable :: ii(:),jyj(:),vv(:),iiom_chan(:)
      integer, allocatable :: nord1(:)

      real*8, allocatable :: potRg_fixed(:,:,:,:)

      double precision,allocatable :: Hmat(:,:),Eigen(:),T(:,:)
      double precision,allocatable :: Hmattot(:,:)
      
      contains
!------------------------------------

      subroutine ordering_channels
      ! orders channels by asymptotic energies
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      implicit none
      include "mpif.h"

      integer :: ichan,i,j,j1,iv,ielec,iom
      real*8, allocatable :: assal(:)

      nchan_dim=nelecmax*(jmax-jini+1)*(nvmax-nvini+1)
     &          *(iommax-iommin+1)
      allocate(Easymp(nchan_dim),assal(nchan_dim)
     &     ,ielec_chan(nchan_dim),nv_chan(nchan_dim)
     &     ,j_chan(nchan_dim),iom_chan(nchan_dim)
     &     ,ii(nchan_dim),vv(nchan_dim)
     &     ,jyj(nchan_dim),iiom_chan(nchan_dim)
     &     ,potRg_fixed(npun1,nanguproc,nelecmax,nelecmax)
     &     ,nord1(nchan_dim))
      
! channels energy
      ichan=0
      do ielec=1,nelecmax
         do j=jini,jmax,inc
            do iom=iommax,iommin,-1
            if(iom.le.j)then
               do iv=nvini,noBCstates(j,ielec)-1
                  ichan=ichan+1
                  Easymp(ichan)=ediat(iv,j,ielec)
                  Assal(ichan)=ediat(iv,j,ielec)
                  ielec_chan(ichan)=ielec
                  nv_chan(ichan)=iv
                  j_chan(ichan)=j
                  iom_chan(ichan)=iom
                  ii(ichan)=ielec
                  jyj(ichan)=j
                  vv(ichan)=iv
                  iiom_chan(ichan)=iom
               enddo
            endif
            enddo
         enddo
      enddo
      nchan=ichan
!     ordering energies by

      do i=1,nchan
         nord1(i)=i
      enddo
               
 40   continue
         ichan=0
         do i=1,nchan-1
            J=NORD1(I)
            J1=NORD1(I+1)
            IF(Easymp(J).GT.Easymp(J1))THEN
               NORD1(I)=J1
               NORD1(I+1)=J
               iCHAN=iCHAN+1
            ENDIF
         enddo
      IF(iCHAN.GT.0)GO TO 40

!     assigned ordered channels

      
      do ichan=1,nchan
         j=nord1(ichan)
         ielec_chan(ichan)=ii(j)
         j_chan(ichan)=jyj(j)
         nv_chan(ichan)=vv(j)
         iom_chan(ichan)=iiom_chan(j)
         Easymp(ichan)=assal(j)
      enddo
         
!     printing channels

      OPEN(10,FILE='capture.dat',STATUS='unknown')
      write(6,*)' '
      write(6,*)' ordered channels'
      write(6,*)' ________________'
      write(6,*)nchan,'  ncanmax: i, v, j, ielec, Omg, Energy(eV) '
      write(10,*)nchan,' ncanmax: i, v, j, ielec , Omega, E(eV)'
      write(6,*)' _______________________________________________'
     
      do ichan=1,nchan
        
         write(10,'(5(1x,i5),2x,e15.7)')ichan
     &    ,nv_chan(ichan),j_chan(ichan),ielec_chan(ichan)
     &    ,iom_chan(ichan),Easymp(ichan)/conve1/8065.5d0

         write(6,'(5(1x,i5),2x,d15.7)')ichan
     &    ,nv_chan(ichan),j_chan(ichan),ielec_chan(ichan)
     &    ,iom_chan(ichan),Easymp(ichan)/conve1/8065.5d0

      enddo
      close(10)

!     deallocation

      deallocate(assal,ii,vv,jyj,nord1,iiom_chan)

      return
      end subroutine ordering_channels
      

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
      integer :: i,j,jp,k,ican0,jcan0
      real*8 :: xj,xom,r1,yyy,x1,x2

      real*8 :: xkinj(nangu,nangu)

      if(nangu.gt.1)then
         jjjmax=(nangu-1)*inc

         xj2(:,:,:,:)=0.d0

         do icanp=1,ncanproc(idproc)
            ican=ibasproc(icanp,idproc)
            ielec=nelebas(ican)
            iom=iombas(ican)
            iomdi=iomdiat(nelebas(ican))
            iomat=iomatom(nelebas(ican))

**>> Evaluation of the matrix for j^2

            do i=1,nojbas(ican)
            do j=1,nojbas(ican)
               xkinj(i,j)=0.d0
               if(i.eq.j)then
                  xj=dble(jbas(i,ican))
                  xom=dble(iomdi)
               xkinj(i,j)=(xj*(xj+1.d0)-xom*xom)*hbr*hbr*0.5d0
               endif
            enddo
            enddo

            do ir1=1,npun1
               r1=rmis1+dble(ir1-1)*ah1
               do k=1,nojbas(ican)
                  if(jbas(k,ican).le.jjjmax)then
                     yyy=xkinj(k,k)/(xm1red*r1*r1)
                     if(yyy.gt.rotcutmax)then
                         yyy=rotcutmax
                     endif
                     do i=1,nangu
                     do jp=1,nanguproc
                      j=indangreal(jp,idproc)
                      x1=Djmmp(i,jbas(k,ican),ielec,iom)
                      x2=Djmmp(j,jbas(k,ican),ielec,iom)
                      xj2(ir1,i,jp,icanp)=xj2(ir1,i,jp,icanp)+yyy*x1*x2
                     enddo
                     enddo
                  endif
                enddo
            enddo
!     adding projection up correction

            if(iom.gt.0)then
               ican0=-1
               do jcan0=1,ncan
                  if(iombas(jcan0).eq.0)ican0=jcan0
               enddo
               write(6,*)' ican0= ',ican0,iombas(ican0)
            
               do ir1=1,npun1
                  r1=rmis1+dble(ir1-1)*ah1
                  do k=1,nojbas(ican0)
                     if(jbas(k,ican0).le.jjjmax)then
                        yyy=rotcutmax
                        do i=1,nangu
                        do jp=1,nanguproc
                         j=indangreal(jp,idproc)
                         x1=Djmmp(i,jbas(k,ican0),ielec,iom)
                         x2=Djmmp(j,jbas(k,ican0),ielec,iom)
                      xj2(ir1,i,jp,icanp)=xj2(ir1,i,jp,icanp)+yyy*x1*x2
                        enddo
                        enddo
                     endif
                   enddo
               enddo
               do ir1=1,npun1
                  r1=rmis1+dble(ir1-1)*ah1
                  do k=1,nojbas(ican)
                     if(jbas(k,ican).le.jjjmax)then
                        yyy=rotcutmax
                        do i=1,nangu
                        do jp=1,nanguproc
                         j=indangreal(jp,idproc)
                         x1=Djmmp(i,jbas(k,ican),ielec,iom)
                         x2=Djmmp(j,jbas(k,ican),ielec,iom)
                       xj2(ir1,i,jp,icanp)=xj2(ir1,i,jp,icanp)-yyy*x1*x2
                        enddo
                        enddo
                     endif
                   enddo
               enddo
            end if !iom>0
            
         enddo                  ! icanp
      else
         xj2(:,:,:,:)=0.d0
      endif
         
      return
      end subroutine angkinj2
      
!--------------------------------------------------
!--------------------------------------------------
      subroutine  bfgridl2mat
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
      integer :: ielec,j,iom0,iom1,nnn,i,jj,jom,ir2
      real*8 :: yylmin,yylmax,yylmean,r2,r2eff,emineff,emeaneff
      real*8 :: eee,emaxeff,x11,x22,yyy,facmass,erotlmax
      integer :: iom,iterm,jjcanp,jjcan,jjelec,jjom,iq,iang,jang,jangp
      integer :: sigelec
      
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

      xl2p(:,:,:,:,:)=0.d0

*  starting loops

      do ielec=1,nelec
         sigelec=1
         if(sigdiat(ielec).eq.-1.d0.or.sigatom(ielec).eq.-1.d0)then
            sigelec=-1
         endif
         do j=j0,(nangu-1)*inc,inc
            iom0=iommin
            if(iparity.ne.sigelec*int((-1.d0)**Jtot).and
     &          .iommin.eq.0)iom0=1
            if(iom0.le.j)then

               iom1=min0(j,iommax)
               nnn=iom1-iom0+1
               write(6,*)' in l2: ',j,iom0,iom1,nnn
               if(nnn.gt.0)then
                  call l2mat(bfl2mat,facmass,iom0,iom1,j,Jtot
     &                         ,iommin,iommax)
              
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

                     do iom=iom0,iom1
                     do jom=iom0,iom1
                        xl2mat(iom,jom)=0.d0
                        if(dabs(bfl2mat(iom,jom)).gt.1.d-8
     &                       .and.r2eff.gt.1.d-3)then
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
                              xl2p(ir2,iang,jangp,jjcanp,iq)=
     &                     xl2p(ir2,iang,jangp,jjcanp,iq)+x11*x22*yyy
                           enddo
                           enddo
                           endif

                        endif
                     enddo  ! ican1p

                  enddo  ! jom
                  enddo  ! iom

               enddo ! ir2
            endif ! nnn > 0

            endif  ! iom0 < j
         enddo  ! j=0,jmax
      enddo  ! ielec=1,nelec

!     deallocating auxiliary matrices
      deallocate(bfl2mat
     &       , xl2mat
     &       , hmat,T,eigval
     &       , ind)
      
      return
      end subroutine  bfgridl2mat
            
!--------------------------------------------------
      subroutine HphiR
*********************************************************************
**                  Subroutine HphiR                               **
**                                                                 **
**                Evaluates  H Phi(t)                              **
**            using Reactant Jacobi coordinates                    **
**         at frozen R2=Rg--> ir2_adiabatic,    r2_adiabatic       **
**                                                                 **
*********************************************************************

      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      implicit none
      include "mpif.h"
      
      integer :: ierror
      
      integer*8 plan,flags
      INTEGER :: req, status(MPI_STATUS_SIZE),ierr
      
      real*8 :: xp2r2p2r1(npun2,npun1),rphi0(npun2,npun1)
      real*8 :: auxR(npun1*npun2),raux(npun2,npun1)
      integer :: icanp,ican,ielec,iom,jangp,jang,ir1,ir2,ir
      integer :: i,nn2,nn1,iii,i1,i2,ican2p,ican2,ielec2,iom2
      integer :: itotp,jdproc,iangp,iang,jdprocc,nsend,nreceived
      integer :: semaforo
      integer*8 :: max_semaforo,icount,i0
      real*8 :: divi,ekin
      integer :: iang_proc,ican_proc,ir1p,itot

      rHpaqproc(:,:)=0.d0
      rHpaqrec(:)=0.d0

      
******************
**  1) T Phi(t) **
******************

      do iang_proc=1,nanguproc
         iang=indangreal(iang_proc,idproc)
         do ican_proc=1,ncanproc(idproc)
            ican=ibasproc(ican_proc,idproc)

            rphi0(:,:)=0.d0
            raux(:,:)=0.d0
            
            do ir1=1,npun1
               i= indtotRproc(ir1,ican_proc,iang_proc,idproc)
               rphi0(ir2_adiabatic,ir1)=rpaqproc(i)
            enddo               ! ir1
            do ir1=1,npun1
                do ir1p=1,npun1
                   raux(ir2_adiabatic,ir1)=raux(ir2_adiabatic,ir1)
     &                 +r1d2(ir1,ir1p)*rphi0(ir2_adiabatic,ir1p)
                enddo
            enddo
            do ir1=1,npun1
               i= indtotRproc(ir1,ican_proc,iang_proc,idproc)
               rHpaqrec(i)=rHpaqrec(i)+raux(ir2_adiabatic,ir1)
            enddo               ! ir1

         enddo               ! ican
      enddo ! iang_proc
             
******************
**  2) V Phi(t) **
******************
      do icanp=1,ncanproc(idproc)
         ican=ibasproc(icanp,idproc)
         ielec=nelebas(ican)
         iom=iombas(ican)
         do jangp=1,nanguproc
            jang=indangreal(jangp,idproc)


            do ican2p=1,ncanproc(idproc)
                ican2=ibasproc(ican2p,idproc)
                ielec2=nelebas(ican2)
                iom2=iombas(ican2)
                if(iom.eq.iom2)then
                   do ir1=1,npun1 
                      itot=indtotRproc(ir1,icanp,jangp,idproc)
                      itotp=indtotRproc(ir1,ican2p,jangp,idproc)
                      rHpaqrec(itotp)=rHpaqrec(itotp)
     &              + rpaqproc(itot)*potRg_fixed(ir1,jangp,ielec,ielec2)
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
                     do ir1=1,npun1
                        itot=indtotRproc(ir1,icanp,jangp,idproc)
                        itotp=indtotRproc(ir1,ican2p,iangp,jdproc)
                     rHpaqproc(itotp,jdprocc)=rHpaqproc(itotp,jdprocc)
     &                   +rpaqproc(itot)*(xj2(ir1,iang,jangp,icanp))
!     &                   +rpaqproc(itot)*(xj2(ir1,iang,jangp,icanp)
!     &                       +xl2p(ir2_adiabatic,iang,jangp,icanp,0))
                     enddo
!                  elseif(iom.eq.iom2+1.and.ielec.eq.ielec2)then
!                     do ir1=1,npun1
!                        itot=indtotRproc(ir1,icanp,jangp,idproc)
!                        itotp=indtotRproc(ir1,ican2p,iangp,jdproc)
!                    rHpaqproc(itotp,jdprocc)=rHpaqproc(itotp,jdprocc)
!     &          +rpaqproc(itot)*xl2p(ir2_adiabatic,iang,jangp,icanp,1)
!                        
!                     enddo
!                  elseif(iom.eq.iom2-1.and.ielec.eq.ielec2)then
!                     do ir1=1,npun1
!                        itot=indtotRproc(ir1,icanp,jangp,idproc)
!                        itotp=indtotRproc(ir1,ican2p,iangp,jdproc)
!                       rHpaqproc(itotp,jdprocc)=rHpaqproc(itotp,jdprocc)
!     &          +rpaqproc(itot)*xl2p(ir2_adiabatic,iang,jangp,icanp,-1)
!                       
!                     enddo
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

            do itotp=1,ntotRproc(idproc)
              rHpaqrec(itotp)=rHpaqrec(itotp)+rHpaqproc(itotp,jdprocc)
            enddo

         else

             nsend=ntotRproc(jdproc)
             nreceived=ntotRproc(idproc)
             do itotp=1,ntotRproc(jdproc)
                rHpaqsend(itotp)=rHpaqproc(itotp,jdprocc)
             enddo
          
             call MPI_IRECV(recibo,nreceived,MPI_REAL8,jdproc,jdproc
     &                  ,MPI_COMM_WORLD,req,ierr)
             call MPI_SEND(rHpaqsend,nsend,MPI_REAL8,jdproc,idproc
     &                  ,MPI_COMM_WORLD,ierr)
             call MPI_WAIT(req,status, ierr)
         
             do itotp=1,ntotRproc(idproc)
                rHpaqrec(itotp)=rHpaqrec(itotp)+ recibo(itotp)
             enddo

          endif
      enddo  ! jdprocc
     
**end mpi communication
     
      return
      end subroutine HphiR
!-------------------------------------------------------------
      subroutine Tradial1
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
!   only for small r, H at fixed R_big              !

      implicit none
      include "mpif.h"
      integer :: ierr,ip,iang_proc,iang,ican_proc,ican,ir1,imem,iii

!     indexes for only one Rg value
      
      allocate(indtotRproc(npun1,ncanprocdim,nanguprocdim,0:nproc-1)
     &    ,    ntotRproc(0:nproc-1))

      indtotRproc(:,:,:,:)=0.d0
      nbastotRproc=0
      do ip=0,nproc-1
         ntotRproc(ip)=nanguproc*ncanproc(ip)*npun1
         write(6,*)' nbastotRproc(',ip,')= ',ntotRproc(ip)
         nbastotRproc=max0(nbastotRproc,ntotRproc(ip))
         iii=0
         do iang_proc=1,nanguproc
            iang=indangreal(iang_proc,ip)
            do ican_proc=1,ncanproc(ip)
               ican=ibasproc(ican_proc,ip)
               do ir1=1,npun1
                  iii=iii+1
                  indtotRproc(ir1,ican_proc,iang_proc,ip)=iii
               enddo
            enddo
         enddo
      enddo

! paqs and Hpaqs

      allocate(rpaqproc(nbastotRproc),recibo(nbastotRproc)
     &         ,rflanz0proc(nbastotRproc)
     &         ,rpaq0(nbastotRproc)
     &         ,rHpaqsend(nbastotRproc),rHpaqrec(nbastotRproc)
     &         ,rHpaqproc(nbastotRproc,ncouprocmax)
     &         ,stat=ierr)

! centrifugal terms

      allocate(xl2p(npun2,nangu,nanguprocdim,ncanprocdim,-1:1)
     &       , xj2(npun1,nangu,nanguprocdim,ncanmax)
     &       , stat=ierr)
      if(ierr.ne.0)then
        write(*,*)" error in initmem for centrifugal "
        stop
      else
         norealproc_mem=norealproc_mem+imem
         imem=
     &    npun2*nangu*nanguprocdim*ncanprocdim*3
     &    +npun1*nangu*nanguprocdim*ncanmax
          write(6,*)' xl2p,xj2 in set_vectors_Hphi, proc= '
     &       ,idproc
     &       ,' imem=',imem
     &       ,' memory(Gb)= '
     &       ,dble(imem*8)*1.d-9
         call flush(6)
      endif
      
**>> Angular kinetic terms 

      call  angkinj2
      call bfgridl2mat
* Kinetic term for r1=rp
* for radial derivatives
* represented in a grid, using sinc functions to start

         write(6,*)' Using sinc DVR for radial kinetic terms '
         write(6,*)

         allocate(r1d2(npun1,npun1)  ,r1d1(npun1,npun1)  ,r1m2(npun1)
     &       , stat=ierr)

        if(npun1.gt.1)then
            call sincbasnew(rmis1,rfin1,npun1,xm1red,hbr
     &                   ,r1d2,r1d1,r1m2)  ! second derivative, first derivative and r^{-2}
         elseif(npun1.eq.1)then
            r1m2(1)=0.5d0*hbr*hbr/(xm1red*rmis1*rmis1)
         endif
      return
      end subroutine Tradial1

!------------------------------------------------------

      subroutine set_rpaq_funchan(ichan)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      implicit none
      include "mpif.h"
      
      integer :: ichan,ican_proc,ican,ielec,iom,iang_proc,iang,ir1,itot
      
       rpaqproc(:)=0.d0
       do ican_proc=1,ncanproc(idproc)
          ican=ibasproc(ican_proc,idproc)
          ielec=nelebas(ican)
          iom=iombas(ican)
          if(ielec.eq.ielec_chan(ichan)
     &             .and.iom.eq. iom_chan(ichan))then
               
              do iang_proc=1,nanguproc
                 iang=indangreal(iang_proc,idproc)
                 do ir1=1,npun1
                    itot= indtotRproc(ir1,ican_proc,iang_proc,idproc)
                    rpaqproc(itot)=rpaqproc(itot)+
     &                        Djmmp(iang,j_chan(ichan),ielec,iom)
     &                    *fd(ir1,nv_chan(ichan),j_chan(ichan),ielec)
                 enddo
              enddo
            endif
         enddo               ! icanproc
         return
         end subroutine set_rpaq_funchan
!--------------------------------------------------
      subroutine V_Rg_fixed
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      implicit none
      include "mpif.h"
      
      integer :: ielec,jelec,iang_proc,iang,ir,ir2,ir1

      do ielec=1,nelecmax
         do jelec=1,nelecmax
            if(ielec==jelec)then
               potRg_fixed(:,:,ielec,jelec)=Vcutmax
            else
               potRg_fixed(:,:,ielec,jelec)=0.d0
            end if
         end do
      end do

      do ielec=1,nelecmax
         do jelec=1,nelecmax
            do iang_proc=1,nanguproc
               iang=indangreal(iang_proc,idproc)
               do ir=1, npunreal(iang)
                  ir2=indr2(ir,iang)
                  if( ir2 == ir2_adiabatic)then
                     ir1=indr1(ir,iang)
                     potRg_fixed(ir1,iang_proc,ielec,jelec)
     &                    =vvv(ir,iang_proc,ielec,jelec)
                  end if
               end do
            end do
         end do
      end do
               
      return
      end subroutine V_Rg_fixed
!--------------------------------------------------

!=======================================================
      end module mod_HphiR_01y2
