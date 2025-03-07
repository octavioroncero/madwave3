      module mod_coortrans01_02
!
!     to transform wvp to products and 02 products to intermediate coordinates
!     to calculate individual fluxes on each v,j of 02 products
!         v. denotes vibrational and electronic states
!
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      use mod_absorcion_01y2
      use mod_photoini_01y2
      use mod_colini_01y2
      use mod_flux_01y2

      implicit none

      save
      real*8,allocatable :: S2prodfac(:,:,:) 
     &     ,Cvjprod(:,:,:),Cvjprodp(:,:,:)
      real*8, allocatable :: transRgjac(:),prodwf(:,:,:,:,:,:)
      real*8, allocatable :: dtrans(:,:,:,:)
      complex*16,allocatable :: zS2prod(:,:,:,:)

      contains
********************************************
*   functions of   mod_flux_01y2.f    *
********************************************
!=======================================================

!--------------------------------------------------
      subroutine ini_transcoor

      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      use mod_absorcion_01y2
      use mod_photoini_01y2
      use mod_colini_01y2
      use mod_flux_01y2
      implicit none
      integer :: ierror,ie,iv,j,indt,it0,iloop0
      integer*8 ::imem
      real*8 :: ekinini,pini,paqini,ekinfin,pfin
      
      write(6,'(40("="),/,10x,"mod_coortrans01_02",/,40("="))')
      write(6,*)'    Initializing variables for Smatprod '
      call flush(6)
      imem=      netot*(nvmaxprod-nviniprod+1)*(jmaxprod-jiniprod+1)
     &           +(iommaxprod-iomminprod+1)*(nvmaxprod-nviniprod+1)
     &           *(jmaxprod-jiniprod+1)*(2+netot*2)
     &    +nbastot/nproc
     &    +(n2prod1-n2prod0+1)*(nangproj1-nangproj0+1)*nelecmax
     &          *(nvmaxprod-nviniprod+1)*(jmaxprod-jiniprod+1)
     &          *(iommaxprod-iomminprod+1)
     &    +(n2prod1-n2prod0+1)*(nangproj1-nangproj0+1)
     &             *(iommaxprod-iomminprod+1)*ncanmax

      write(6,*)' memory to be allocated in ini_transcoor '
     &    ,'  in proc= ',idproc,' is= '
     &    ,8.d0*(dble(imem)*1.d-9),' Gb'
      call flush(6)
      
      
         allocate(S2prodfac(netot,nviniprod:nvmaxprod,jiniprod:jmaxprod)
     &  ,Cvjprod(nviniprod:nvmaxprod,jiniprod:jmaxprod
     &                                         ,iomminprod:iommaxprod)
     &  ,Cvjprodp(nviniprod:nvmaxprod,jiniprod:jmaxprod
     &                                         ,iomminprod:iommaxprod)
     &  ,zS2prod(netot,nviniprod:nvmaxprod,jiniprod:jmaxprod
     &                                        ,iomminprod:iommaxprod)
     &       , stat=ierror)
         if(ierror.ne.0)then
           write(*,*)" error in initmem for paraS2prod "
           call flush(6)
           stop
         else
            imem=
     &           netot*(nvmaxprod-nviniprod+1)*(jmaxprod-jiniprod+1)
     &           +(iommaxprod-iomminprod+1)*(nvmaxprod-nviniprod+1)
     &           *(jmaxprod-jiniprod+1)*(2+netot*2)
            norealproc_mem=norealproc_mem+imem
            write(6,*)' Cvjprod,S2prod in ini_transcoor, proc= ',idproc
     &       ,' imem=',imem
     &       ,' memory(Gb)= '
     &       ,dble(imem*8)*1.d-9
            call flush(6)
         endif
         S2prodfac(:,:,:)=0.d0
         Cvjprod(:,:,:)=0.d0
         Cvjprodp(:,:,:)=0.d0
         zS2prod(:,:,:,:)=dcmplx(0.d0,0.d0)

!-------------------------------------------------------
!  S2 matrix factor for collisions and photodissociation
!-------------------------------------------------------
! iphoto= 0  --> state-to-state in collisions
      if(iphoto.eq.0)then
            
         do ie=1,netot
            ekinini=etotS2(ie)
            pini=0.d0
            if(ekinini.gt.0.d0)pini=dsqrt(ekinini*xmasa0*2.d0/hbr/hbr)
            paqini=apaqini(ie)
            do iv=nviniprod,nvmaxprod
               do j=jiniprod,jmaxprod
                  ekinfin=etotS2(ie)-ediatprod(iv,j)
                  if(ekinfin.gt.0.d0)then
                     pfin=dsqrt(ekinfin*xm2prod*2.d0/hbr/hbr)
                  else
                     pfin=0.d0
                  endif
                  if(paqini.lt.1.d-10)then
                     S2prodfac(ie,iv,j)=0.d0
                  else
                     S2prodfac(ie,iv,j)=
     &                    hbr*hbr*pfin*pini/(2.d0*paqini)/xmasa0/xm2prod
                  endif
               enddo
            enddo
         enddo

      else
! iphoto > 0  --> state-to-state in photodissociation

         do ie=1,netot
            do iv=nviniprod,nvmaxprod
               do j=jiniprod,jmaxprod
 
                  ekinfin=etotS2(ie)-ediatprod(iv,j)
                  if(ekinfin.gt.0.d0)then
                     pfin=dsqrt(ekinfin*xm2prod*2.d0/hbr/hbr)
                  else
                     pfin=0.d0
                  endif
!pru                  S2prodfac(ie,iv,j)=2.d0*pi*pfin/xm2prod
                  S2prodfac(ie,iv,j)=4.d0*pfin/(pi*xm2prod)
               enddo
            enddo
         enddo
      endif
!---------------------------------------------------------
!  Transformation matrices of final states of 02 products
!---------------------------------------------------------

      write(6,*)
      write(6,*)' Preparing matrices for final states of 02 products'
      write(6,*)

      call changeini

!----------------------------------------------
!reading information for continuing propagation
!----------------------------------------------      
      it0=0
      indt=0
      iloop0=0
      open(3,file='cont.data',status='old',err=1)
      read(3,*)it0,indt,iloop0
      close(3)
 1    continue
      if(it0.eq.0)then
            zS2prod(:,:,:,:)=dcmplx(0.d0,0.d0)
      else
         if(ncontfile.eq.1)then
            write(name,'("cont",i2.2,".S2prod")')idproc
            open(4,file=name,status='old',form='unformatted')
            read(4)zS2prod
            write(6,*)' zS2prod read in ',name
            close(4)
         endif
      endif

      write(6,*)' ending ini_transcoor in proc= ',idproc
     &     ,' norealproc_mem=',norealproc_mem
     &     ,' nointegerproc_mem=',nointegerproc_mem
     &     ,' memory(Gb)= '
     &   ,dble(nointegerproc_mem*4+norealproc_mem*8)*1.d-9
      call flush(6)
       
      return     
      end subroutine ini_transcoor

!--------------------------------------------------

      subroutine changeini
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      use mod_absorcion_01y2
      use mod_photoini_01y2
      use mod_colini_01y2
      use mod_flux_01y2
      implicit none
      integer :: lmax,jmax2,ierror,icanp,ican,ielec,iom,iang,jangp,jang
      integer :: ir2,ir1,ir,itot,jinip,j,iv,m,irga
      real*8 :: box1,box2,Rgprod,xt1,xt2,xt3,xt4,A,B,cg,r2,xxx1,xxx2,xxx
      real*8 :: r1inter,r1,xk,xjacobian1,xjacobian,Rgb,xn,ang
      real*8 :: cgama,sqr,rpbmin,rpbmax,rga,rpa,rpb
      real*8 :: x1p(npun1),x1p2(npun1),x2p(npun2),x2p2(npun2)
      real*8 :: prodaux(n2prod0:n2prod1,nangproj0:nangproj1,nelecmax)
      real*8 :: pinfun(jiniprod:jmaxprod),raux(npun2,npun1)
      integer*8 :: imem

      include "mpif.h"
*
* Preparing matrices for coordinate transformation:
*           a) from Reactant Jacobi --> Product Jacobi 
*                             if    iprod.gt.1
* The transformation is done change one variable in each step
*              to involve the minimum quantity of nested loops
*              and the first one is towads the Jacobi Rg_infty
*              to reduce to a single point this variable.
*

      box1=dble(npun1-1)*ah1
      box2=dble(npun2-1)*ah2
      call sinmom(box1,npun1,npun1,xm1red,hbr,x1p,x1p2)
      call sinmom(box2,npun2,npun2,xm2red,hbr,x2p,x2p2)
      lmax=(nangu-1)*inc
      jmax2=lmax*2

      Rgprod=Rbalinprod

      xt1=0.d0
      xt2=0.d0
      xt3=0.d0
      xt4=0.d0
      allocate(transRgjac(nbastot/nproc)
     &      ,prodwf(n2prod0:n2prod1,nangproj0:nangproj1,nelecmax
     &          ,nviniprod:nvmaxprod,jiniprod:jmaxprod
     &           ,iomminprod:iommaxprod)
     &       ,dtrans(n2prod0:n2prod1,nangproj0:nangproj1
     &             ,iomminprod:iommaxprod,ncanmax)
     &     ,stat=ierror)

      imem=
     &     nbastot/nproc
     &    +(n2prod1-n2prod0+1)*(nangproj1-nangproj0+1)*nelecmax
     &          *(nvmaxprod-nviniprod+1)*(jmaxprod-jiniprod+1)
     &          *(iommaxprod-iomminprod+1)
     &    +(n2prod1-n2prod0+1)*(nangproj1-nangproj0+1)
     &             *(iommaxprod-iomminprod+1)*ncanmax
      if(imem.lt.0)then
          write(6,*)' Be carefull: imem = ',imem,' < 0' 
          write(6,*)' nbastot/nproc= ',nbastot/nproc
          write(6,*)' second term = ',
     &     (n2prod1-n2prod0+1)*(nangproj1-nangproj0+1)*nelecmax
     &          *(nvmaxprod-nviniprod+1)*(jmaxprod-jiniprod+1)
     &          *(iommaxprod-iomminprod+1)
          write(6,*)' third  term = ',
     &     (n2prod1-n2prod0+1)*(nangproj1-nangproj0+1)
     &             *(iommaxprod-iomminprod+1)*ncanmax
      endif


      norealproc_mem=norealproc_mem+imem
      write(6,*)' transRgjac,prodwf,dtrans in changeini, proc= '
     &       ,idproc
     &       ,' imem=',imem
     &       ,' memory(Gb)= '
     &       ,dble(imem*8)*1.d-9
      call flush(6)
      
      if(ierror.ne.0)then
         write(*,*)" error in changeini for mod_coortrans "
         call flush(6)
         stop
      endif
      transRgjac(:)=0.d0
      prodwf(:,:,:,:,:,:)=0.d0
      dtrans(:,:,:,:)=0.d0
********************************************************************
***> matrices for Reactant Jacobi --> product Jacobi transformation
********************************************************************

      B=(xm0+xm1)/xm1 - xm2/(xm0+xm2)
      A=(xm0+xm2)/xm2 - xm1/(xm0+xm1)

*1a) preparation to transform the wvp to R_prod , R_reac, Gama_reac
*             to project in product channel

      do icanp=1,ncanproc(idproc)
         ican=ibasproc(icanp,idproc)
         ielec=nelebas(ican)
         iom=iombas(ican)
         do iang=nangproj0,nangproj1
         do jangp=1,nanguproc
            jang=indangreal(jangp,idproc)
            if(jang.eq.iang)then      
               cg=cgamma(iang)

               do ir2=1,npun2
               do ir1=1,npun1
                  raux(ir2,ir1)=0.d0
               enddo
               enddo
               
               do ir2=1,npun2
                  r2=rmis2+dble(ir2-1)*ah2
                  xxx1=(xm0+xm2)*Rgprod/xm2
                  xxx1=xxx1*xxx1
                  xxx2=r2*r2*(cg*cg-1.d0)
                  if(xxx1+xxx2.gt.0.d0)then
                     xxx=dsqrt(xxx1+xxx2)
                     r1inter=r2*cg+xxx
                     r1inter=r1inter/A
                     if(r1inter.gt.rmis1.and.r1inter.lt.rfin1)then
                         do ir1=1,npun1
                            r1=rmis1+dble(ir1-1)*ah1
                            xxx=0.d0
                            do m=1,npun1
                               xk=x1p(m)
                               xxx=xxx+dsin(xk*r1)*dsin(xk*r1inter)
                            enddo
                            raux(ir2,ir1)=xxx*2.d0/dble((npun1+1))
                            xjacobian1=dsqrt(r2*r2*(cg*cg-1.d0)
     &                                 +(Rgprod*(xm0+xm2)/xm2)**2)
                                 
                            xjacobian=Rgprod*(((xm0+xm2)/xm2)**2)
                            xjacobian=dsqrt(xjacobian/(A*xjacobian1))
                            raux(ir2,ir1)=raux(ir2,ir1)*xjacobian
                            raux(ir2,ir1)=raux(ir2,ir1)/dsqrt(ah1)
                         enddo
                     endif
                  endif
               enddo

               do ir=1,npunreal(iang)
                  ir1=indr1(ir,iang)
                  ir2=indr2(ir,iang)
                  itot=indtotproc(ir,icanp,jangp,idproc)

                
                  if(itot.lt.1.or.itot.gt.nbastot/nproc)then
                  write(6,*)' transRg: itot= ',itot,icanp,jangp,idproc
                     call flush(6)
                     call MPI_BARRIER(MPI_COMM_WORLD, ierror)
                     stop
                  endif
                  transRgjac(itot)=raux(ir2,ir1)
                  xt1=xt1+transRgjac(itot)*transRgjac(itot)
               enddo
            endif  ! jang=iang
         enddo  ! jangp
         enddo  !  iang
      enddo !  icanp 

      if(idproc.eq.0)then
        write(6,*)' xt1= ',xt1
***> Product asymptotic wv in R_prod, R_reac, Gam_reac

         write(6,*) 
         write(6,*)"--> Norm of product wf's for projection"
         write(6,*)
      endif

      write(6,*)'    rp of products sampled in transformation between'
      write(6,*)'        n2prod0,n2prod1 = ',n2prod0,n2prod1
         rpbmin=1.d10
         rpbmax=0.d0
         Rgb=Rbalinprod
         do iRga=n2prod0,n2prod1
            Rga=rmis2+dble(iRga-1)*ah2
            do iang=nangproj0,nangproj1
               cgama=cgamma(iang)
               xxx=0.d0
               sqr=Rga*Rga*(cgama*cgama-1.d0)+(Rgb*(xm0+xm2)/xm2)**2
               if(sqr.gt.0.d0)then
                  rpa=Rga*cgama+dsqrt(sqr)
                  rpa=rpa/A
               
                  rpb=(rpa*xm1/(xm0+xm1))**2+Rga*Rga
     &            +2.d0*rpa*Rga*cgama*xm1/(xm0+xm1)
                  rpb=dsqrt(rpb)

                  rpbmin=dmin1(rpbmin,rpb)
                  rpbmax=dmax1(rpbmax,rpb)
               endif
            enddo
         enddo
         write(6,*)'  rpbmin,rpbmax= ',rpbmin,rpbmax
         write(6,*)
      call flush(6)
         
      do iom=iomminprod,iommaxprod
         jinip=max0(jiniprod,iom)
         do j=jinip,jmaxprod
            do iv=nviniprod,nvmaxprod
               Rgb=Rbalinprod
               call prodwftrans(Rgb,prodaux,iv,j,iom)
*norm check
               xn=0.d0
               do ielec=1,nelecmax
                  do iang=nangproj0,nangproj1
                     do ir2=n2prod0,n2prod1
                        xxx=prodaux(ir2,iang,ielec)
                        prodwf(ir2,iang,ielec,iv,j,iom)=
     &                     prodwf(ir2,iang,ielec,iv,j,iom)+xxx             
                        xn=xn+xxx*xxx
                     enddo
                  enddo
               enddo
               if(idproc.eq.0)write(6,"(
     &             ' norm (iom,j,v=',i3,',',i3,',',i3
     &                   ,' )= ',e15.7)")iom,j,iv,xn

            enddo ! iv=nviniprod,nvmaxprod
         enddo ! j=jinip,jmaxprod
               
               
         if(idproc.eq.0.and.iwrt_pot.eq.1)then
            do iv=nviniprod,nvmaxprod
            do ielec=1,nelecmax
               write(name,'("prodwv.v",i2.2,".e",i2.2,".Omg",i2.2)')
     &                  iv,ielec,iom
               open(43,file=name,status='unknown')
               do iang=nangproj0,nangproj1
                  ang=dacos(cgamma(iang))*180.d0/pi
                  do ir2=n2prod0,n2prod1
                     r2=rmis2+dble(ir2-1)*ah2
                     do j=jiniprod,jmaxprod
                        pinfun(j)=0.d0
                           pinfun(j)=pinfun(j)
     &                          +prodwf(ir2,iang,ielec,iv,j,iom)**2
                        if(dabs(pinfun(j)).lt.1d-20)pinfun(j)=0.d0
                     enddo
                
                     write(43,'(100(1x,e15.7))')r2,ang
     &                         ,(pinfun(j),j=jiniprod,jmaxprod)
                  enddo
                  write(43,'()')
               enddo             
               close(43)            
            enddo ! ielec
            enddo ! iv
         endif
      enddo ! iom  

      write(6,*)' --> prod functions transformed to intermediate coor.'
      call flush(6)
***> calculating bf reac-2-prod transformation

      call bfreac2prod
      
      write(6,*)' --> end of react-2-prod body-fixed transf. '
       
!      if(idproc.eq.0)
      write(6,*)' xt1,xt2,xt3,xt4= ',xt1,xt2,xt3,xt4

     
      call flush(6)
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)

      return
      end  subroutine changeini
!--------------------------------------------------
      subroutine prodwftrans(Rgb,prodaux,iv,j,iom)
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      use mod_absorcion_01y2
      use mod_photoini_01y2
      use mod_colini_01y2
      use mod_flux_01y2
      implicit none
      real*8 :: prodaux(n2prod0:n2prod1,nangproj0:nangproj1,nelecmax)
      real*8 :: ff(npun2,2),xx(npun2)
      integer :: iv,j,iom,ielec,ir2,lmax,jmax2,iRga,iang,iold,m2,mp2
      integer :: j2
      real*8 :: Rgb,facnorm,A,B,cgama,xxx,sqr,rpa,rpb,fdiatom,angjac,xn
      real*8 :: angpart,xjac1,xjaco,Rga,xn0

      include "mpif.h"
      
* Generates the prod wavefunctions in R_reac,Gam_reac coordinates
      prodaux(:,:,:)=0.d0

**> initialization of variables for interpolation
      xn0=0.d0
      do ielec=1,nelecmax
         
         do ir2=1,n2prod1
            xx(ir2)=rmis2+dble(ir2-1)*ah2
            ff(ir2,1)=fdprod(ir2,iv,j,ielec)
        
            xn0=xn0+ ff(ir2,1)* ff(ir2,1)
            ff(ir2,2)=0.d0
         enddo
         call splset(ff,xx,n2prod1,npun2)

         lmax=(nangu-1)*inc
         jmax2=lmax*2

         facnorm=1.d0
         if(inc.eq.2)facnorm=1.d0/dsqrt(2.d0)

         A=(xm0+xm2)/xm2-xm1/(xm0+xm1)
         B=(xm0+xm1)/xm1-xm2/(xm0+xm2)

*---  > (ra,Ra,gama --> Rb Ra,gama)
         do iRga=n2prod0,n2prod1
            Rga=rmis2+dble(iRga-1)*ah2
            do iang=nangproj0,nangproj1
               cgama=cgamma(iang)
               xxx=0.d0
               sqr=Rga*Rga*(cgama*cgama-1.d0)+(Rgb*(xm0+xm2)/xm2)**2
               if(sqr.gt.0.d0)then
                  rpa=Rga*cgama+dsqrt(sqr)
                  rpa=rpa/A
               
                  rpb=(rpa*xm1/(xm0+xm1))**2+Rga*Rga
     &            +2.d0*rpa*Rga*cgama*xm1/(xm0+xm1)
                  rpb=dsqrt(rpb)

                  iold=2
                  fdiatom=0.d0
                  if(rpb.ge.rmis2.and.rpb.le.rfin2)then
                     call splinqq(ff,xx,iold,n2prod1,rpb,npun2,fdiatom)
                  endif

                  angjac=Rgb*Rgb+B*B*rpb*rpb-(Rga*(xm0+xm1)/xm1)**2
                  angjac=angjac*0.5d0/(B*rpb*Rgb)

                  if(dabs(angjac).le.1.d0)then
                     angjac=dacos(angjac)

                     m2=2*(iom-iomatom(ielec))
                     mp2=2*iomdiat(ielec)
                     j2=2*j
                     call dwigner(pm,jmax2,m2,mp2,angjac,ndim)
                     xn=dsqrt((2.d0*dble(j)+1.d0)*0.5d0)
                     angpart=pm(j2)*xn*dsqrt(weight(iang))*facnorm

                  xjac1=(Rga**2)*(cgama**2-1.d0)+ (Rgb*(xm0+xm2)/xm2)**2
                     xjaco=Rgb*( ((xm0+xm2)/xm2)**2)/(A*dsqrt(xjac1))
                 
                     xxx=fdiatom*angpart*Rga*rpa*dsqrt(xjaco)/(rpb*Rgb)
                  endif
               endif
            
               prodaux(iRga,iang,ielec)=xxx

            enddo
         enddo
      enddo                     ! ielec

!      xn=0.d0
!      do ielec=1,nelecmax
!         do iRga=n2prod0,n2prod1
!            do iang=nangproj0,nangproj1
!               xn=xn+ prodaux(iRga,iang,ielec)**2
!            enddo
!         enddo
!      enddo
!      write(6,*)xn0,xn
!      call flush(6)
      
      return
      end subroutine prodwftrans

!--------------------------------------------------

      subroutine bfreac2prod
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      use mod_absorcion_01y2
      use mod_photoini_01y2
      use mod_colini_01y2
      use mod_flux_01y2
      implicit none
      real*8 :: Rgb,A,B,cgama,Rga,xxx,sqr,rpa,rpb,angjac,cgamb,angb
      real*8 :: deno,ctrans,angtrans,x1,ss
      integer :: iang,ir,ican,iomb,ifail,iRga,iomreac,ielec,iomprod,m2
      integer :: jtot2,mp2,ierror

      include "mpif.h"

      Rgb=Rbalinprod
      A=(xm0+xm2)/xm2-xm1/(xm0+xm1)
      B=(xm0+xm1)/xm1-xm2/(xm0+xm2)
      if(iommaxprod.gt.Jtot)then
         write(6,*)' iommaxprod= ',iommaxprod
     &        ,'  > Jtot= ',Jtot
         call flush(6)
         call MPI_BARRIER(MPI_COMM_WORLD, ierror)
         stop
      endif

      dtrans(:,:,:,:)=0.d0
      do iang=nangproj0,nangproj1
         cgama=cgamma(iang)
         do ir=n2prod0,n2prod1

            if(Jtot.eq.0)then

               do ican=1,ncan
               do iomb=iomminprod,iommaxprod
                  dtrans(ir,iang,iomb,ican)=1.d0
               enddo
               enddo
            elseif(Jtot.gt.0)then
               ifail=0
               iRga=ir
               Rga=rmis2+dble(iRga-1)*ah2
               xxx=0.d0
               sqr=Rga*Rga*(cgama*cgama-1.d0)
     &             +(Rgb*(xm0+xm2)/xm2)**2
               if(sqr.lt.0.d0)ifail=1
                   rpa=Rga*cgama+dsqrt(sqr)
                   rpa=rpa/A
               if(ifail.eq.0)then
                  rpb=(rpa*xm1/(xm0+xm1))**2+Rga*Rga
     &               +2.d0*rpa*Rga*cgama*xm1/(xm0+xm1)
                  rpb=dsqrt(rpb)

                  angjac=Rgb*Rgb+B*B*rpb*rpb-(Rga*(xm0+xm1)/xm1)**2
                  angjac=angjac*0.5d0/(B*rpb*Rgb)

                  cgamb=angjac
                  if(dabs(angjac).le.1.d0)then
                     angb=dacos(angjac)

                     deno=(B*rpb)**2 +Rgb*Rgb 
     &                   - 2.d0*B*rpb*Rgb*dcos(angb)

                     ctrans=(B*rpb*dcos(angb) -Rgb)/dsqrt(deno)
                     if(dabs(ctrans).le.1.d0)then
                        angtrans=dacos(ctrans)
                        do ican=1,ncan
                           iomreac=iombas(ican)
                           ielec=nelebas(ican)
                           do iomprod=iomminprod,iommaxprod
                             m2=iomreac*2
                             jtot2=jtot*2
                             mp2=iomprod*2
                             call dwigner(pm,jtot2,m2,mp2,angtrans,ndim)
                             x1=pm(jtot2)
                             if(iomreac.ne.0)then
                            call dwigner(pm,jtot2,-m2,mp2,angtrans,ndim)
                                ss = dble(iparity)
     &                             * ((-1.d0)**(iabs(jtot+iomreac)))
                                x1=x1+ss*pm(jtot2)
                                x1=x1/dsqrt(2.d0)
                             endif
                             if(iomprod.ne.0)x1=x1*dsqrt(2.d0)
                             dtrans(ir,iang,iomprod,ican)=x1
                           enddo
                        enddo          
                     endif ! angtrans
                  endif ! angb

               endif ! ifail
            endif ! Jtot
         enddo  ! ir (irpa o iRga)
      enddo ! iang

    
      return
      end  subroutine bfreac2prod

!--------------------------------------------------

      subroutine prodpaq
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      use mod_absorcion_01y2
      use mod_photoini_01y2
      use mod_colini_01y2
      use mod_flux_01y2
      implicit none

      integer :: iomb,j,iv,icanp,ican,ielec,iom,jang,ir2,itot,iang,jinip
      integer :: nsend,ierr,ir,jangp
      real*8 :: xxx,x2

      real*8 :: funprod(n2prod0:n2prod1,nangu2,ncanmax)
      real*8 :: funprodp(n2prod0:n2prod1,nangu2,ncanmax)
      real*8 :: xproj(jiniprod:jmaxprod)
      real*8 :: aux(npun2)
      real*8 :: aux2(n2prod0:n2prod1,nangproj0:nangproj1)
      include "mpif.h"

      Cvjprod(:,:,:)=0.d0
      funprod(:,:,:)=0.d0
      funprodp(:,:,:)=0.d0

* transformation: from reactant to product Jacobi coordinates
*1)  from rp1,Rg1,gamma1 --> Rg2,Rg1,gamma1
      do icanp=1,ncanproc(idproc)
         ican=ibasproc(icanp,idproc)
         ielec=nelebas(ican)
         iom=iombas(ican)
         do iang=nangproj0,nangproj1
         do jangp=1,nanguproc
         jang=indangreal(jangp,idproc)
         if(jang.eq.iang)then

            aux(:)=0.d0

            do ir=1,npunreal(iang)
               ir2=indr2(ir,iang)
               itot=indtotproc(ir,icanp,jangp,idproc)
               aux(ir2)=aux(ir2)+transRgjac(itot)*rpaqproc(itot)
            enddo

            do ir2=n2prod0,n2prod1
               funprodp(ir2,iang,ican)=aux(ir2)
            enddo
      
         endif  ! jangp = iang
         enddo  ! jangp
         enddo ! iang
      enddo ! icanp

      nsend=(n2prod1-n2prod0+1)*nangu2*ncanmax
      call MPI_REDUCE(funprodp,funprod,nsend,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)

*2) Projection in product state

      if(IDPROC.EQ.0)THEN
         do ican=1,ncan
           iom=iombas(ican)
           ielec=nelebas(ican)
           do iomb=iomminprod,iommaxprod
              jinip=max0(jiniprod,iomb)
              do iang=nangproj0,nangproj1
              do ir2=n2prod0,n2prod1
                 aux2(ir2,iang) = funprod(ir2,iang,ican)
     &                       * dtrans(ir2,iang,iomb,ican)
              enddo
              enddo
*bf transformation angle
              do j=jinip,jmaxprod
              do iv=nviniprod,nvmaxprod
                 xxx=0.d0
                 
                 do iang=nangproj0,nangproj1
                 do ir2=n2prod0,n2prod1
                    x2=prodwf(ir2,iang,ielec,iv,j,iomb)
                    xxx=xxx+x2*aux2(ir2,iang)
                 enddo
                 enddo
                    
                 Cvjprod(iv,j,iomb)=Cvjprod(iv,j,iomb)+xxx
               enddo
               enddo
            enddo  ! iomb

      enddo ! ican

      ENDIF ! idproc=0

      return
      end subroutine prodpaq
!--------------------------------------------------

      subroutine prodcvj
      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_Hphi_01y2
      use mod_absorcion_01y2
      use mod_photoini_01y2
      use mod_colini_01y2
      use mod_flux_01y2
      implicit none

      integer :: iomprod,jinip,j,iv,ie
      complex*16 :: zzz
      
      include "mpif.h"
     
*b) Balint-Kurti method to extract probabilities
      if(idproc.eq.0)then
         do iomprod=iomminprod,iommaxprod
            jinip=max0(jiniprod,iomprod)
            do j=jinip,jmaxprod
            do iv=nviniprod,nvmaxprod
               zzz=dcmplx(Cvjprod(iv,j,iomprod),0.d0)
               do ie=1,netot
                  zS2prod(ie,iv,j,iomprod)=
     &              zS2prod(ie,iv,j,iomprod)+zzz*zCkcheby(iE)
               enddo
            enddo
            enddo
         enddo

      endif ! idproc = 0

      return
      end subroutine prodcvj
!--------------------------------------------------
!=======================================================
      end module mod_coortrans01_02
