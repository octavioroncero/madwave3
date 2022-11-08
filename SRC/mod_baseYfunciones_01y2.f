      module mod_baseYfunciones_01y2
      use mod_gridYpara_01y2
      use mod_pot_01y2
      implicit none 

!---------------------------------------------------------------------------------------!
! Determination of  base for the AB + C  in a grid on R=R2, r=R1,gam                    !
!                                01 + 2                                                 !
!                                                                                       !
!     basis: determines basis set quantum number                                        !
!     angular_functions: determines the d^j_mmp functions for reactants and products    !
!     radial_functions_write: phi_vj(r=R1) functions describing 01=AB functions         !
!                     _read: read the functions calculated in previous rutine           !
!     product_radialf_write: phi^j_k(r,electr.) calculated in several electronic states !
!                    _read: reads products radial functions                             !
!_______________________________________________________________________________________!

      save
*     public ! 
      integer, allocatable :: invbc(:),j00(:)
      integer, allocatable :: jbas(:,:),nojbas(:)
      integer, allocatable :: noBCstates(:,:)
      integer, allocatable :: noBCprod(:)

      real*8 :: ediatref
      real*8, allocatable :: Djmmp(:,:,:,:),Djprod(:,:,:,:),pm(:)
      real*8, allocatable :: ediat(:,:,:),fd(:,:,:,:)
      real*8, allocatable :: ediatprod(:,:),fdprod(:,:,:,:)
      
      integer :: ndim
 !--------------------------------------------------
 !     interface
 !        function splinq(ff,xx,iold,npunt,r1,npuntdim)
 !           integer :: iold,npunt,npuntdim
 !           real*8 :: ff(npuntdim,2),xx(npuntdim),r1
 !        end function splinq
 !     end interface
 !
 !--------------------------------------------------     
      contains
********************************************
*   functions of  mod_baseYfunciones_01y2  *
********************************************

!--------------------------------------------------
      subroutine basis
!--------------------------------------------------
      use mod_gridYpara_01y2
      use mod_pot_01y2 
      
      implicit none
      include "mpif.h"
      integer :: ierror,lmax,i,ip,iparbc,ielec,iom,iomdi,iomat
      integer :: isignexp,isign,ipar,iomtot,ifail,jmin,isi
      integer :: iomind,j,jjj,jmax2,ijump

      allocate(invbc(ncanmax),iombas(ncanmax),j00(ncanmax)
     &       ,jbas(nangu,ncanmax),nojbas(ncanmax),nelebas(ncanmax)
     &     , stat=ierror)

      nointegerproc_mem=nointegerproc_mem+ncanmax*(4+nangu)
      write(6,*)'norealproc_mem= ',norealproc_mem
      write(6,*)'nointegerproc_mem= ',nointegerproc_mem
      call flush(6)
 
**>> electronic/Omega basis set

      if(inc.eq.2)then
         if((-1.d0)**jref.ne.(-1.d0)**j0)then
            write(6,*)' be carefull: j0 do not macth with jref'
            write(6,*)'    is j0  of the same parity as jref???'
            call MPI_BARRIER(MPI_COMM_WORLD, ierror)
            stop
         endif
      endif

      if(idproc.eq.0)then
         write(6,*)
         write(6,*)'    ** Electronic/Rotational basis set **'
         write(6,*)
      endif

      ijump=0
      open(10,file='../pot/cont.pot',status='old',err=1)
      ijump=1
 1    continue
      if(ijump.eq.0)then
          open(10,file='pot/cont.pot',status='old')
      endif
       
      read(10,*)nelec
      do ielec=1,nelec
         read(10,*)iomdiat(ielec),iomatom(ielec)
     &           ,sigdiat(ielec),sigatom(ielec)
      enddo
      close(10)

      if(iommax.lt.iommin.or.iommax.gt.Jtot)then
         write(6,*)' !!Be carefull with iommin= ',iommin,
     &             '                     iommax= ',iommax,
     &        '             and     Jtot  = ',Jtot
         call flush(6)
         call MPI_BARRIER(MPI_COMM_WORLD, ierror)
         stop
      endif

      ncan=0
      lmax=inc*(nangu-1)
      jmax2=lmax*2

      iparbc=(-1)**j0
      do ielec=1,nelec
      do iom=iommin,iommax
         iomdi=iomdiat(ielec)
         iomat=iomatom(ielec)
         isignexp=(iom-iomdi+Jtot)
         isign=1
         if(mod(isignexp,2).ne.0)isign=-1
         ipar=sigdiat(ielec)*sigatom(ielec)*isign
         iomtot=iom*iomdi*iomat

         ifail=0
c         if(iomtot.eq.0.and.ipar.ne.ipartot)ifail=1
         if(ifail.eq.0)then

           jmin=j0
           if(inc.eq.2)then
           isi=1
c             isi=(-1.d0)**(iomdi)
           if(isi*sigdiat(ielec).eq.iparbc)then
              jmin=0
           else
              jmin=1
           endif
         endif
             
 10      continue
         if(jmin.lt.iabs(iom-iomat).or.jmin.lt.iabs(iomdi))then
               jmin=jmin+inc
               go to 10
         endif

         ncan=ncan+1
         invbc(ncan)=(-1)**(jmin+iom-iomat)
         iombas(ncan)=iom
         nelebas(ncan)=ielec
         j00(ncan)=jmin

         if(idproc.eq.0)then
             write(6,"(' ielec= ',i3,'  Omeg_diat= ',i3,'  Omeg_at= '
     &                    ,i3, ' Omega= ',i3)")ielec,iomdi,iomat,iom
         endif
            
         do j=1,nangu
            jjj=jmin+(j-1)*inc
            if(jjj.le.inc*(nangu-1))then 
               jbas(j,ncan)=jjj
               nojbas(ncan)=j
            endif
         enddo
         if(idproc.eq.0)write(6,*)"        j's from  ",j00(ncan)
     &               ,"  to jmax",jbas(nojbas(ncan),ncan)
     &               ,"  in intervals of ",inc
         endif                     ! ifail

         if(ielec.eq.ielecref)then
         if(iomref.lt.0.or.iom.eq.iomref)then
         if(jref.lt.jmin)then
         write(6,*)' jref < jini for ielecref,iomref=',ielecref,iomref
         call flush(6)
         endif
         endif
         endif
      enddo  ! iom
      enddo  ! ielec

      if(ncan.eq.0)THEN
           write(6,*)' !!NCAN= ',ncan,' in INPUT !!'
           call flush(6)
           call MPI_BARRIER(MPI_COMM_WORLD, ierror)
           stop
      endif
      write(6,'(//,10x,"ncana= ",i4,//)')ncan

      do ip=0,nproc-1
         write(6,*)' ip=',ip,'  ncan_proc= ',ncanproc(ip)
         if(ncanproc(ip).gt.ncanprocdim)then
            write(6,*)' for ip= ',ip,' ncan= '
     &           ,ncanproc(ip),' > ncanprocdim'
            call flush(6)
            call MPI_BARRIER(MPI_COMM_WORLD, ierror)
            stop
         endif
      enddo
      
      return
      end subroutine basis
      
!--------------------------------------------------
      subroutine angular_functions
!--------------------------------------------------
      use mod_gridYpara_01y2
      use mod_pot_01y2 
      implicit none
      include "mpif.h"

      integer :: ierror,ielec,iom,m2,mp2,mmax,lmax,jprin,j,iang
      integer :: ican,i,iip,jp,jmax2,j2

      real*8 :: y,gamma,xn,xxx,uno,criterj

      criterj=1.d-13
      lmax=(nangu-1)*inc
      jmax2=lmax*2 
      ndim=3*jmax2
      allocate(Djmmp(nangu,0:nangu2,nelecmax,0:Jtot)
     &       , Djprod(nangu2,0:nangu2,nelecmax,0:Jtot)
     &       , pm(0:ndim)
     &     , stat=ierror)
      norealproc_mem=norealproc_mem
     &     + nangu*(nangu2+1)*nelecmax*(Jtot+1)
     &     + nangu2*(nangu2+1)*nelecmax*(Jtot+1)
     &     +ndim+1
      write(6,*)'norealproc_mem= ',norealproc_mem
      write(6,*)'nointegerproc_mem= ',nointegerproc_mem
      call flush(6)

**>> Evaluation of angular functions: d^j_(Omega-Omg_at),Omg_diat (gamma)

      do ielec=1,nelec
      do iom=0,Jtot
         m2=2*(iom-iomatom(ielec))
         mp2=2*iomdiat(ielec)
         mmax=max0(m2,mp2)/2
         jprin=max0(mmax,jini)

         if(nangu.eq.1)then
           do j=jprin,lmax
              j2=2*j
              Djprod(1,j,ielec,iom)=1.d0
              Djmmp(1,j,ielec,iom)=1.d0
           enddo
         else
            do iang=1,nangu2
               y=cgamprod(iang)
               gamma=dacos(y)
               call dwigner(pm,jmax2,m2,mp2,gamma,ndim)
              
               do j=jprin,lmax
                  j2=2*j
                  xn=dsqrt((2.d0*dble(j)+1.d0)*0.5d0)
                  xxx=xn*dsqrt(weiprod(iang))
                  Djprod(iang,j,ielec,iom)=pm(j2)*xxx
                  if(iang.le.nangu)then
                     xxx=xn*dsqrt(wreal(iang))
                     Djmmp(iang,j,ielec,iom)=pm(j2)*xxx
                  endif
               enddo
            enddo
         endif
      enddo
      enddo

**>> check of normalization

      if(idproc.eq.0)then
      write(6,*)'  ** checking norm. of wigner rot. mat. d^j_m m" **'
      write(6,*)'           only prints errors > ',criterj

         do ican=1,ncan
            write(6,*)'  Omg= ',iombas(ican)
     &               ,'  Omg_at= ',iomatom(nelebas(ican))
     &               ,'  Omg_di= ',iomdiat(nelebas(ican))
            ielec=nelebas(ican)
            iom=iombas(ican)
            do i=1,nojbas(ican)
               j=jbas(i,ican)
               if(j.le.lmax)then
               do iip=i,nojbas(ican)
                  jp=jbas(iip,ican)
                  if(jp.le.lmax)then
                     xxx=0.d0
                     do iang=1,nangu
               xxx=xxx+Djmmp(iang,j,ielec,iom)*Djmmp(iang,jp,ielec,iom)
                     enddo
                     uno=0.d0
                     if(i.eq.iip)uno=1.d0
                  if(dabs(uno-xxx).gt.criterj)write(6,*)j,jp,ican,xxx
                  endif
               enddo
               endif
            enddo
         enddo
      endif
          
      return
      end subroutine angular_functions
!--------------------------------------------------
      subroutine radial_functions01_write
!--------------------------------------------------
      use mod_gridYpara_01y2
      use mod_pot_01y2 
      
      implicit none
      include "mpif.h"
      integer, parameter:: npunt=2000
     
      real*8 :: xx(npunt),V(npunt),potmat(nelecmax,nelecmax)
      real*8 :: rrefA(nelecmax),vv(npunt),beta(npunt),alpha(npunt)
      real*8 :: p(npunt),ff(npunt,2)
      real*8 :: rmis,rfin,ah,rg,ctet,xm,xz,r,rref,be,vasin,xl,eps,spl
      real*8 :: e0,e2,vmean,xnorm,r1,xxx,a,b,vmin,xz1
      integer :: ielec,irmin,ir,j,iv,maxit,ifail,iold,ir1,itry,nchan,kv
      integer :: ierror,ie,je

* temporal allocation of matrices

      allocate(noBCstates(jini:jmax,nelecmax)
     &       ,ediat(nvini:nvmax,jini:jmax,nelecmax)
     &       ,fd(npun1,nvini:nvmax,jini:jmax,nelecmax)
     &       , stat=ierror)

      nointegerproc_mem=nointegerproc_mem+(jmax-jini+1)*nelecmax
      norealproc_mem=norealproc_mem
     &    +(nvmax-nvini+1)*(jmax-jini+1)*nelecmax *(1+npun1)
      write(6,*)'norealproc_mem= ',norealproc_mem
      write(6,*)'nointegerproc_mem= ',nointegerproc_mem
      call flush(6)
      
      noBCstates(:,:)=0
      ediat(:,:,:)=0.d0
      fd(:,:,:,:)=0.d0
      xx(:)=0.d0
      v(:)=0.d0
      vv(:)=0.d0
      beta(:)=0.d0
      alpha(:)=0.d0
      rrefA(:)=0.d0
      p(:)=0.d0
      ff(:,:)=0.d0
               eps=1.d-10
               maxit=20

      write(6,*)
      write(6,*)" ** reactant AB^j(k) (01 fragments) wvf energies",
     &              " expanded in all elec. states **"
      write(6,*)
*a)  grid determination (in a.u.: de rest is in zots)

      rmis=rmis1/convl
      rfin=(rmis1+dble(npun1-1)*ah1)/convl
      ah=(rfin-rmis)/dble(npunt-1)
      rg=100.d0
      ctet=0.d0
      xm=xm1reac*convm
      xz=0.5d0/xm

*b)  potential determination

      do ielec=1,nelecmax

         write(6,*)'---------> Ielectronic = ',ielec
         if(idproc.eq.0)then   
            write(name,'("potr.e",i2.2)')ielec
            open(23,file=name,status='unknown')
         endif
         vmin=1.d10
         irmin=1
         do ir=1,npunt
            r=rmis+dble(ir-1)*ah
            xx(ir)=r
            call potelebond(r,rg,ctet,potmat,nelec,nelecmax)
            V(ir)=potmat(ielec,ielec)
            if(v(ir).lt.vmin)then
                vmin=v(ir)
                irmin=ir
            endif
            if(idproc.eq.0)then
              write(23,'(100(1x,e15.7))')r*convl
     &      ,((potmat(ie,je)/conve/ev2cm,ie=je,nelec),je=1,nelec)   !,ielec
            endif
         enddo
         if(idproc.eq.0)close(23)

         rrefA(ielec)=rmis+dble(irmin-1)*ah
         rrefA(ielec)=rrefA(ielec)*convl
         rref=rrefA(ielec)/convl
         be=0.5d0/rref/rref/xm1reac/convm
         be=be/conve
         vasin=v(npunt)
         if(idproc.eq.0)then   
            if(irmin.eq.1)then
               write(6,*)'  Repulsive state: no bound state calc.'
               write(6,*)
            else
               write(6,*)'           Be (eV) = ',be/8065.5d0
     &               ,'  r_e(Angs.)= ',rrefA(ielec)
               write(6,*)'           Vmin(eV)= ',(vmin)/conve/8065.5d0
     &         ,'    Vinfty(eV)= ',(vasin)/conve/8065.5d0
               write(6,*)
            endif
         endif
         if(irmin.ne.1)then
         
*c)  Loop in rotational quantum number

         do j=jini,jmax
            if(idproc.eq.0)write(6,*)'  ** j= ',j,' **'
*d)   eigenvalues

            do ir=1,npunt
               r=rmis+dble(ir-1)*ah
               xl=0.5d0*dble(j*(j+1))/xm
               vv(ir)=v(ir)+xl/r/r
               beta(ir)=1.d0
               alpha(ir)=-2.d0-2.d0*xm*vv(ir)*ah*ah
            enddo

            call tqli(alpha,beta,npunt,npunt)

            do iv=1,npunt
               alpha(iv)=-xz*alpha(iv)/ah/ah
            enddo

**>> orderig of eigenvalues
 
 44         continue
               nchan=0
               do ir=1,npunt-1
                  if(alpha(ir).gt.alpha(ir+1))then
                     a=alpha(ir)
                     b=alpha(ir+1)
                     alpha(ir)=b
                     alpha(ir+1)=a
                     nchan=nchan+1
                  endif
               enddo
            if(nchan.gt.0)go to 44

*e)    eigenfunctions

            do iv=nvini,nvmax
               xz1=2.d0*xm

               ifail=-1
               if(alpha(iv+1).lt.v(npunt))then
                  e0=alpha(iv+1)*xz1
                  do ir=1,npunt
                     r=rmis+dble(ir-1)*ah
                     xl=0.5d0*dble(j*(j+1))/xm
                     vv(ir)=v(ir)+xl/r/r
                     vv(ir)=vv(ir)*xz1
                     p(ir)=0.d0                     
                  enddo

           call schr(e0,rmis,rfin,npunt,maxit,eps,e2,kv,itry,vv,p,npunt)

                  e2=e2/xz1
                  ediat(iv,j,ielec)=e2/conve*conve1
                  noBCstates(j,ielec)=iv
                  ifail=0
*f)     interpolation to adapt eigenvalues to the WP grid

                  vmean=0.d0
                  do ir=1,npunt
                     r=rmis+dble(ir-1)*ah
                     xx(ir)=r
                     ff(ir,1)=p(ir)/dsqrt(convl)
                     ff(ir,2)=0.d0
                     vmean=vmean+p(ir)*p(ir)*v(ir)
                  enddo
                  vmean=vmean*ah

                  call splset(ff,xx,npunt,npunt)
                  iold=2
                  
                  xnorm=0.d0
                  do ir1=1,npun1
                     r1=rmis1+dble(ir1-1)*ah1
                     r1=r1/convl
                     fd(ir1,iv,j,ielec)=0.d0
                     if(r1.ge.xx(1).and.r1.le.xx(npunt))then
c                  fd(ir1,iv,j,ielec)=splinq(ff,xx,iold,npunt,r1,npunt)
                        call splinqq(ff,xx,iold,npunt,r1,npunt,spl)
                        fd(ir1,iv,j,ielec)=spl
                     endif
                         
                     xnorm=xnorm+fd(ir1,iv,j,ielec)**2
                  enddo
                  xnorm=xnorm*ah1
                  xxx=dsqrt(ah1)/dsqrt(xnorm)
                  do ir1=1,npun1
                     fd(ir1,iv,j,ielec)=fd(ir1,iv,j,ielec)*xxx
                  enddo
      
                  write(6,'(2(2x,i4),4(2x,e15.7))')iv,kv
     &               ,alpha(iv+1)/conve/eV2cm,e2/conve/eV2cm
     &                 ,xnorm,vmean/conve/eV2cm
                  call flush(6)
                  endif
            enddo
            if(noBCstates(j,ielec).eq.0.and.ifail.eq.-1)then
               if(idproc.eq.0)write(6,*)' elec. stat. doesn"t support'
     &                   ,'  bound states for higher j values' 
               go to 1234
            endif
         enddo  ! j
         endif
 1234    continue               ! in case no bound for higher j's
         
      enddo  ! ielec

! reference energy for iprod.ne.1
 
      ediatref=0.d0
      if(iprod.ne.1.and.npun1.gt.1)then

         ediatref=ediat(nvref,jref,ielecref)
           
         write(6,*)' reference energy in reactants= '
     &                                  ,ediatref/conve1/ev2cm

         if(idproc.eq.0)then
         open(10,file='func/eref',status='new')
         write(10,*)ediatref
         close(10)
         endif
      endif

! writing reactants functions

         if(idproc.eq.0)then
           open(10,file='func/bcwf',status='new')
              do ielec=1,nelecmax
              do j=jini,jmax
                 write(10,*)ielec,j,noBCstates(j,ielec)
                 do iv=nvini,noBCstates(j,ielec)
                    write(10,*)iv,ediat(iv,j,ielec)
                    do ir1=1,npun1
                       write(10,*)fd(ir1,iv,j,ielec)
                    enddo
                 enddo
              enddo
              enddo
           close(10)
         endif

! dealocating 

!      deallocate(noBCstates,ediat,fd)

      return
 1100 format(2x,'v= ',i2,2x,'n.nodes= ',i2,2x,'E(approx)= ',d15.7,
     & 2x,'E(exact)= ',d15.7)
 9994 format(2x,'v= ',i2,2x,'E= ',d15.7,2x,'B= ',d15.7)

      end subroutine radial_functions01_write
!--------------------------------------------------
      subroutine radial_functions01_read
!--------------------------------------------------
      use mod_gridYpara_01y2
      use mod_pot_01y2 
      
      implicit none
      include "mpif.h"
      integer :: ierror,ielec,j,iv,ir1,iielec,jj,iiv

      allocate(noBCstates(jini:jmax,nelecmax)
     &       ,ediat(nvini:nvmax,jini:jmax,nelecmax)
     &       ,fd(npun1,nvini:nvmax,jini:jmax,nelecmax)
     &     , stat=ierror)

      nointegerproc_mem=nointegerproc_mem+(jmax-jini+1)*nelecmax
      norealproc_mem=norealproc_mem
     &    +(nvmax-nvini+1)*(jmax-jini+1)*nelecmax *(1+npun1)
      write(6,*)'norealproc_mem= ',norealproc_mem
      write(6,*)'nointegerproc_mem= ',nointegerproc_mem
      call flush(6)

      noBCstates(:,:)=0
      ediat(:,:,:)=0.d0
      fd(:,:,:,:)=0.d0

*** Reading   R1 wavefunctions ()
********************************
         
      write(6,*)
      write(6,*)" *******************************************"
      write(6,*)'  Reading 01 wavefunctions in file func/bcwf'
      write(6,*)" ********************************************"
      write(6,*)

      open(10,file='../func/bcwf',status='old')
         do ielec=1,nelecmax
         do j=jini,jmax
            read(10,*)iielec,jj,noBCstates(j,ielec)
            do iv=nvini,noBCstates(j,ielec)
               read(10,*)iiv,ediat(iv,j,ielec)
               do ir1=1,npun1
                  read(10,*)fd(ir1,iv,j,ielec)
               enddo
            enddo
         enddo
         enddo
      close(10)
              
!     substracting ref energy

      ediatref=0.d0
!      if(npun1.gt.1)then
         open(10,file='../func/eref',status='old')
         read(10,*)ediatref
         close(10)
!      endif                 
      write(6,*)' substracting ediatref to the diatomic energies'
 
      do ielec=1,nelec
      do j=jini,jmax
         do iv=nvini,noBCstates(j,ielec)
            ediat(iv,j,ielec)=ediat(iv,j,ielec)-ediatref
         enddo
      enddo
      enddo

!     writting energies

      write(6,*)
      write(6,*)'ie  j   v   E(eV) '
      write(6,*)'-------------------'
      write(6,*)
      do ielec=1,nelecmax
      do j=jini,jmax
         if(noBCstates(j,ielec).ge.nvini)then
            do iv=nvini,noBCstates(j,ielec)
c            if(ediat(iv,j,ielec).gt.ediatref)then
               write(6,*)ielec,j,iv,ediat(iv,j,ielec)/(conve1*eV2cm)
               call flush(6)
c            endif
            enddo
         endif
      enddo
      enddo
      
      return
      end subroutine radial_functions01_read
!--------------------------------------------------
      subroutine radial_functions01eq_write
!--------------------------------------------------
      use mod_gridYpara_01y2
      use mod_pot_01y2 
      
      implicit none
      include "mpif.h"
      integer, parameter:: npunt=2000
     
      real*8 :: potmat(nelecmax,nelecmax)
      real*8 :: rrefA(nelecmax)
      real*8 :: rmis,rfin,ah,rg,ctet,xm,xz,r,rref,be,vasin,xl,eps,spl
      real*8 :: e0,e2,vmean,xnorm,r1,xxx,a,b,vmin,xz1
      integer :: ielec,irmin,ir,j,iv,maxit,ifail,iold,ir1,itry,nchan,kv
      integer :: ierror,ie,je

* temporal allocation of matrices

      allocate(noBCstates(jini:jmax,nelecmax)
     &       ,ediat(nvini:nvmax,jini:jmax,nelecmax)
     &       ,fd(npun1,nvini:nvmax,jini:jmax,nelecmax)
     &       , stat=ierror)

      nointegerproc_mem=nointegerproc_mem+(jmax-jini+1)*nelecmax
      norealproc_mem=norealproc_mem
     &    +(nvmax-nvini+1)*(jmax-jini+1)*nelecmax *(1+npun1)
      write(6,*)'norealproc_mem= ',norealproc_mem
      write(6,*)'nointegerproc_mem= ',nointegerproc_mem
      call flush(6)
      
      noBCstates(:,:)=0
      ediat(:,:,:)=0.d0
      fd(:,:,:,:)=0.d0
      rrefA(:)=0.d0

      write(6,*)
      write(6,*)" ** reactant AB^j(k) (01 fragments) wvf energies",
     &              " expanded in all elec. states **"
      write(6,*)
*a)  grid determination (in a.u.: de rest is in zots)

      r=rmis1/convl
      rg=100.d0
      ctet=0.d0
      xm=xm1reac*convm

      if(nvmax.gt.nvini)then
         write(6,*)'  nvmax= ',nvmax,' changed to nvini= ',nvini
         nvmax=nvini
      endif

*b)  potential determination

      do ielec=1,nelecmax

         write(6,*)'---------> Ielectronic = ',ielec
            call potelebond(r,rg,ctet,potmat,nelec,nelecmax)

         rrefA(ielec)=r
         rrefA(ielec)=rrefA(ielec)*convl
         rref=rrefA(ielec)/convl
         be=0.5d0/rref/rref/xm1reac/convm
         be=be/conve
         if(idproc.eq.0)then   
            write(6,*)'           Be (eV) = ',be/8065.5d0
     &               ,'  r_e(Angs.)= ',rrefA(ielec)
     &               ,'  V_e (eV)= ',potmat(ielec,ielec)/conve/ev2cm  
         endif
         
*c)  Loop in rotational quantum number

         do j=jini,jmax
            if(idproc.eq.0)write(6,*)'  ** j= ',j,' **'
*d)   eigenvalues and eigenfunctions
            noBCstates(j,ielec)=nvmax
            do iv=nvini,nvmax
                  xl=0.5d0*dble(j*(j+1))/xm
                  e2=xl/r/r+potmat(ielec,ielec)
                  ediat(iv,j,ielec)=e2/conve*conve1
                  do ir1=1,npun1
                     fd(ir1,iv,j,ielec)=1.d0
                  enddo
      
                  write(6,'(2(2x,i4),4(2x,e15.7))')iv,kv
     &               ,e2/conve/eV2cm
                  call flush(6)
            enddo ! iv=nvini=nvmax
         enddo  ! j
         
      enddo  ! ielec

 
      ediatref=0.d0

         ediatref=ediat(nvref,jref,ielecref)
           
         write(6,*)' reference energy in reactants= '
     &                                  ,ediatref/conve1/ev2cm

         if(idproc.eq.0)then
         open(10,file='func/eref',status='new')
         write(10,*)ediatref
         close(10)
         endif

! writing reactants functions

         if(idproc.eq.0)then
           open(10,file='func/bcwf',status='new')
              do ielec=1,nelecmax
              do j=jini,jmax
                 noBCstates(j,ielec)=nvmax
                 write(10,*)ielec,j,noBCstates(j,ielec)
                 do iv=nvini,noBCstates(j,ielec)
                    write(10,*)iv,ediat(iv,j,ielec)
                    do ir1=1,npun1
                       write(10,*)fd(ir1,iv,j,ielec)
                    enddo
                 enddo
              enddo
              enddo
           close(10)
         endif

! dealocating 

!      deallocate(noBCstates,ediat,fd)

      return
 1100 format(2x,'v= ',i2,2x,'n.nodes= ',i2,2x,'E(approx)= ',d15.7,
     & 2x,'E(exact)= ',d15.7)
 9994 format(2x,'v= ',i2,2x,'E= ',d15.7,2x,'B= ',d15.7)

      end subroutine radial_functions01eq_write
!--------------------------------------------------
!--------------------------------------------------
      subroutine product_radialf_write
      use mod_gridYpara_01y2
      use mod_pot_01y2
      
      implicit none
      include "mpif.h"

      real*8 :: potmatrix(n2prod1,nelecmax,nelecmax)
     &     ,potmat(nelecmax,nelecmax)
     &     ,fun(n2prod1,nelecmax,nviniprod:nvmaxprod)
     &     ,eval(nviniprod:nvmaxprod)
     &     ,T(nelecmax,nelecmax),eigen(nelecmax)

      real*8 :: rmis,rfin,ah,rg,ctet,xm,xz,r,xnorm
      integer :: ierror,ir,ie,je,j,iv,ielec,ir2

      allocate(ediatprod(nviniprod:nvmaxprod,jiniprod:jmaxprod)
     &  ,fdprod(n2prod1,nviniprod:nvmaxprod,jiniprod:jmaxprod,nelecmax)
     &     ,noBCprod(jiniprod:jmaxprod)
     &       , stat=ierror)

      nointegerproc_mem=nointegerproc_mem+(jmaxprod-jiniprod+1)
      norealproc_mem=norealproc_mem
     &  +(nvmaxprod-nviniprod+1)*(jmaxprod-jiniprod+1)
     &  +(nvmaxprod-nviniprod+1)*(jmaxprod-jiniprod+1)*n2prod1*nelecmax
      write(6,*)'norealproc_mem= ',norealproc_mem
      write(6,*)'nointegerproc_mem= ',nointegerproc_mem
      call flush(6)
      
      write(6,*)
      write(6,*)" ** product BC^j(k) (02 fragments) wvf energies",
     &              " expanded in all elec. states **"
      write(6,*)
!a)  grid determination (in a.u.: de rest is in zots)

         rmis=rmis2/convl
         ah=ah2/convl
         rg=100.d0
         ctet=0.d0
         xm=xm1prod*convm

!b)  potential determination

         if(idproc.eq.0)then
            open(11,file='potadi-prod02.dat',status='unknown')
            open(10,file='potdia-prod02.dat',status='unknown')
         endif
         
         do ir=1,n2prod1
            r=rmis+dble(ir-1)*ah
            call potelebond(rg,r,ctet,potmat,nelec,nelecmax)
            do ie=1,nelec
            do je=1,nelec
               potmatrix(ir,ie,je)=potmat(ie,je)
            enddo
            enddo

            call diagon(potmat,nelec,nelecmax,T,eigen)

            if(idproc.eq.0)then
                  write(10,'(100(1x,e15.7))')r*convl
     &       ,(potmatrix(ir,ie,ie)/conve/ev2cm,ie=1,nelec)
                  write(11,'(100(1x,e15.7))')r*convl
     &       ,(eigen(ie)/conve/ev2cm,ie=1,nelec)
            endif
            
         enddo
         rfin=rmis+dble(n2prod1-1)*ah
         if(idproc.eq.0)then
            close(10)
            close(11)
         endif
         
!c)  Loop in rotational quantum number

         do j=jiniprod,jmaxprod
            if(idproc.eq.0)write(6,*)'  ** j= ',j,' **'
!d)   eigenvalues
            call  bndbcele(Eval,fun,potmatrix,xm
     &         ,rmis,rfin,nviniprod,nvmaxprod,j,n2prod1,nelec)

!e)    eigenfunctions

            do iv=nviniprod,nvmaxprod
               xnorm=0.d0
               ediatprod(iv,j)=eval(iv)/conve*conve1
               do ielec=1,nelec
               do ir2=1,n2prod1
                   fdprod(ir2,iv,j,ielec) = fun(ir2,ielec,iv)
     &                                    * dsqrt(ah2/convl)
                    xnorm=xnorm+fdprod(ir2,iv,j,ielec)**2
               enddo
               enddo
               xnorm=xnorm
               write(6,'(1(2x,i4),4(2x,e15.7))')iv
     &              ,ediatprod(iv,j)/conve1/eV2cm,xnorm
               call flush(6)
            enddo ! iv

            ielec=1
            noBCprod(j)=nvmaxprod
         enddo  ! j

!     reference energy for iprod.eq.1
 

         ediatref=ediatprod(nvref,jref)
         if(iprod.eq.1)then
            write(6,*)' reference energy in products= '
     &                   ,ediatref/conve1/8065.5d0
            if(idproc.eq.0)then
            open(10,file='func/eref',status='new')
            write(10,*)ediatref
            close(10)
            endif
         endif

! writing products functions
            
         if(idproc.eq.0)then
            open(10,file='func/prodwf',status='new')
            do j=jiniprod,jmaxprod
               write(10,*)j,noBCprod(j),n2prod1,nelecmax
               do iv=nviniprod,noBCprod(j)
                  write(10,*)iv,ediatprod(iv,j)
                  do ir2=1,n2prod1
                  write(10,*)(fdprod(ir2,iv,j,ielec),ielec=1,nelecmax)
                  enddo
               enddo
            enddo
            close(10)

            do j=jiniprod,jiniprod
               do iv=nviniprod,noBCprod(j)
                  write(name,'("wv02-v",i2.2,"j0.dat")')iv
                  open(10,file=name,status='new')
                  write(10,*)iv,ediatprod(iv,j)/conve/ev2cm 
                  do ir2=1,n2prod1
                     r=rmis2+dble(ir2-1)*ah2
                     write(10,*)r
     &               ,(fdprod(ir2,iv,j,ielec),ielec=1,nelecmax)
                  enddo
                  close(10)
               enddo
            enddo
            
         endif
         
!     deallocate matrices
         
!      deallocate(noBCprod,ediatprod,fdprod)

      return
 1100 format(2x,'v= ',i2,2x,'E(eV)= ',d15.7,
     & 2x,'Norm= ',d15.7)
 9994 format(2x,'v= ',i2,2x,'E= ',d15.7,2x,'B= ',d15.7)
      end subroutine product_radialf_write
!--------------------------------------------------
      subroutine product_radialf_read
      use mod_gridYpara_01y2
      use mod_pot_01y2 
      
      implicit none
      include "mpif.h"

      integer :: ierror,ielec,j,iv,ir2,iielec,jj,iiv,nn2prod1,nnelec
      real*8 :: r2,xnorm

      allocate(ediatprod(nviniprod:nvmaxprod,jiniprod:jmaxprod)
     &  ,fdprod(n2prod1,nviniprod:nvmaxprod,jiniprod:jmaxprod,nelecmax)
     &    ,noBCprod(jiniprod:jmaxprod)
     &       , stat=ierror)
      
      nointegerproc_mem=nointegerproc_mem+(jmaxprod-jiniprod+1)
      norealproc_mem=norealproc_mem
     &  +(nvmaxprod-nviniprod+1)*(jmaxprod-jiniprod+1)
     &  +(nvmaxprod-nviniprod+1)*(jmaxprod-jiniprod+1)*n2prod1*nelecmax
      write(6,*)' in  product_radialf_read '
     &   ,' norealproc_mem=',norealproc_mem
     &     ,' nointegerproc_mem=',nointegerproc_mem
     &     ,' memory(Gb)= '
     &   ,dble(nointegerproc_mem*4+norealproc_mem*8)*1.d-9
      call flush(6)

***   Reading   02 wavefunctions ()
********************************
      write(6,*)
      write(6,*)" *************************************"
      write(6,*)" ** Reading 02 product BC^j(k) wvf  **"
      write(6,*)" **  expanded in all elec. states.  **"
      write(6,*)" *************************************"
      write(6,*)

      fdprod(:,:,:,:)=0.d0
      open(10,file='../func/prodwf',status='old')
         do j=jiniprod,jmaxprod
            read(10,*)jj,noBCprod(j),nn2prod1,nnelec
            if(nvmaxprod.lt.noBCprod(j))then
               write(6,*)' nvmaxprod = ',nvmaxprod
     &              ,'<noBCprod(j)=',noBCprod(j)
     &              ,' in ../func/prodwf: please increase it input.dat'
               call flush(6)
               call MPI_BARRIER(MPI_COMM_WORLD, ierror)
               stop
            endif
            if(nn2prod1.ne.n2prod1)then
               write(6,*)' in product_radialf_read '
               write(6,*)'  in file ../func/prodwf '
               write(6,*)'     n2prod1= ',nn2prod1
               write(6,*)' while in input.dat = ',n2prod1
               call flush(6)
               call MPI_BARRIER(MPI_COMM_WORLD, ierror)
               stop
            endif
            do iv=nviniprod,noBCprod(j)
               read(10,*)iiv,ediatprod(iv,j)
               xnorm=0.d0
               do ir2=1,nn2prod1
                  read(10,*)(fdprod(ir2,iv,j,ielec),ielec=1,nnelec)
                  do ielec=1,nnelec
                     xnorm=xnorm+fdprod(ir2,iv,j,ielec)**2
                  enddo
               enddo
!               if(j == jiniprod)then
!                  write(name,'("checkprod.v",i2.2,".dat")')iv
!                  open(11,file=name,status='unknown')
!                  write(11,*)iv,xnorm,ediatprod(iv,j)
!                  do ir2=1,nn2prod1
!                     r2=rmis2+dble(ir2-1)*ah2
!                     write(11,*)r2,(fdprod(ir2,iv,j,ielec),ielec=1,nnelec)
!                  enddo
!                  close(11)
!               endif
            enddo
         enddo
      close(10)

!     substracting ref energy

      open(10,file='../func/eref',status='old')
      read(10,*)ediatref
      close(10)
         
      write(6,*)' substracting ediatref to the products energies '
       
      do j=jiniprod,jmaxprod
      do iv=nviniprod,noBCprod(j)
         ediatprod(iv,j)=ediatprod(iv,j)-ediatref
      enddo
      enddo

      
!     writting energies

      write(6,*)
      write(6,*)' j   v   E(eV) '
      write(6,*)'-------------------'
      write(6,*)
      do j=jiniprod,jmaxprod
         if(noBCprod(j).gt.0)then
         do iv=nviniprod,noBCprod(j)
            write(6,*)j,iv,ediatprod(iv,j)/(conve1*eV2cm)
         enddo
         endif
      enddo
      
      return
      end subroutine product_radialf_read
!--------------------------------------------------
!--------------------------------------------------
      end module mod_baseYfunciones_01y2
