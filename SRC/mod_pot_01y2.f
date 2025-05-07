      module mod_pot_01y2
      
!---------------------------------------------------------------------------------!
! calculation of the AB + C potential in a grid on R=R2, r=R1,gam                 !
!                    01 + 2                                                       !
!  pot0 initialices                                                               !
!  pot1 calculate and write the potential for those values < Vcutmax with indexes !
!  pot2 reads the potential and indexes to perform  calculations                  !
!_________________________________________________________________________________!

      use mod_gridYpara_01y2
      implicit none 

      save
*     public ! for input data in namelist, potential and reduced grid
      
      real*8 :: vmintot,emaxtot,emindlt,delta2  ! for cheby propagations
      real*8 :: vmaxtot                         ! for cheby propagations

* mata potential and kinetic terms
      real*8 ::  vcutmaxeV,radcutmaxeV,rotcutmaxeV
      real*8 ::  vcutmax  ,radcutmax  ,rotcutmax
* pot
      integer :: nomaxV,maxpoint
      integer :: nelec,radau,NO_ref_energy
      integer,allocatable :: iomdiat(:),iomatom(:)
      real*8, allocatable :: sigdiat(:),sigatom(:)
      real*8, allocatable :: VVV(:,:,:,:)
* masses
      real*8 :: xm1red, xm2red, xm0red,xmtot
     &                ,xm1reac,xm2reac,xm1prod,xm2prod
     &     ,xm1,xm0,xm2
      real*8 :: R1inf_radial_functions, R2inf_radial_functions

      character*70 :: name
      character*10 :: system

*     indexes and dimensions

      integer, allocatable :: npunreal(:)
      integer :: npuntot,ngridtot,nbastot,nbastotproc
      integer, allocatable :: indtotproc(:,:,:,:)
     &     ,indcanproc(:),indangproc(:),indrgridproc(:)
     &     ,indr1(:,:),indr2(:,:)

      contains
**********************************
*   functions of mod_pot_01y2    *
**********************************
!--------------------------------------------------
      subroutine pot0
!--------------------------------------------------      
      use  mod_gridYpara_01y2
      implicit none
      integer :: ierror,ir1,ir2
*********************************************************
      namelist /inputpotmass/system,xm1,xm0,xm2
     &     ,VcutmaxeV,radcutmaxeV,rotcutmaxeV
     &     ,radau,R1inf_radial_functions,R2inf_radial_functions
     &     ,No_ref_energy


         write(6,'(40("_"),/,10x,"Pot_mod",/,40("_"))')
         open(10,file='input.dat',status='old')
         radau=0
         R1inf_radial_functions=100.d0
         R2inf_radial_functions=100.d0
         NO_ref_energy=1
         read(10,nml = inputpotmass)
         write(6,'(80("-"),/,10x
     &      ,"Mass and pot determination for 01+2= ",a20
     &      ,/80("-"),/)')system
         write(6,nml = inputpotmass)
         close(10)

         nelec=nelecmax
!**>> masses
!* reduced masses (amu are the mass units using zots, in this program

      xmtot=xm0+xm1+xm2

      if(radau.eq.0)then
            write(6,*)' Using A + BC body fixed Jacobi coordinates'
* Reactant Jacobi
         xm1reac = (xm0*xm1)/(xm0+xm1)
         xm2reac = (xm2*(xm0+xm1))/xmtot
* Product Jacobi
         xm1prod = (xm0*xm2)/(xm0+xm2)
         xm2prod = (xm1*(xm0+xm2))/xmtot
* Reduced masses for dynamical calculations

         xm1red = (xm0*xm1) / (xm0+xm1)
         xm2red = (xm2*(xm0+xm1)) / xmtot
         xm0red = 0.d0
      elseif(radau.eq.1)then
         write(6,*)' Using body fixed Radau(2+1) coordinates'
         xm1red=xm1
         xm2red=xm2
         xm1reac=xm1red
         xm2reac=xm2red
         xm1prod=xm1red
         xm2prod=xm1red
      endif
! converting quantities to zots, masses are in amus

      write(6,*)' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(6,*) '  Pot is cut at E(eV)= ',vcutmaxeV
      write(6,*) '      above energy of ref. state'
      write(6,*)' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      vcutmax=vcutmaxeV*ev2cm*conve1
      radcutmax=radcutmaxeV*ev2cm*conve1
      rotcutmax=rotcutmaxeV*ev2cm*conve1

**>> PES initialization, and electronic states parameters
*  iomdiat: diatomic electronic angular momentum projection in the diatomic axis
*  Iomatom:  atomic electronic angular momentum projection in the triatomic axis
*  sigdiat: eigenvalue of the sigma_xz of the diatomic electronic f.
*  sigatom: eigenvalue of the sigma_xz of the atomic electronic f.

      allocate( iomdiat(nelecmax),iomatom(nelecmax)
     &        ,sigdiat(nelecmax),sigatom(nelecmax)
     &        ,npunreal(nangu)
     &       , stat=ierror)

      nointegerproc_mem=nointegerproc_mem+2*nelecmax+nangu
      norealproc_mem=norealproc_mem+2*nelecmax

      return
      end subroutine pot0

!--------------------------------------------------
      subroutine pot1
!--------------------------------------------------
      use mod_gridYpara_01y2
      implicit none
      include "mpif.h"
      
      real*8 :: vmin,x1min,x2min,angmin,ctet,r1,r2,rpeq,Rg
      integer :: ie,ir1,ir2,iang,maxpoint,ielec,jelec,icount,maxpointang
      integer :: je,ke,nsend,ierror,icorte,ieleccut
      real*8 :: x1,gam,x2,calpha,pp,vref,coefmax
      real*8 :: vmat(nelecmax,nelecmax),potmat(nelecmax,nelecmax)
      real*8 :: vdiagon(nelecmax,nelecmax)
      real*8 :: Tpot(nelecmax,nelecmax),eigenpot(nelecmax)

      vmaxtot=vcutmax
      npunreal(:)=0
      if(idproc.eq.0)then
         write(6,*)
         write(6,*)'   ** Calculating Potential **'
         write(6,*)
         write(6,*)'         only in idproc=0 '
         open(10,file='func/eref',status='old')
         read(10,*)vref
         close(10)
         write(6,*)'       substracting Eref(vref)= ',vref/conve1/eV2cm 

         vmin=1.d10
         x1min=1.d10
         x2min=1.d10
         angmin=1.d10
         maxpoint=0

         do iang=1,nangu
            ctet=cgamma(iang)
            write(6,*)'iang=', iang
            call flush(6)
            write(name,'("pot/pot.",i3.3,".dat")')iang
            write(6,*)' writting file = ',name
            call flush(6)
            open(10,file=name,status='new')
            npunreal(iang)=0
            do ir2=1,npun2
            r2=rmis2+dble(ir2-1)*ah2
            do ir1=1,npun1
               r1=rmis1+dble(ir1-1)*ah1

               rpeq=r1/convl  ! to call pot in a.u.
               Rg=r2/convl
               ctet=cgamma(iang)
               if(radau.eq.0)then

                  X1=rpeq
                  gam=xm1/(xm0+xm1)
                  X2=Rg*Rg+gam*gam*rpeq*rpeq+2.d0*gam*Rg*rpeq*ctet
                  X2=dsqrt(X2)
                  calpha=(Rg*ctet+gam*rpeq)/X2
                  if(dabs(calpha).gt.1.d0)then
                     calpha=calpha/dabs(calpha)
                  endif
               elseif(radau.eq.1)then
                  x1=rpeq
                  x2=Rg
                  calpha=ctet
               endif
               call potelebond(x1,x2,calpha,potmat,nelecmax,nelecmax)

               do ie=1,nelecmax
               do je=1,nelecmax
                  vmat(ie,je)=potmat(ie,je)/conve*conve1
                  if(ie.eq.je)vmat(ie,je)=vmat(ie,je)-vref
               enddo
               enddo

               do ielec=1,nelecmax
                  if(vmat(ielec,ielec).gt.vcutmax)then
                     do jelec=1,nelec
                        vmat(ielec,jelec)=0.d0
                        vmat(jelec,ielec)=0.d0
                     enddo
                     vmat(ielec,ielec)=vcutmax
                  endif
               enddo
              
            
               if(nelecmax.eq.1)then
                  Tpot(1,1)=1.d0
                  eigenpot(1)=vmat(1,1)
               else
                  Tpot(:,:)=0.d0
                  eigenpot(:)=0.d0
                  vdiagon=vmat
                  call DIAGON(vdiagon,nelecmax,nelecmax,Tpot,EIGENpot)
               endif

               icorte=0
               do ielec=1,nelecmax
                  if(eigenpot(ielec).gt.vcutmax)icorte=icorte+1
                  if(eigenpot(ielec).gt.vmaxtot)vmaxtot=eigenpot(ielec)
                  
               enddo

                  
                  if(icorte.ge.1)then
!----                     
                       write(6,*)' vcutmax= ',vcutmax,iang,r2,r1
                        write(6,*)' eigenpot= ',eigenpot
                         write(6,*)' Vmat '
                        do ielec=1,nelecmax
                         write(6,*)(vmat(jelec,ielec),jelec=1,nelecmax)
                      enddo
!----                      
                    
                  endif
                  
               icount=0
               do ielec=1,nelecmax              
                  if(eigenpot(ielec).lt.vcutmax)icount=icount+1
                  if(eigenpot(ielec).lt.vmin)then
                      vmin=eigenpot(ielec)
                      x1min=r1
                      x2min=r2
                      angmin=dacos(cgamma(iang))*180.d0/pi
                  endif
               enddo
               
               if(icount.gt.0)then
                   npunreal(iang)=npunreal(iang)+1
                   if(idproc.eq.0)then
                write(10,*)ir1,ir2,((vmat(ielec,jelec),ielec=1,nelecmax)
     &                                                ,jelec=1,nelecmax)
                   endif
               endif ! icount > 0 

               if(npunreal(iang).gt.npun1*npun2
     &          .or.npunreal(iang).lt.0)then
                  write(6,*)' problem with npunreal for iang= ',iang
                  write(6,*)ir1,ir2,npunreal(iang)
               endif
            enddo ! ir1
            enddo ! ir2
 
            close(10)
            maxpoint=maxpoint+npunreal(iang)
         enddo ! iang                
       
         open(10,file='pot/cont.pot',status='unknown')
         write(10,*)nelecmax
         do ielec=1,nelecmax
            write(10,*)iomdiat(ielec),iomatom(ielec)
     &           ,sigdiat(ielec),sigatom(ielec)
         enddo  
         write(10,*)vmin,maxpoint,vmaxtot
         write(10,*)rmis1,rfin1,npun1
         write(10,*)rmis2,rfin2,npun2
         write(10,*)nangu,inc
         write(10,'(1000(1x,i12))')(npunreal(iang),iang=1,nangu)
         close(10)

         write(6,*)' No. of points per angle'
         do iang=1,nangu
            write(6,*)iang,npunreal(iang)
         enddo
         call flush(6)
         
         if(nelecmax.eq.0)then
             write(6,*)'  !!! nelec= ',nelecmax
             call flush(6)
             call MPI_BARRIER(MPI_COMM_WORLD, ierror)
             stop
         endif

         write(6,*)'         Vmax (eV)= ',vmaxtot/conve1/eV2cm
     &        ,' while Vcutmax(eV)= ',vcutmax/conve1/eV2cm
         write(6,*)
         write(6,*)'         Vmin (eV)= ',vmin/conve1/eV2cm
         write(6,*)
         if(radau.eq.0)then
            write(6,*)'      at r1=rpeq (angstroms) = ',x1min
            write(6,*)'         r2=Rgran (angstroms)= ',x2min
         elseif(radau.eq.1)then
            write(6,*)'      at r1 (angstroms) = ',x1min
            write(6,*)'         r2 (angstroms) = ',x2min
         endif
            
         write(6,*)'         gam (degrees)       = ',angmin 
         write(6,*)
         write(6,*)'    ** end potential calculation **'
         call flush(6)
      endif ! idproc.eq.0

      nsend=1
      call MPI_BCAST(vmin,nsend,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

      return
      end subroutine pot1

!--------------------------------------------------
      subroutine pot2
!--------------------------------------------------
      use mod_gridYpara_01y2
      implicit none
      include "mpif.h"
      
      integer :: mpun1,mpun2,mangu,minc,ie,je,ir
      integer :: iang,iang_proc,ip,ican,ican_proc,ijump
      integer :: ierror,nsend
      real*8 :: vmin,xmis1,xfin1,xmis2,xfin2
      real*8 :: vaux(nelecmax,nelecmax)
      integer*8 :: imem
      
!     *> reading initialization data for dimensions
      npunreal(:)=0

      ijump=0
      open(10,file='../pot/cont.pot',status='old',err=1)
      ijump=1
 1    continue
      if(ijump.eq.0)then
         open(10,file='pot/cont.pot',status='old')
      endif
       
      read(10,*)nelec
      do ie=1,nelec
         read(10,*)iomdiat(ie),iomatom(ie)
     &           ,sigdiat(ie),sigatom(ie)
      enddo
      read(10,*)vmin,maxpoint,vmaxtot
      read(10,*)xmis1,xfin1,mpun1
      read(10,*)xmis2,xfin2,mpun2
      read(10,*)mangu,minc
      read(10,*)(npunreal(iang),iang=1,nangu)
      close(10)

      if(npun1.ne.mpun1.or.npun2.ne.mpun2
     &  .or.nangu.ne.mangu.or.inc.ne.minc)then
          write(6,*)' Problem in input no. of points and incj values'
          write(6,*)'   no correspond to those used to generate pot'
          call flush(6)
          call MPI_BARRIER(MPI_COMM_WORLD, ierror)
          stop
      endif
      if(dabs(rmis1-xmis1).gt.1.d-5.or.dabs(rfin1-xfin1).gt.1.d-5
     & .or.dabs(rmis2-xmis2).gt.1.d-5.or.dabs(rfin2-xfin2).gt.1.d-5)then
     
          write(6,*)' Problem in input in radial intervals'
          write(6,*)'   no correspond to those used to generate pot'
          call flush(6)
          call MPI_BARRIER(MPI_COMM_WORLD, ierror)
          stop
      endif

!*> allocating memory
      
      npuntot=0
      do iang=1,nangu
         npuntot=max0(npuntot,npunreal(iang))         
      enddo
      ngridtot=npuntot*nangu
      nbastot=ngridtot*ncanmax
      nbastotproc=nbastot/nproc
      
        allocate(indcanproc(nbastotproc)
     &    ,indangproc(nbastotproc)
     &    ,indrgridproc(nbastotproc)
     &    ,indR1(npuntot,nangu),indR2(npuntot,nangu)
     &    ,indtotproc(npuntot,ncanprocdim,nanguprocdim,0:nproc-1)
     &    ,VVV(npuntot,nanguprocdim,nelecmax,nelecmax)
     &       , stat=ierror)

        nointegerproc_mem=nointegerproc_mem+3*nbastotproc
     &       +npuntot*nangu*2
     &    +npuntot*ncanprocdim*nanguprocdim*nproc
        imem=
     &    npuntot*nanguprocdim*nelecmax*nelecmax
        norealproc_mem=norealproc_mem+imem
      
      if(ierror.ne.0)then
        write(*,*)" error in initmem for potential in sub. pot2"
        stop
      endif
      indcanproc(:)=0
      indangproc(:)=0
      indrgridproc(:)=0
      indR1(:,:)=0
      indR2(:,:)=0
      indtotproc(:,:,:,:)=0
      VVV(:,:,:,:)=0.d0

      do ip=0,nproc-1
         ntotproc(ip)=0
         do iang_proc=1,nanguproc
            iang=indangreal(iang_proc,ip)
            do ican_proc=1,ncanproc(ip)
               ican=ibasproc(ican_proc,ip)
               do ir=1,npunreal(iang)
                  ntotproc(ip)=ntotproc(ip)+1
                  if(ntotproc(ip).le.nbastotproc)then
                     if(ip.eq.idproc)then
                        indcanproc(ntotproc(ip))=ican_proc
                        indangproc(ntotproc(ip))=iang_proc
                        indrgridproc(ntotproc(ip))=ir
                     endif
                     indtotproc(ir,ican_proc,iang_proc,ip)=ntotproc(ip)
                  endif
               enddo
            enddo
         enddo
         write(6,*)'  in processor ',ip,'  ntotproc= ',ntotproc(ip)
      enddo  ! ip = 0,nproc-1

!     reading potential
      do ip=0,nproc-1
         if(ip.ne.idproc)then
            do iang_proc=1,nanguproc
               iang=indangreal(iang_proc,ip)
               ijump=0
               write(name,'("../pot/pot.",i3.3,".dat")')iang
               open(10,file=name,status='old',err=2)
               ijump=1
 2             continue
               if(ijump.eq.0)then
                  write(name,'("pot/pot.",i3.3,".dat")')iang
                  open(10,file=name,status='old')
               endif
               
               do ir=1, npunreal(iang)
                  read(10,*)indr1(ir,iang),indr2(ir,iang)
     &          ,((vaux(ie,je),ie=1,nelecmax),je=1,nelecmax)
               enddo
               close(10)
            enddo ! iang_proc
         else
            do iang_proc=1,nanguproc
               iang=indangreal(iang_proc,ip)
               ijump=0
               write(name,'("../pot/pot.",i3.3,".dat")')iang
               open(10,file=name,status='old',err=3)
               ijump=1
 3             continue
               if(ijump.eq.0)then
                  write(name,'("pot/pot.",i3.3,".dat")')iang
                  open(10,file=name,status='old')
               endif

               do ir=1, npunreal(iang)
                  read(10,*)indr1(ir,iang),indr2(ir,iang)
     &          ,((vvv(ir,iang_proc,ie,je),ie=1,nelecmax),je=1,nelecmax)
               enddo
               close(10)
            enddo ! iang_proc
         endif
      enddo ! ip
            
!     energy interval

      vmintot=vmin
      emaxtot=vmaxtot+radcutmax+3.d0*rotcutmax
      delta2=0.5d0*(emaxtot-vmintot)
      emindlt=vmin+delta2

      write(6,*)
      write(6,'(40("-"))')
      write(6,*)'  -- Energy interval  in the calculation (eV) '

      write(6,*)'     Emin= ',vmintot/conve1/ev2cm
     &               ,'  Emax= ',emaxtot/conve1/ev2cm
      write(6,*)
      write(6,'(40("-"))')
      write(6,*)
      call flush(6)
      if( iwrt_pot == 1)then
         write(6,*)
         write(6,*)'  Write potential files to plot'
         write(6,*)
         call flush(6)
         
         call write_pot
      endif
     
      write(6,*)' ending pot2 in proc= ',idproc
     &     ,' norealproc_mem=',norealproc_mem
     &     ,' nointegerproc_mem=',nointegerproc_mem
     &     ,' memory(Gb)= '
     &   ,dble(nointegerproc_mem*4+norealproc_mem*8)*1.d-9
      call flush(6)
      return
      end subroutine pot2
!--------------------------------------------------------------

      subroutine write_pot
      use mod_gridYpara_01y2
      implicit none
      include "mpif.h"
      double precision ::  vang(npun1,nangu),vangtot(npun1,nangu)
      integer :: ifile,ie,je,jr2,ir1,ir2,iang,iang_proc,ir,nnn,ierr
      double precision :: r1,r2,ctet,angulo

      ifile=10
      do ie=1,nelec
      do je=ie,nelec
         if(idproc.eq.0.and.npun1.gt.1)then
            write(name,"('potr1r2.e',i1,'.'i1)")ie,je
            open(ifile,file=name,status='unknown')
         elseif(idproc.eq.0.and.npun1.eq.1)then
            write(name,"('potr2gam.e',i1,'.'i1)")ie,je
            open(ifile,file=name,status='unknown')
         endif

         do jr2=1,npun2
            r2=rmis2+dble(jr2-1)*ah2
!         if(r2.lt.absr2)then
            do ir1=1,npun1
            do iang=1,nangu
               vang(ir1,iang)=0.d0
               vangtot(ir1,iang)=0.d0
            enddo
            enddo
         
            do iang_proc=1,nanguproc
               iang=indangreal(iang_proc,idproc)
               do ir=1, npunreal(iang)
                  ir1=indr1(ir,iang)
                  ir2=indr2(ir,iang)
                  if(ir2.eq.jr2)then
                     vang(ir1,iang)=vvv(ir,iang_proc,ie,je)
                  endif
               enddo
            enddo

            nnn=npun1*nangu
            call MPI_REDUCE(vang,vangtot,nnn,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)

            if(idproc.eq.0.and.npun1.gt.1)then
               do ir1=1,npun1,n1plot
                  r1=rmis1+dble(ir1-1)*ah1
!                  if(r1.le.absr1)then
                     do iang=1,nangu
                        if(dabs(vangtot(ir1,iang)).lt.1.d-90)then
                           vangtot(ir1,iang)=0.d0
                        endif
                     enddo
                     write(ifile,'(500(1x,e20.12))')r1,r2
     &      ,(vangtot(ir1,iang)/conve1/ev2cm,iang=1,nangu,nangplot)
!                  endif
               enddo

               write(ifile,'()')
            elseif(idproc.eq.0.and.npun1.eq.1)then
               do ir1=1,npun1,n1plot
                  r1=rmis1+dble(ir1-1)*ah1
!                  if(r1.le.absr1)then
                     do iang=1,nangu
                        if(dabs(vangtot(ir1,iang)).lt.1.d-90)then
                           vangtot(ir1,iang)=0.d0
                        endif
                     enddo
                     do iang=1,nangu,nangplot
                        ctet=cgamma(iang)
                        angulo=dacos(ctet)*180.d0/pi
                     write(ifile,'(500(1x,e15.7))')r2,angulo
     &                       ,(vangtot(ir1,iang)/conve1/ev2cm)
                     enddo
!                  endif
               enddo

               write(ifile,'()')

            endif  ! idproc==0, npun1>1 or 1

!         endif ! r2.le.r2abs
         enddo ! ir2

         if(idproc==0)close(ifile)
      enddo  ! jele
      enddo  ! iele
  
      return
      end subroutine write_pot

!--------------------------------------------------------------
      subroutine indiproc(i,icanp,ielec,iom,iangp,ir,ir1,ir2)
      use mod_gridYpara_01y2

      implicit none

      include "mpif.h"
      integer ::i,icanp,ielec,iom,iangp,ir,ir1,ir2,ican,iang

      if(i.le.ntotproc(idproc).and.i.gt.0)then
         ir=indrgridproc(i)
         iangp=indangproc(i)
         iang=indangreal(iangp,idproc)
         icanp=indcanproc(i)
         ican=ibasproc(icanp,idproc)
         ielec=nelebas(ican)
         iom=iombas(ican)
         ir1=indr1(ir,iang)
         ir2=indr2(ir,iang)
         if(ir2.gt.npun2)write(6,*)' in indiproc ir2 > npun2, ',ir2
         if(ir1.gt.npun1)write(6,*)' in indiproc ir1 > npun1, ',ir1
      else
         write(6,*)'  i= ',i,'  > ntot= ',ntotproc(idproc)
         call flush(6)
      endif

      return
      end subroutine indiproc
!--------------------------------------------------------------
!==================================================================
      end module mod_pot_01y2
     
