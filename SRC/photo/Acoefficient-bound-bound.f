!!      program Aeinstein3

      implicit real*8(a-h,o-z)
      character*40 :: filegrid,filebnd
      
      include "madwave3.parameter.h"

      parameter(nv1A=1,nv2A=119,JtotA=0)
      parameter(nv1X=1,nv2X=132,JtotX=1)
      parameter(EShiftA_eV=13.054051d0)    ! in eV  

      dimension weight(nangu2),cgamma(nangu2)
      dimension wreal(nangu2)
      dimension dipbnd(npun1,npun2,nangu,nelecmax,0:Jtot)
      dimension bnd(npun1,npun2,nangu,nelecmax,0:Jtot)
      dimension EA(nv1A:nv2A),cangleA(nv1A:nv2A)
      dimension EX(nv1X:nv2X),cangleX(nv1X:nv2X)
      dimension Apartial(nv1X:nv2X)
      
*********************************************************************
*     *  Calculate Einstein coefficients between bound states
*     *     of a triatomic molecule
*     *       Using  Jacobi coordinates                                 **
*              Program works internally in a.u.                              **
*********************************************************************

      open(6,file='sal',status='unknown')
      open(7,file='sal-HCN',status='unknown')
      open(8,file='sal-HNC',status='unknown')
      
**> constants

      zero = dcmplx(0.d0,0.d0)
      zeye = dcmplx(0.d0,1.d0)
      pi = dacos(-1.d0)
      conve1 = 1.197d-4
      CONVL=.52917726D0
      CONVE=4.55633538D-6
      CONVM=.182288853D4
      hbr=convl/dsqrt(convm*conve/conve1) ! in zots
      
       ev2cm=8065.5d0
       au2eV=27.21113957d0
       zot2au=(1.d0/conve1)*conve

       cluz_au=137.d0
       epsilon0_au=0.07957747d0 

       xxx=(cluz_au)**3
       Aconstant_au= 1.d0/(3.d0*pi*xxx*epsilon0_au)    ! = 1/(3 pi hbar^4 Epsilon_0 c^3) 
       CSconstant_au= 1.d0/(cluz_au*epsilon0_au)    ! = 1/(hbar^2 Epsilon_0 c) 
       EShiftA=EshiftA_eV/au2eV
**>>   radial grids (in angstroms)

      ah2 = (rfin2-rmis2)/dble(npun2-1)
      if(npun1.gt.1)then
        div=dble(npun1-1)
        ah1 = (rfin1-rmis1)/div
        steptot=ah1*ah2
      else
        ah1=0.d0
        steptot=ah2
      endif

**>> angular grid
      
c         call gauleg(weight,cgamma,nangu2)

         do iang=1,nangu
            wreal(iang)=weight(iang)*dble(inc)
         enddo

*  Reading energies
      open(10,file='bnd/energies.dat',status='old')
      do ivX=nv1X,nv2X
         read(10,*)iii,EX(ivX),cangleX(ivX)
         Ex(ivX)=Ex(ivX)/au2eV
      enddo
      close(10)
      open(10,file='dipbnd/energies.dat',status='old')
      do ivA=nv1A,nv2A
         read(10,*)iii,EA(ivA),cangleA(ivA)
         EA(ivA)=EA(ivA)/au2eV
      enddo
      close(10)

*     Loop over initial vibrational state in the excited electyronic state
*     reading d.e | Psi_v^A >  and < Psi_v^X |

      do ivA=nv1A,nv2A
         dipbnd(:,:,:,:,:)=0.d0
         Atot=0.d0
         XCShcn=0.d0
         XCShnc=0.d0
         Ahcn=0.d0
         Ahnc=0.d0
         do iprocbnd=0,nprocdim-1
            write(filegrid,'("dipbnd/grid.ip",i3.3)')iprocbnd
            write(filebnd,'("dipbnd/dipbnd.iv",i3.3,".ip",i3.3)')
     &                       ivA,iprocbnd
            open(110,file=filebnd,status='old')
            open(111,file=filegrid,status='old')
            call readbnd(dipbnd)
            close(110)
            close(111)
         enddo  ! reading in nprocbnd files

         do ivX=nv1X,nv2X
            bnd(:,:,:,:,:)=0.d0
            do iprocbnd=0,nprocdim-1
               write(filegrid,'("bnd/grid.ip",i3.3)')iprocbnd
               write(filebnd,'("bnd/bnd.iv",i3.3,".ip",i3.3)')
     &                       ivX,iprocbnd
               open(110,file=filebnd,status='old')
               open(111,file=filegrid,status='old')
               call readbnd(bnd)
               close(110)
               close(111)
            enddo  ! reading in nprocbnd files

**>> electric dipole matrix elements

            sola=0.d0
            xmat=0.d0
            do iom=0,Jtot
            do ielec=1,nelecmax
            do iang=1,nangu
            do ir2=1,npun2
            do ir1=1,npun1
               x= bnd(ir1,ir2,iang,ielec,iom) 
               y= dipbnd(ir1,ir2,iang,ielec,iom) 
               sola=sola+x*x
               xmat=xmat+x*y
            enddo
            enddo
            enddo
            enddo
            enddo


            Ephoton=Ea(ivA)+EshiftA-Ex(ivX)
         Apartial(ivX)=xmat*xmat*Ephoton*Ephoton*Ephoton*Aconstant_au
            write(17,'(2(1x,i3),3(d15.7))')
     &           ivA,ivX,sola,xmat,Apartial(ivX)

            Atot=Atot+Apartial(ivX)
            if(cangleX(ivX).lt.0.d0)then
               Ahcn=Ahcn+Apartial(ivX)
            else
               Ahnc=Ahnc+Apartial(ivX)
            endif

            if(ivX.eq.1)then
               XCShcn=xmat*xmat*Ephoton*CSconstant_au
               EphotonHCN=Ephoton
            elseif(ivX.eq.19)then
               XCShnc=xmat*xmat*Ephoton*CSconstant_au
               EphotonHNC=Ephoton
            endif

         enddo  ! ivX
         
         write(6,'(1x,i5,20(1x,e15.7))')ivA,Atot,Ephoton,Ahcn/Atot
     &                                        ,XCShcn,xCShnc
         if(cangleA(ivA).lt.0.d0)then  ! HCN
              write(7,'(1x,i5,20(1x,e15.7))')ivA,Atot,Ahcn/Atot
     &                              ,EphotonHCN,XCShcn
         else  ! HNC
              write(8,'(1x,i5,20(1x,e15.7))')ivA,Atot,hhnc/Atot
     &                              ,EphotonHNC,XCShnc
         endif
      enddo  ! ivA

      stop
      end

*******************************************************

      subroutine readbnd(f)
      implicit real*8(a-h,o-z)
      include "madwave3.parameter.h"

      dimension f(npun1,npun2,nangu,nelecmax,0:Jtot)

      ah1=(rfin1-rmis1)/dble(npun1-1)
      ah2=(rfin2-rmis2)/dble(npun2-1)
      read(111,*)nbnddim

         read(111,*)rmis1bnd,rfin1bnd,npun1bnd
         ah1bnd=(rfin1bnd-rmis1bnd)/dble(npun1bnd-1)
         if(dabs(rmis1-rmis1bnd).gt.1.d-4
     &      .or.dabs(ah1-ah1bnd).gt.1.d-4)then
               write(6,*)'  grid in r1 for bnd non equal '
               write(6,*)'  rmis1bnd,ah1bnd= '
     &                    ,rmis1bnd,ah1bnd
               stop
         endif

         read(111,*)rmis2bnd,rfin2bnd,npun2bnd
         ah2bnd=(rfin2bnd-rmis2bnd)/dble(npun2bnd-1)
         if(dabs(rmis2-rmis2bnd).gt.1.d-4
     &      .or.dabs(ah2-ah2bnd).gt.1.d-4)then
               write(6,*)'  grid in r2 for bnd non equal '
               write(6,*)'  rmis2bnd,ah2bnd= '
     &                    ,rmis2bnd,ah2bnd
               stop
         endif

         read(111,*)nangubnd,incbnd
         if(nangu.ne.nangubnd.or.inc.ne.incbnd)then
            write(6,*)' grig in angle for bnd non equal '
            write(6,*)'  nangubnd,incbnd= ',nangubnd,incbnd
            stop
         endif

      do ibnd=1,nbnddim
            read(111,*)iii,iiicanp,ielecbnd,iombnd,iangp
     &            ,iiir,ir1bnd,ir2bnd,iangbnd
            
            read(110,*)bndvec

            f(ir1bnd,ir2bnd,iangbnd,ielecbnd,iombnd)=bndvec
      enddo

      return
      end
