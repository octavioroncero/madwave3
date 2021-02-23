!!      program Aeinstein3

      implicit real*8(a-h,o-z)
      character*40 :: filegrid,filebnd,name,system
      
      include "madwave3.parameter.h"

      parameter(nv1A=1,nv2A=103,JtotA=0)
      parameter(nv1X=1,nv2X=11,JtotX=0)
      parameter(EShiftA_eV=(93.2820599d0-92.81154726)*27.211d0)    ! in eV  para el HCN(A')
!      parameter(EShiftA_eV=(93.2820599d0-92.78133357)*27.211d0)    ! in eV  para el HCN(A'')

      dimension weight(nangu2),cgamma(nangu2)
      dimension wreal(nangu2)
      dimension bndA(npun1,npun2,nangu,nelecmax,0:Jtot)
      dimension bndX(npun1,npun2,nangu,nelecmax,0:Jtot)
      dimension EA(nv1A:nv2A),cangleA(nv1A:nv2A)
      dimension EX(nv1X:nv2X),cangleX(nv1X:nv2X)
      dimension Apartial(nv1X:nv2X)
      
**********************************************************************************************
*     *  Calculatethe overlap  between bound states                                          *
*     *              belonging to different electronic states X and A                        *
*     *     of a triatomic molecule                                                          *
*     *       Using  Jacobi coordinates                                                      *
*              Program works internally in a.u.                                              *
**********************************************************************************************
*     *    Ej. HCN photoionization

      open(6,file='sal',status='unknown')
      
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
      open(10,file='bndX/energias.dat',status='old')
      do ivX=nv1X,nv2X
         read(10,*)iii,EX(ivX),cangleX(ivX)
!         Ex(ivX)=Ex(ivX)/au2eV
         Ex(ivX)=Ex(ivX)*conve
      enddo
      close(10)
      open(10,file='bndA/energias.dat',status='old')
      do ivA=nv1A,nv2A
         read(10,*)iii,EA(ivA),cangleA(ivA)
         EA(ivA)=EA(ivA)/au2eV
      enddo
      close(10)

*     Loop over initial vibrational state in the excited electyronic state
*     reading  | Psi_v^A >  and < Psi_v^X |
      do ivX=nv1X,nv2X
         bndX(:,:,:,:,:)=0.d0
         write(name,'("Abs.vX",i3.3)')ivx
         open(17,file=name,status='unknown')
         do iprocbnd=0,nprocdim-1
            write(filegrid,'("bndX/grid.ip",i3.3)')iprocbnd
            write(filebnd,'("bndX/bnd.iv",i3.3,".ip",i3.3)')
     &                       ivX,iprocbnd
            open(110,file=filebnd,status='old')
            open(111,file=filegrid,status='old')
            call readbnd(bndX)
            close(110)
            close(111)
         enddo  ! reading in nprocbnd files

         do ivA=nv1A,nv2A
            bndA(:,:,:,:,:)=0.d0
            Atot=0.d0
            XCShcn=0.d0
            XCShnc=0.d0
            Ahcn=0.d0
            Ahnc=0.d0
            do iprocbnd=0,nprocdim-1
               write(filegrid,'("bndA/grid.ip",i3.3)')iprocbnd
               write(filebnd,'("bndA/bnd.iv",i3.3,".ip",i3.3)')
     &                       ivA,iprocbnd
               open(110,file=filebnd,status='old')
               open(111,file=filegrid,status='old')
               call readbnd(bndA)
               close(110)
               close(111)
            enddo  ! reading in nprocbnd files


**>> overlaps, neglecting the contribution from the electronic part

            sola=0.d0
            xmat=0.d0
            do iom=0,Jtot
            do ielec=1,nelecmax
            do iang=1,nangu
            do ir2=1,npun2
            do ir1=1,npun1
               x= bndX(ir1,ir2,iang,ielec,iom) 
               y= bndA(ir1,ir2,iang,ielec,iom) 
               sola=sola+x*x
               xmat=xmat+x*y
            enddo
            enddo
            enddo
            enddo
            enddo


            Ephoton=Ea(ivA)+EshiftA-Ex(ivX)
            Apartial(ivX)=xmat*xmat*Ephoton*CSconstant_au
            write(17,'(2(e15.7),2(1x,i3),3(e15.7))')
     &           Ephoton,Apartial(ivX),ivA,ivX,sola,xmat,Apartial(ivX)


         enddo  ! ivA
         close(17)
      enddo  ! ivX

      stop
      end

*******************************************************

      subroutine readbnd(f)
      implicit real*8(a-h,o-z)
      character*40 system
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
