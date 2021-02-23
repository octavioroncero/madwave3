

!     heads to be used with miQCT

      subroutine setpoten
      implicit real*8(a-h,o-z)

      write(6,*)' H3 B K M P 2     P E S'


      return
      end
      
      subroutine poten(rbc,rca,rab,e,der)
      implicit real*8(a-h,o-z)
      dimension der(3)

      rref=1.4d0
      rinf=100.d0
      
         eref=4.71791677d0 /27.2113966413d0
c     id=0
c     call vh3(rref,rinf,rinf,eref,der,id)

      id=1
      call vh3(rbc,rca,rab,e,der,id)
      e=e+eref
      return
      end



!     heads to be used with MadWave3
      
      subroutine setxbcpotele(iomdiat,iomatom,sigdiat,sigatom
     &                       ,nelec,nelecmax) 
      implicit real*8(a-h,o-z)
      dimension iomdiat(nelecmax),iomatom(nelecmax)
      dimension sigdiat(nelecmax),sigatom(nelecmax)

      call setxbcpot

      if(nelecmax.ne.1)then
           write(6,'("  This PES is prepared for a single state ")')
           stop
      endif
      nelec=1
      iomdiat(1)=0
      iomatom(1)=0
      sigdiat(1)=+1
      sigatom(1)=+1

      return
      end

      subroutine potelebond(r1,r2,costet,potmat,nelec,nelecmax)
      implicit real*8(a-h,o-z)
      dimension potmat(nelecmax,nelecmax)
      dimension der(3)
      rref=1.4d0
      rasi=100.d0
      id=0
      call  vh3(rref,rasi,rasi,pot,der,id)
      vref=pot

      r3=r1*r1+r2*r2-2.d0*r1*r2*costet
      r3=dsqrt(r3)
      rr1=r1
      rr2=r2
      rr3=r3
      if(rr1.lt.0.5d0)rr1=0.5d0
      if(rr2.lt.0.5d0)rr2=0.5d0
      if(rr3.lt.0.5d0)rr3=0.5d0
      id=0
      call vh3(rr1,rr2,rr3,pot,der,id)
      potmat(1,1)=pot-vref

      return
      end

      subroutine potelexbc(rpeq,rgran,costet,potmat,nelec,nelecmax)
      implicit real*8(a-h,o-z)
      dimension potmat(nelecmax,nelecmax)

c      call potxbc(rpeq,Rgran,costet,pot)
      potmat(1,1)=pot

      return
      end
       
******************************************  SETXBCPOT  ********************
      subroutine setxbcpot
      implicit real*8(a-h,o-z)

      write(6,"(/,40('-'),//
     & ,10x,'H_3  (B K M P 2     P E S) '
     & ,//,15x,'',//,40('-'),//)")

      return
      end
********************************************************

      subroutine potxbc(rpeq,R,ctet,pot)
      implicit real*8(a-h,o-z)
      dimension der(3)

**  calculatates the potential for H_3^+( 1 ^1 A') in reactant Jacobi coordinates
** distances and energy given in a.u.
*          H
*           ^
*            \
*             \  R
*              \
*          tet  \
*              ( \
*    H  <----------------- H
*            r
      iopt=0

**    ref energy

      rref=1.4d0
      rasi=100.d0
c      call  potv(pot,rref,rasi,ctet)
      vref=pot

**    Jacobi coordinates transformation
c      xm0 = 15.99491463d0
c      xm1 = 1.00782503512d0
c      gam0=xm0/(xm0+xm1)
c      gam1=xm1/(xm0+xm1)
c      y=ctet
c      r12=rpeq
c      r13=r*r+gam1*gam1*rpeq*rpeq+2*rpeq*gam1*r*y
c      r13=dsqrt(r13)
c      r23=r*r+gam0*gam0*rpeq*rpeq-2*rpeq*gam0*r*y
c      r23=dsqrt(r23)
      

c      if(r12.gt.0.1d0.and.r13.gt.0.1d0.and.r23.gt.0.1d0)then
c         call potv(pot,rpeq,r,ctet)
c      else
c         pot=1.d10
c      endif
      pot=pot-vref

      return
      end
