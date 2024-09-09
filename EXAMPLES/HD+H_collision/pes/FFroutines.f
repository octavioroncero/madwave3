!   A,B,C,D  --> points in cartesian coordinates
!-----------------------------------------------
! Functions :
!   DISTANCE(A,B,derA,derB) --> distance between A and B 
!                               and derivative in cartesian (convention (A-B))
!   ANGLE(A,B,C)  -->  angle between AB and BC
!   DIHEDRAL(A,B,C,D) --> dihedral angle
!   dipinter(dip1,dip2,vecCM)  --> dipole interaction between dip1 and dip2 at relative position vectCM
!                                   and derivatives respect to  dip1, dip2,vectCM 
!   DELTAKRO(i,) --> Kronecker Delta
!   mydamp       --> tanh based damping function
!-----------------------------------------------
! Rountines
! getVector(A,B,AB)  --> returns vector AB (B-A)     
! ddotprod(A,B,prod) --> prod = A . B 
! vvecprod(A,B,C)    --> C = A ^ B
! atommass(At,mass)  --> returns atomic mass of atom At
! dipolepointcharge  --> return dipole vector calculated from given point-charges distribution
! DipDipPC(pA,pB,natA,natB,atlA,atlB,qA,qB)  --> giving position and charge distribution of molecule A and B, returns
!                                              dipole interaction and derivatives respect to atom positions
! VECPROJ(A,B,C)     --> C is the vector A rotated on vector B axis
! 
!-----------------------------------------------------
      subroutine getVECTOR(A,B,AB)
       implicit none
       integer i
       real*8 A(3), B(3), AB(3)
!
       do i=1,3
         AB(i)=B(i)-A(i)
       enddo
       return
      end subroutine getVECTOR
!-----------------------------------------------------
      real*8 function DISTANCE(A,B,derA,derB)
! convention for derivatives : A - B
       implicit none
       integer i
       real*8 A(3), B(3), derA(3), derB(3)
! 
       distance=0.d0
       do i= 1,3
        distance= distance+ (A(i)-B(i))*(A(i)-B(i))
       enddo
        distance= dsqrt(distance)

! derviatives with respect to  cartesian coordinates
       do i=1,3
        derA(i)= (A(i)-B(i))/distance
        derB(i)= -derA(i)
       enddo
       return
      end function
!-----------------------------------------------------
      real*8 function ANGLE(A,B,C,derA,derB,derC)
       implicit none
       integer i
       real*8 A(3), B(3), C(3), BA(3), BC(3),derA(3), derB(3), derC(3)
       real*8 derdA(3),derdB1(3),derdB2(3),derdC(3)
       real*8 distance,dab,dcb,prod1,prod0,calpha,calpha2
!        
       call getVECTOR(B,A,BA)
       call getVECTOR(B,C,BC)
       call ddotprod(BA,BC,prod1)
       dab = DISTANCE(A,B,derdA,derdB1)       
       dcb = DISTANCE(C,B,derdC,derdB2)       
       prod0= dab * dcb
       if(dabs(prod0).lt.1.d-10) then
          calpha = 1.d0
       else
          calpha = prod1/prod0
          if(calpha.gt.1.d0)  calpha=1.d0
          if(calpha.lt.-1.d0) calpha=-1.d0
       endif      
       calpha2 = calpha*calpha
       angle= dacos(calpha)
! derviatives with respect to  cartesian coordinates
       if((1.d0-dabs(calpha)).lt.1.d-10) then
        do i=1,3
         derA(i)=0.d0
         derB(i)=0.d0
         derC(i)=0.d0
        enddo
       else
        do i=1,3
         derA(i)=(calpha*BA(i)/dab-BC(i)/dcb)/(dab*dsqrt(1.d0-calpha2)) 
         derC(i)=(calpha*BC(i)/dcb-BA(i)/dab)/(dcb*dsqrt(1.d0-calpha2))
         derB(i)= -derA(i)-derC(i)
        enddo
       endif
       return
      end function ANGLE
!-----------------------------------------------------
      real*8 function DIHEDRAL(A,B,C,D,derA,derB,derC,derD)
! convention :
!
!       A        D       |   C-------B     
!        \      /        |    \     /
!          B---C         |     \   A
!                        |      \ 
!                        |       D 
!        normal          |    improper  
!                        |
!     call order         |   call order
!   A <-- B <-- C --> D  |  A <-- C <-- B  --> D
!
       implicit none
       integer i
       real*8 A(3), B(3), C(3), D(3)
       real*8 derA(3), derB(3), derC(3), derD(3)
       real*8 F(3), G(3), H(3), angle,distance,der(3)
       real*8 dAB, dBC, dDC,ABC,BCD,ssign 
       real*8 AA(3),BB(3), ZZ(3),prod0, prod1,normAA,normBB
!
       call getVECTOR(B,A,F)
       call getVECTOR(C,B,G)
       call getVECTOR(C,D,H)
       dAB = DISTANCE(A,B,der,der)
       dBC = DISTANCE(B,C,der,der)
       dDC = DISTANCE(D,C,der,der)
       ABC = ANGLE(A,B,C,der,der,der)
       BCD = ANGLE(B,C,D,der,der,der)
       if((1.d0-dabs(dcos(ABC))).lt.1.d-4) then
!       write(*,*) 'ohoh', dacos(ABC)
          DIHEDRAL =0.d0
          do i=1,3
           derA(i)= 0.d0
           derB(i)= 0.d0
           derC(i)= 0.d0
           derD(i)= 0.d0
          enddo 
          goto 100
        endif
       if((1.d0-dabs(dcos(BCD))).lt.1.d-4) then
!      write(*,*) 'ahah', dacos(BCD)
          DIHEDRAL =0.d0
          do i=1,3
           derA(i)= 0.d0
           derB(i)= 0.d0
           derC(i)= 0.d0
           derD(i)= 0.d0
          enddo
          goto 100
        endif

!
       call vvecprod(F,G,AA)
       call vvecprod(H,G,BB)
       call ddotprod(AA,BB,prod1) 
       call vvecprod(AA,BB,ZZ)
       call ddotprod(G,ZZ,ssign)
       ssign=-sign(1.d0,ssign)
       normAA=0.d0
       normBB=0.d0
       do i=1,3
        normAA= normAA + AA(i)*AA(i) 
        normBB= normBB + BB(i)*BB(i) 
       enddo
       normAA= dsqrt(normAA)
       normBB= dsqrt(normBB)
!       prod0 = dsin(ABC)*dAB*dBC*dBC*dDC*dsin(BCD)
       if(abs(normAA).lt.1.d-4) normAA= 1.d-4
       if(abs(normBB).lt.1.d-4) normBB= 1.d-4
       prod0 = normAA*normBB 
       DIHEDRAL = prod1/prod0 
       if (DIHEDRAL.lt.-1.d0) DIHEDRAL = -1.d0
       if (DIHEDRAL.gt.1.d0) DIHEDRAL = 1.d0
       DIHEDRAL = dacos(DIHEDRAL)
!  derviatives with respect to  cartesian coordinates
       if(dBC.lt.1.d-10) then
        do i=1,3
         derA(i)= 0.d0
         derB(i)= 0.d0
         derC(i)= 0.d0
         derD(i)= 0.d0
        enddo
       else
        call ddotprod(F,G,prod0)
        call ddotprod(H,G,prod1)
        do i=1,3
         derA(i)= -AA(i)*dBC/(normAA*normAA) 

         derD(i)=  BB(i)*dBC/(normBB*normBB) 
         derB(i)= -derA(i) + AA(i)*prod0/(normAA*normAA*dBC) 
     &            - BB(i)*prod1/(normBB*normBB*dBC)
         derC(i)= -derD(i) - AA(i)*prod0/(normAA*normAA*dBC)
     &            + BB(i)*prod1/(normBB*normBB*dBC) 
        enddo
        do i=1,3
         derA(i)= derA(i)*ssign
         derB(i)= derB(i)*ssign
         derC(i)= derC(i)*ssign
         derD(i)= derD(i)*ssign
        enddo
       endif       
 100   continue
       return
      end function DIHEDRAL
!-----------------------------------------------------
      subroutine atommass(atom,mass)
      character*1 atom
      real*8 mass
      SELECT CASE (atom)
      case ('C')
        mass = 12.d0
      case ('H')
        mass = 1.007825035d0
      case ('O')
        mass = 15.99491463d0
      END SELECT
      end subroutine
!-----------------------------------------------------
      subroutine getCM(posit,natom,atoms,CM,der)
!  atoms -- > array with list of atoms 
       implicit none
       integer i, j, natom
       character*1 atoms(natom)
       real*8 posit(3*natom), mass(natom), masstot,CM(3),der(natom,3)
       real*8 coor(3),r12,dr12(3),kk(3)
!  
       masstot=0.d0
       do i=1,3
        CM(i)=0.d0
       enddo
       do i=1,natom
        call atommass(atoms(i),mass(i))
        masstot= masstot + mass(i)
        CM(1) = CM(1) + mass(i)*posit(3*(i-1)+1)
        CM(2) = CM(2) + mass(i)*posit(3*(i-1)+2)
        CM(3) = CM(3) + mass(i)*posit(3*(i-1)+3)
       enddo
        
        CM(1) = CM(1)/masstot
        CM(2) = CM(2)/masstot
        CM(3) = CM(3)/masstot
! derivatives respect to cartesian coordinates
       do j=1,natom
        do i =1,3
         der(j,i)= mass(j)/masstot
        enddo
       enddo
       return
      end subroutine getCM
!-----------------------------------------------------
      subroutine vec_distance(vec, dis, d_vec)
      ! Calculates module of a vector. 
      ! Input:
      ! vec: Cartesian coordinates of vector which module is calculated.
      !---------------------------------------------- 
      ! Output:
      ! dis: Calculated distance.
      ! d_vec: Derivative of the module wrt each coordinate of vec.
      
      implicit none
      real(8), intent(in) :: vec(3)
      real(8), intent(out) :: dis, d_vec(3)
      
      integer :: i
      
      dis = 0e0
      d_vec = 0e0
      
      do i=1,3
         dis = dis + (vec(i) * vec(i))
      enddo
      dis = sqrt(dis)
      
      d_vec = vec / dis
      end subroutine

!-----------------------------------------------------
      subroutine derunitvector(A,B,dA,dB)
        implicit none
        integer i,j
        real*8 A(3), B(3), dA(3,3), dB(3,3), derA(3), derB(3), r,r2
        real*8 vec(3),distance,pp(3,3),deltakro
!
        r=DISTANCE(B,A,derB,derA)
        r2=r*r
        call getVECTOR(A,B,vec)
!
        do i=1,3
!         dA(1,i)= -r*deltakro(1,i)+
!     &            derA(1)*vec(1)+derA(1)*vec(2)+derA(1)*vec(3) 
!         dA(2,i)= -r*deltakro(2,i)+
!     &            derA(2)*vec(1)+derA(2)*vec(2)+derA(2)*vec(3) 
!         dA(3,i)= -r*deltakro(3,i)+
!     &            derA(3)*vec(1)+derA(3)*vec(2)+derA(3)*vec(3) 
        dA(1,i)=  -r*deltakro(1,i)-derA(i)*vec(1)
        dA(2,i)=  -r*deltakro(2,i)-derA(i)*vec(2)
        dA(3,i)=  -r*deltakro(3,i)-derA(i)*vec(3)
        dB(1,i)=  r*deltakro(1,i)-derB(i)*vec(1)
        dB(2,i)=  r*deltakro(2,i)-derB(i)*vec(2)
        dB(3,i)=  r*deltakro(3,i)-derB(i)*vec(3)
        enddo
!
         dA=dA/r2
         dB=dB/r2

         return
      end subroutine derunitvector  
!-----------------------------------------------------
      subroutine DipDipPC(pA,pB,natA,natB,atlA,atlB,qA,qB,pot,dpA,dpB)
       implicit none
       integer i, j, natA, natB
       character*2 atlA(natA), atlB(natB)
       real*8 pA(3*natA), pB(3*natB),dpA(natA*3),dpB(natB*3) 
       real*8 dpAtmp(natA,3),dpBtmp(natB,3)
       real*8 cmA(3), cmB(3), dercmA(natA,3), dercmB(natB,3)
       real*8 qA(natA), qB(natB), pot, RRab(3),dRRA(3),dRRB(3)
       real*8 dipA(3), dipB(3),ddipA(3),ddipB(3),ddipRAB(3), dipinter
       real*8 derdipA(natA,3), derdipB(natB,3)
!
       call getCM(pA,natA,atlA,cmA,dercmA)
       call getCM(pB,natB,atlB,cmB,dercmB)
       call getVECTOR(cmA,cmB,RRab)
!
       call DIPOLEpointcharge(pA,natA,atlA,qA,dipA,derdipA)
       call DIPOLEpointcharge(pB,natB,atlB,qB,dipB,derdipB)
!
       pot=dipinter(dipA,dipB,RRab,ddipA,ddipB,ddipRAB)
! derivative respect to cartersina coordinates 
! dpot/dR = dpot/ddipA*ddipA/DR + dpot/ddipB*ddipB/DR + dpot/ddipRAB*ddipRAB/DR
       do j=1,natA
        do i=1,3
         dpAtmp(j,i) = ddipA(i)*derdipA(j,i) - ddipRAB(i)  *dercmA(j,i)
         dpA(3*(j-1)+i) = dpAtmp(j,i)
        enddo
       enddo
       do j=1,natB
        do i=1,3
         dpBtmp(j,i) = ddipB(i)*derdipB(j,i) + ddipRAB(i)  *dercmB(j,i)
         dpB(3*(j-1)+i) = dpBtmp(j,i)
        enddo
       enddo
       return
      end subroutine DipDipPC
!-----------------------------------------------------
      subroutine DIPOLEpointcharge(pos,natom,atoms,charge,dip,derdip)
       implicit none
       integer i, j, natom
       character*2 atoms(natom)
       real*8 pos(3*natom),mass(natom),charge(natom),CM(3),dip(3)
       real*8 derCM(natom,3), derdip(natom,3),r12,dr12(3),kk(3)
       real*8 coor(3)
!
       call getCM(pos,natom,atoms,CM,derCM)
       do i=1,3
         dip(i)=0.d0
       enddo
       do j=1,natom
        do i=1,3
         dip(i) = dip(i) + charge(j)*(pos(3*(j-1)+i)-CM(i))
        enddo
       enddo
! derivatives of dipole respect to positions
       do j=1,natom
        do i= 1,3
         derdip(j,i)= charge(j)   
        enddo
       enddo   
       return
      end subroutine DIPOLEpointcharge
!-----------------------------------------------------
      real*8 function dipinter(dip1,dip2,vectCM,ddip1,ddip2,dvectCM) 
       implicit none
       integer i
       real*8 dip1(3),dip2(3),cm1(3),cm2(3),ddip1(3),ddip2(3)
       real*8 vectCM(3), dvectCM(3), sp1, sp2, sp3
       real*8 n1,n2,nCM, pi, eps0,pre
       real*8 DdipDsp1,DdipDsp2,DdipDsp3, DdipDncm, DncmDvectcm(3)
       real*8 Dsp1Ddip1(3),Dsp2Ddip1(3),Dsp3Ddip1(3) 
       real*8 Dsp1Ddip2(3),Dsp2Ddip2(3),Dsp3Ddip2(3)
       real*8 Dsp1DvectCM(3),Dsp2DvectCM(3),Dsp3DvectCM(3)
       pi=dacos(-1.d0)
       Eps0=0.079577471545947687d0         
       pre=4.d0*pi*eps0
!
       nCM= vectCM(1)*vectCM(1)+vectCM(2)*vectCM(2)+vectCM(3)*vectCM(3)
       nCM= dsqrt(nCM)
       n1 = dip1(1)*dip1(1)+dip1(2)*dip1(2)+dip1(3)*dip1(3)
       n1 =dsqrt(n1)
       n2 = dip2(1)*dip2(1)+dip2(2)*dip2(2)+dip2(3)*dip2(3)
       n2 =dsqrt(n2)
! 
       call ddotprod(vectCM,dip1,sp1)
       call ddotprod(vectCM,dip2,sp2)
       call ddotprod(dip1,dip2,sp3)
!        
       dipinter= pre*(sp3-3.d0*sp1*sp2/nCM**2)/nCM**3
!
!  derviatives of dipinter respect to sp1, sp2, sp3, nCM
       DdipDsp1=-3.d0*pre*sp2/nCM**5
       DdipDsp2=-3.d0*pre*sp1/nCM**5
       DdipDsp3= pre/nCM**3
       DdipDncm= (15.d0*sp1*sp2/nCM**6-3.d0*sp3/nCM**4)/pre
! derivatives of sp1, sp2, sp3 respect to dip1,dip2 and vectCM 
!
       do i=1,3
        Dsp1Ddip1(i) = vectCM(i) 
        Dsp1DvectCM(i) = dip1(i) 
!    
        Dsp2Ddip2(i)= vectCM(i)
        Dsp2DvectCM(i) = dip2(i)
!
        Dsp3Ddip1(i)= dip2(i)
        Dsp3Ddip2(i)= dip1(i)
       enddo
! derivatives of dipinter respect to dip1,dip2 and vectCM, chain rule
       do i=1,3
        ddip1(i)= DdipDsp1*Dsp1Ddip1(i) + DdipDsp3*Dsp3Ddip1(i)
        ddip2(i)= DdipDsp2*Dsp2Ddip2(i) + DdipDsp3*Dsp3Ddip2(i)
        dvectCM(i)= DdipDsp1*Dsp1DvectCM(i) + DdipDsp2*Dsp2DvectCM(i)
     &            + DdipDncm*vectCM(i)/nCM 
       enddo
       return
      end function dipinter
!-----------------------------------------------------
      real*8 function deltaKro(i,j)
       implicit none
       integer i,j

       if(i.eq.j) then
          deltaKro= 1.d0
       else
          deltaKro= 0.d0
       endif
       end function deltaKro
!-----------------------------------------------------
       subroutine VECPROJ(A,B,C)
! vector A projected (rotation) on vector B gives vector C
        implicit none
        integer i, j
        real*8 A(3), B(3), C(3), rot(3,3), skew(3,3), skew2(3,3)
        real*8 vectemp(3), A1(3), nA, B1(3), nB, cc, deltaKro
! normalize A and B
        nA=0.d0
        nB=0.d0
        skew= 0.d0
        do i=1,3
         nA= nA+ A(i)*A(i)
         nB= nB+ B(i)*B(i)
        enddo
         nA= dsqrt(nA)
         nB= dsqrt(nB)
        do i=1,3
         A1(i)= A(i)/ nA
         B1(i)= B(i)/ nB
        enddo
! obtain rotation matrix
        call vvecprod(A1,B1,vectemp)
        call ddotprod(A1,B1,cc)
        skew(1,2)= -vectemp(3)
        skew(1,3)= vectemp(2)
        skew(2,1)= vectemp(3)
        skew(2,3)= -vectemp(1)
        skew(3,1)= -vectemp(2)
        skew(3,2)= vectemp(1)
        skew2=MATMUL(skew,skew)
        do i=1,3
        do j=1,3
         if((1.d0+cc).gt.1.d-15) then
           rot(i,j)= deltaKro(i,j) + skew(i,j) + skew2(i,j)/(1.d0+cc)
         else
           rot(i,j)= -deltaKro(i,j)
         endif
        enddo
        enddo
        C = MATMUL(rot,A)
        return
       end subroutine VECPROJ 
!-----------------------------------------------------
       real*8 function mydamp(x,a,b,dx)
       implicit none
       real*8 a,b,x,dx

        mydamp=0.5d0- 0.5d0*tanh(a*(x-b))
        dx= 2.d0*cosh(a*(b-x))*cosh(a*(b-x)) 
        dx= -a/dx
       end function mydamp
!-----------------------------------------------------
      real*8 function Vmorse(x,D,r,alp,dx)
      implicit none
      real*8 x,D,r,alp,dx,kk
!      Vmorse=D*(dexp(-2.d0*alp*(x-r))-2.d0*dexp(-alp*(x-r)))
      kk = dexp(-alp*(x-r))
      Vmorse = D*(kk*kk-2.d0*kk)
      dx = -2.d0*D*(kk-1.d0)*alp*kk
      return
      end function Vmorse 
!-----------------------------------------------------
      subroutine ddotprod(a,b,prod)
      implicit real*8(a-h,o-z)
      dimension a(3),b(3)

!* scalar product:  c= a . b

      prod=0.d0
      do i=1,3
         prod=prod+a(i)*b(i)
      enddo

      return
      end subroutine ddotprod
!-----------------------------------------------------
      subroutine vvecprod(a,b,c)
      implicit real*8(a-h,o-z)
      dimension a(3),b(3),c(3)

!        vector product:  c= a x b

      c(1)=a(2)*b(3)-b(2)*a(3)
      c(2)=a(3)*b(1)-b(3)*a(1)
      c(3)=a(1)*b(2)-b(1)*a(2)

      return
      end subroutine vvecprod

      subroutine scale_dis(n_inp, a_exp, min_x, d, sd, dsd)
      ! Exponential scaling of distances.
      implicit none
      integer :: i, n_inp
      real(8), dimension(n_inp) :: d, sd, dsd, a_exp, min_x
      real(8) :: aux_exp, diff
      
      do i=1,n_inp
         if (d(i) .ge. min_x(i)) then
            diff = d(i) - min_x(i)
            aux_exp = dexp(- a_exp(i) * diff * diff)
            sd(i) = 1d0 - aux_exp
            dsd(i) = 2d0 * a_exp(i) * diff * aux_exp
         else
            sd(i) = 0d0
            dsd(i) = 0d0
         endif
      enddo
      end subroutine
      !-----------------------------------------------------
      subroutine scale_dis_lin(n_inp, max_x, min_x, d, sd, dsd)
      ! Linear scaling of distances.
      implicit none
      integer :: n_inp
      real(8), dimension(n_inp) :: d, sd, dsd, max_x, min_x, diff
      
      diff = max_x - min_x
      
      dsd = 2d0 / diff
      sd = dsd * (d - min_x) - 1d0
      end subroutine
      !-----------------------------------------------------
      subroutine scale_dis_sig(n_inp, mean, std, d, sd, dsd)
      ! Sigmoid scaling of distances.
      implicit none
      integer :: i, n_inp
      real(8), dimension(n_inp) :: d, sd, dsd, mean, std
      real(8) :: aux, sg
      
      do i=1,n_inp
            aux = 2d0 / std(i)
            sd(i) = sg(aux * (d(i) - mean(i))) 
            dsd(i) = sd(i) * (1d0 - sd(i)) 
            dsd(i) = aux * dsd(i) 
      enddo
      end subroutine
      !-----------------------------------------------------
      subroutine step(n_inp, n_hid, input, output, w, b, der,activation)
      ! Input:
      !n_inp: size of input vector
      !n_hid: number of hidden neurons in layer.
      !activation: Logical. True if sigmoid activation is performed.
      !input: input of the perceptron.
      !w: weights of the perceptron.
      !b: bias of the perceptron.
      
      !Output:
      !output: output vector
      !der: derivative of output wrt input.
      implicit none
      integer, intent(in) :: n_inp, n_hid
      real(8), intent(in) :: input(n_inp), w(n_inp, n_hid), b(n_hid)
      logical, intent(in) :: activation
      real(8), intent(out) :: output(n_hid), der(n_inp, n_hid)
      real(8) :: auxder, sg
      integer :: i, j
      
      output = b
      der = w
      do i=1,n_hid
         do j=1, n_inp
            output(i) = output(i) + input(j) * w(j,i)
         enddo
         if (activation) then
            output(i) = 1e0 / (1e0 + exp(-output(i)))
            auxder = output(i) * (1d0 - output(i))
            der(:,i) = der(:,i) * auxder
         endif
      enddo
      end subroutine
      !-----------------------------------------------------
      subroutine distance_vec(vec, dis, d_vec)
      ! Calculates module of a vector. 
      ! Input:
      ! vec: Cartesian coordinates of vector which module is calculated.
      !---------------------------------------------- 
      ! Output:
      ! dis: Calculated distance.
      ! d_vec: Derivative of the module wrt each coordinate of vec.
      
      implicit none
      real(8), intent(in) :: vec(3)
      real(8), intent(out) :: dis, d_vec(3)
      
      integer :: i
      
      dis = 0d0
      d_vec = 0d0
      
      do i=1,3
         dis = dis + (vec(i) * vec(i))
      enddo
      dis = dsqrt(dis)
      
      d_vec = vec / dis
      end subroutine
      !-----------------------------------------------------
      real(8) function gauss(x, mean, sd)
      implicit none
      real(8), intent(in) :: x, mean, sd
      
      gauss = (x - mean) / sd
      gauss = gauss * gauss
      gauss = - 5d-1 * gauss
      gauss = dexp(gauss)
      end function
      !-----------------------------------------------------

      subroutine dipAB(vec1, vec2, q, a, dip, d_dip1, d_dip2)
      ! Calculates the contribution of a bond to the total dipole:
      ! vec = vec1 - vec2
      ! dip = (q + a / |vec|) * vec
      ! Input:
      ! vec*: Position of atoms for which the dipole contribution is calculated
      ! q, a: Fitted parameters to reproduce the dipole.
      !---------------------------------------------- 
      ! Output:
      ! dip: Dipole contribution of the bond.
      ! d_dip*: d(dip)_x/d(vec*)_y with x=(x,y,z) and y=(x,y,z)
      
      implicit none
      real(8), intent(in) :: vec1(3), vec2(3), q, a
      real(8), intent(out) :: dip(3), d_dip1(3,3), d_dip2(3,3)
      
      integer :: i, j
      real(8) :: vec(3), d_vec(3), vec_dis
      real(8) :: aux, vec_dis2
      
      dip = 0d0
      d_dip1 = 0d0
      d_dip2 = 0d0
      
      ! Calculate vector and distance of the bond.
      vec = vec1 - vec2
      call distance_vec(vec, vec_dis, d_vec)
      
      ! Calculate dipole moment.
      aux = q + (a / vec_dis)
      dip = aux * vec
      
      ! Square of bond distance.
      vec_dis2 = vec_dis * vec_dis
      
      ! Derivative of inverse of distance with respect to each coordinate of vec.
      ! a * d(1/|r|)/d(r_a) = - a * (d(|r|)/d(r_a)) / (|r|)**2
      d_vec = - a * d_vec / vec_dis2
      
      ! Calculate derivative of dipole moment wrt each coordinate of "vec" vector:
      d_dip1 = 0d0
      do i=1,3
         do j=i,3
            if (i.eq.j) then
               d_dip1(i,j) = aux
            endif
            d_dip1(i,j) = d_dip1(i,j) + vec(i) * d_vec(j) 
            d_dip1(j,i) = d_dip1(i,j)
         enddo
      enddo
      
      d_dip2 = - d_dip1
      end subroutine

      function sg(x)
      implicit none
      real(8) :: x, sg
      
      sg = 1 + dexp(-x)
      sg = 1 / sg
      end function

!-----------------------------------------------------

      subroutine jacdiagon2(Hmat,n,ndim,T,eigen)
!     diagonalize Hmat  --> providing eigen (eigenvalues) and T (eigenvectors)
!            which are ordered in increasing energy
      implicit none
      integer n,ndim,nrot
      real*8 :: Hmat(ndim,ndim),eigen(ndim),T(ndim,ndim)

      call jacobi2(Hmat,n,ndim,eigen,T,nrot)

      call eigsrt2(eigen,T,n,ndim)

      return
      end
      
!***************************************************
      subroutine jacobi2(a,n,np,d,v,nrot)
      implicit none
      integer,parameter :: NMAX=50000
      INTEGER :: n,np,nrot
      REAL*8 :: a(np,np),d(np),v(np,np)
      INTEGER :: i,ip,iq,j
      REAL*8 :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.d0
11      continue
        v(ip,ip)=1.d0
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.d0
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.d0
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+dabs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.d0)return
        if(i.lt.4)then
          tresh=0.2d0*sm/n**2
        else
          tresh=0.d0
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.d0*dabs(a(ip,iq))
            if((i.gt.4).and.(dabs(d(ip))+
     &  g.eq.dabs(d(ip))).and.(dabs(d(iq))+g.eq.dabs(d(iq))))then
              a(ip,iq)=0.d0
            else if(dabs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(dabs(h)+g.eq.dabs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5d0*h/a(ip,iq)
                t=1.d0/(dabs(theta)+dsqrt(1.d0+theta**2))
                if(theta.lt.0.d0)t=-t
              endif
              c=1.d0/dsqrt(1+t**2)
              s=t*c
              tau=s/(1.d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.d0
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.d0
23      continue
24    continue
      stop 'too many iterations in jacdiagon'
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software *1n#!-013.
!*********************************************************************
      subroutine eigsrt2(d,v,n,np)
!ordering eigenvalues and eigenvectors
      implicit none
      INTEGER:: n,np
      REAL*8 :: d(np),v(np,np)
      INTEGER :: i,j,k
      REAL*8 :: p

      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).le.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software *1n#!-013.
!-----------------------------------------------------
