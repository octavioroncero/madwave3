      subroutine setxbcpotele(iomdiat,iomatom,sigdiat,sigatom
     &                       ,nelec,nelecmax)
      implicit real*8(a-h,o-z)
      dimension iomdiat(nelecmax),iomatom(nelecmax)
      dimension sigdiat(nelecmax),sigatom(nelecmax)

      write(6,"(/,40('-'),//
     & ,10x,'2x2 symetrized from 3x3 PES of H3+ singlets) '
     & ,//,15x,' by Aguado,Roncero Sanz, PCCP 2020',//,40('-'),//)")

      if(nelecmax.ne.3)then
           write(6,'("  This PES is prepared for 3 states ")')
           stop
      endif
      nelec=nelecmax
      iomdiat(1)=0
      iomatom(1)=0
      sigdiat(1)=+1
      sigatom(1)=+1
      iomdiat(2)=0
      iomatom(2)=0
      sigdiat(2)=+1
      sigatom(2)=+1
      iomdiat(3)=0
      iomatom(3)=0
      sigdiat(3)=+1
      sigatom(3)=+1

      return
      end

********************************************************
********************************************************

      subroutine potelebond(r10,r02,costet,potmat,nelec,nelecmax)
      use dimmatrix
      implicit real*8(a-h,o-z)
      parameter(nndim=3,nstates=3)
      dimension potmat(nelecmax,nelecmax)
      dimension potaux(nndim,nndim)
      dimension B(nNDIM,nNDIM),T(nNDIM,nNDIM),EIGEN(nNDIM)

      dimension S(nndim,nndim),aux(nndim,nndim),potsim(nndim,nndim)

      dimension der(3)

* bond coordinates  1-0-2 for:   10 + 2 --> 1 + 02

       vcutmax=10d0/27.211d0
**    ref energy
      rref=1.4d0
      rasi=100.d0
       call poth3ps(rref,rasi,rasi,potaux)
      vref=potaux(1,1)

      r12=r10*r10+r02*r02-2.d0*r10*r02*costet
      r12=dsqrt(r12)

!      if(r10.lt.0.1d0.or.r02.lt.0.1d0.or.r12.lt.0.1d0)then
!         potaux(:,:)=0.d0
!         do i=1,nndim
!            potaux(i,i)=vcutmax
!         enddo
!      else
          call poth3ps(r10,r02,r12,potaux)
!     endif

      do i=1,ndim
         potaux(i,i)=potaux(i,i)-vref
      enddo

!---------> transformation to symmetry adapted basis

      S(:,:)=0.d0
      S(1,1)=1.d0

      S(2,2)=1.d0/dsqrt(2.d0)
      S(3,2)=-1.d0/dsqrt(2.d0)

      S(2,3)=1.d0/dsqrt(2.d0)
      S(3,3)=1.d0/dsqrt(2.d0)

      do i=1,nndim
      do j=1,nndim
         aux(i,j)=0.d0
         pp=0.d0
         do k=1,nndim
            pp=pp+s(i,k)*potaux(k,j)
         enddo
         aux(i,j)=pp
      enddo
      enddo 
      do i=1,nndim
      do j=1,nndim
         pp=0.d0
         do k=1,nndim
            pp=pp+aux(i,k)*s(j,k)
         enddo
         potsim(i,j)=pp
      enddo
      enddo 
 

!         write(69,'(20(1x,f15.10))')r10,r02,costet
!     &    , potaux(1,2) , potaux(1,3) , potaux(2,3)
!     &    , potsim(1,2) , potsim(1,3) , potsim(2,3)!

!         call flush(69)

      do ie=1,nelec
      do je=1,nelec
         potmat(ie,je)=potsim(ie,je)
      enddo
      enddo

      return
      end

