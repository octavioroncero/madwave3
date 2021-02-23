***********************  Lanzval  ***********************************

      subroutine lanzvalpariter(flanz,Hu,eigen,ntot,nbastot
     &                    ,etrue,ntrue,
     &                    nloop,ntimes,nkryini,criter,conve,
     &                    idproc,nproc) 
      parameter(nnnnmax=100000)
      implicit real*8(a-h,o-z)
      character*50 name
      include "mpif.h"

*     ************************************************
*     *            Evaluation of Eigenvalues         *
*     *  using non-orthogonal Lanczos' procedure     *
*     * of Cullum & Willoughby, "Lanczos Algorithms  *
*     *      for large Eigenvalue Computations",     *
*     *            Birkh"auser, Boston,1985          *
*     ************************************************
*     *   Input                                      *
*     *        flanz: double precision matrix        *
*     *               of dimensions (nbastot,0:2)    *
*     *               for lanczos iteration.         *
*     *        Hu,eigen: double precision matrices   *
*     *                  of dimensions (nbastot)     *
*     *        ntot: no. of basis functions          *
*     *        nbastot: dimension .ge.ntot           *
*     *        criter: maximum error for convergency *
*     *                 of eigenvalues.              *
*     ************************************************
*     *   Output                                     *
*     *        etrue: converged eigenvalues          *
*     *        ntrue: no. of converged eigenvalues   *
*     ************************************************
*     *   External subroutines                       *
*     *         - initia: provides initial guess     *
*     *                   in flanz(i,2)              *
*     *         - Hmatrix: evaluates H \phi          *
*     *                input--> eigen(i) = \phi      *
*     *                output--> Hu(i) = H \phi      *
*     *      with i=1,..,ntot                        *
*     *                  with a dimension nbastot    *
*     *                                              *
*     *   uses also VALNKRY subroutine               *
*     *      to calculate eigenvalues for different  *
*     *      nkry values for convergency check.      *
*     ************************************************
*     * O. Roncero, Madrid (Spain)                   *
*     *     Last modification, December 2010         *
*     ************************************************

      dimension alpha(nnnnmax),beta(nnnnmax)
      dimension alpori(nnnnmax),betori(nnnnmax)
      dimension flanz(nbastot,0:2),Hu(nbastot),eigen(nbastot)
      dimension etrue(nloop*ntimes)

      nkrymax=nloop*ntimes
      if(nkrymax.gt.nnnnmax)then
          write(6,*)' !!! Attention in lanzvec !!!!'
          write(6,*)' ***>> nkrymax >> ',nnnnmax
          stop
      endif
     
* Initial guess and second krilov function
      if(nkryini.le.1)then
         nkryini=1

         call initia_lanzvec(flanz,ntot,nbastot)
 
         do i=1,ntot
            Hu(i)=0.d0
            eigen(i)=flanz(i,2)
         enddo

         call Hmatrix(eigen,Hu,ntot,nbastot)
 
         alpha(1)=0.d0
         do i=1,ntot
           alpha(1)=alpha(1)+Hu(i)*flanz(i,2)
         enddo

         a1=alpha(1)
         call MPI_REDUCE(a1,a2,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr) 
         call MPI_BCAST(a2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

         alpha(1)=a2
c      write(6,*)'  alpha1 = ',a2 
         bet=0.d0
         do i=1,ntot
            flanz(i,1)=Hu(i)-alpha(1)*flanz(i,2)
             bet=bet+flanz(i,1)*flanz(i,1)
         enddo

         call MPI_REDUCE(bet,bet2,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(bet2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

         beta(1)=dsqrt(bet2)
c         write(6,*)'  beta1 = ',beta(1)
         do i=1,ntot
            flanz(i,1)=flanz(i,1)/beta(1)
         enddo
 
         do i=1,ntot
            Hu(i)=0.d0
            eigen(i)=flanz(i,1)
         enddo

         call Hmatrix(eigen,Hu,ntot,nbastot)
 
         alpha(2)=0.d0
         do i=1,ntot
            alpha(2)=alpha(2)+Hu(i)*flanz(i,1)
     &           -beta(1)*flanz(i,2)*flanz(i,1)
         enddo

         a1=alpha(2)
         a2=0.d0
         call MPI_REDUCE(a1,a2,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(a2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

         alpha(2)=a2

         ndelta=2
      else ! nkryini.eq.0
         ndelta=1
c        open(10,file='Lanczos.mat',status='unknown')
c        do i=1,nkryini
c           read(10,*)alpha(i),beta(i)
c        enddo
c        close(10)

c        write(name,'("Lanczos.fun.ip",i3.3)')idproc
c        open(10,file=name,status='unknown')
c        do i=1,ntot
c           read(10,*)flanz(i,2),flanz(i,1)
c        enddo
c        close(10)

         do i=1,ntot
            Hu(i)=0.d0
            eigen(i)=flanz(i,1)
         enddo

         call Hmatrix(eigen,Hu,ntot,nbastot)

      endif

c     write(6,*)'  initial guess taken in Lanczos recursion'
c     do i=nkryini,nkryini+ndelta
c        write(6,*)i,alpha(i),beta(i)
c     enddo
                                                                                 
*       Loop over the whole Krilov subspace

 1000 continue

         do ilanc=nkryini+ndelta,nkryini+ntimes
            bet=0.d0
            do i=1,ntot
               flanz(i,0)=Hu(i)-alpha(ilanc-1)*flanz(i,1)
     &                -beta(ilanc-2)*flanz(i,2)
               bet=bet+flanz(i,0)*flanz(i,0)
            enddo

            bet2=0.d0
            call MPI_REDUCE(bet,bet2,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(bet2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

            beta(ilanc-1)=dsqrt(bet2)
            do i=1,ntot
               flanz(i,0)=flanz(i,0)/beta(ilanc-1)
            enddo
            do i=1,ntot
               Hu(i)=0.d0
               eigen(i)=flanz(i,0)
            enddo

            call Hmatrix(eigen,Hu,ntot,nbastot)

            alpha(ilanc)=0.d0
            do i=1,ntot
               alpha(ilanc)=alpha(ilanc)+Hu(i)*flanz(i,0)
     &                       -beta(ilanc-1)*flanz(i,1)*flanz(i,0)
            enddo
 
            a1=alpha(ilanc)
            a2=0.d0
            call MPI_REDUCE(a1,a2,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(a2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        
            alpha(ilanc)=a2

            do i=1,ntot
               flanz(i,2)=flanz(i,1)
               flanz(i,1)=flanz(i,0)
               flanz(i,0)=0.d0
            enddo
         enddo  ! ilanc

         nkry=nkryini+ntimes
         nkryini=nkry
**>> beta (nkry)

         bet=0.d0
         do i=1,ntot
            xx=Hu(i)-alpha(nkry)*flanz(i,1)
     &                 -beta(nkry-1)*flanz(i,2)
            bet=bet+xx*xx
         enddo

          bet2=0.d0
         call MPI_REDUCE(bet,bet2,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(bet2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

         beta(nkry)=dsqrt(bet2)
 
         do ikry=1,nkry
            alpori(ikry)=alpha(ikry)
            betori(ikry)=beta(ikry)   
         enddo

** eigenvalues 

         if(idproc.eq.0)then
            open(10,file='Lanczos.mat',status='unknown')
            do i=1,nkry
               write(10,*)alpha(i),beta(i)
            enddo
            close(10)

         endif

         write(name,'("Lanczos.fun.ip",i3.3)')idproc
         open(10,file=name,status='unknown')
         do i=1,ntot
            write(10,*)flanz(i,2),flanz(i,1)
         enddo
         close(10)

         call valnkry(alpori,betori,etrue,nkry,ntrue,nkrymax)


         if(idproc.eq.0)then
            write(6,'(//,10x,"--> Converged eigenvalues: "
     &       ,2x,"no. val",2x,i5,2x
     &        ,"nkry= ",2x,1(10x,i5),//)')
     &        ntrue,nkry

            write(name,'("Lanczos.E.nkry",i6.6)')nkry
            open(10,file=name,status='unknown')
            write(10,*)ntrue
            do i=1,ntrue
            write(10,*)i,etrue(i)
            enddo
            close(10)
         endif

         ndelta=1
         if(ilanc.lt.nkrymax)goto 1000

         if(idproc.eq.0)then
         open(3,file='cont.data',status='unknown')
         write(3,'(1x,i6,1x,i6,"   ! ilanc, nkrymax")')nkrymax,nkrymax
         write(3,'(1x,i6, "   ! ntrue")')ntrue
         do i=1,ntrue
            write(3,*)i,etrue(i)
         enddo
         write(3,'(" 0  ! ivec")')
         close(3)
         endif
      return
      end
*********************** cgradvec  **********************************

      subroutine cgradvecpar(ediag,flanz,Hu,eigen,ntot,nbastot
     &                   ,nomxiter,criter,conve,idproc,nproc)

      implicit real*8(a-h,o-z)
      dimension flanz(nbastot,0:2),Hu(nbastot),eigen(nbastot)
      include "mpif.h"

*     ************************************************************
*     *   eigenvectors, using the conjugate gradient method      *
*     *        v_i = v_(i-1) + lambda_(i-1) p_(i-1)              *
*     *                      r_i = A v_i                         *
*     *  eps_(i-1)= < r_i | A p_(i-1) > / < p_(i-1)| A p_(i-1)>  *
*     *              p_i = -r_i + eps_(i-1) p_(i-1)              *
*     *                 lambda_i = -  < r_i|p_i>                 *
*     *                   v_i equiv  flanz( , 0)                 *
*     *                   r_i equiv  flanz( , 1)                 *
*     *                   p_i equiv  flanz( , 2)                 *
*     ************************************************************
                                                  
* Initialization

      call initia_lanzvec(flanz,ntot,nbastot)

*  v_0
      do i=1,ntot
         flanz(i,0)=flanz(i,2)
         flanz(i,2)=0.d0
      enddo
 
* r_0 and p_0
      do i=1,ntot
         Hu(i)=0.d0
         eigen(i)=flanz(i,0)
      enddo

      call Hmatrix(eigen,Hu,ntot,nbastot)

      do i=1,ntot
         Hu(i)=Hu(i)-ediag*flanz(i,0)
         flanz(i,1)=Hu(i)
         flanz(i,2)=-flanz(i,1)
      enddo
* lambda_0:  Hu equi A p_0
      do i=1,ntot
         Hu(i)=0.d0
         eigen(i)=flanz(i,2)
      enddo
  
      call Hmatrix(eigen,Hu,ntot,nbastot)

      xlamnum=0.d0
      xlamden=0.d0
      xlamnumtot=0.d0
      xlamdentot=0.d0
      do i=1,ntot
         Hu(i)=Hu(i)-ediag*flanz(i,2)
         xlamnum=xlamnum+flanz(i,1)*flanz(i,2)
         xlamden=xlamden+flanz(i,2)*Hu(i)
      enddo

      call MPI_REDUCE(xlamnum,xlamnumtot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xlamnumtot,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      call MPI_REDUCE(xlamden,xlamdentot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xlamdentot,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)


         xlam=-xlamnumtot/xlamdentot

      iter=0
      ifail=0
 1928 continue
         iter=iter+1                                                         

* v_i : flanz( ,0)
         xnormproc=0.d0
         xnorm=0.d0
         do i=1,ntot
            eigen(i)=flanz(i,0)
            flanz(i,0)=flanz(i,0)+xlam*flanz(i,2)
            xnormproc=xnormproc+flanz(i,0)*flanz(i,0)
         enddo

         call MPI_REDUCE(xnormproc,xnorm,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(xnorm,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

         xnorm=dsqrt(xnorm)
         do i=1,ntot
            flanz(i,0)=flanz(i,0)/xnorm
         enddo
* r_i : flanz( ,1)
         do i=1,ntot
            Hu(i)=0.d0
            eigen(i)=flanz(i,0)
         enddo

         call Hmatrix(eigen,Hu,ntot,nbastot)

         emeanproc=0.d0
         emean=0.d0
         do i=1,ntot
            emeanproc=emeanproc+Hu(i)*flanz(i,0)
            Hu(i)=Hu(i)-ediag*flanz(i,0)
            flanz(i,1)=Hu(i)
         enddo

         call MPI_REDUCE(emeanproc,emean,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(emean,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  
* epsilon_(i-1)
         do i=1,ntot
            Hu(i)=0.d0
            eigen(i)=flanz(i,2)
         enddo

         call Hmatrix(eigen,Hu,ntot,nbastot)

         epsnum=0.d0
         epsden=0.d0
         epsnumtot=0.d0   
         epsdentot=0.d0   
         do i=1,ntot
            Hu(i)=Hu(i)-ediag*flanz(i,2)
            epsnum=epsnum+flanz(i,1)*Hu(i)
            epsden=epsden+flanz(i,2)*Hu(i)
         enddo

         call MPI_REDUCE(epsnum,epsnumtot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(epsnumtot,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(epsden,epsdentot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(epsdentot,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
 
         epsilon=epsnumtot/epsdentot

* p_i : flanz( ,2)
         do i=1,ntot
            flanz(i,2)=epsilon*flanz(i,2)-flanz(i,1)
         enddo
* lambda
         do i=1,ntot
           Hu(i)=0.d0
           eigen(i)=flanz(i,2)
         enddo

         call Hmatrix(eigen,Hu,ntot,nbastot)

         xlamnum=0.d0
         xlamden=0.d0
         xlamnumtot=0.d0
         xlamdentot=0.d0
         do i=1,ntot
            Hu(i)=Hu(i)-ediag*flanz(i,2)
            xlamnum=xlamnum+flanz(i,1)*flanz(i,2)
            xlamden=xlamden+flanz(i,2)*Hu(i)
         enddo

         call MPI_REDUCE(xlamnum,xlamnumtot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(xlamnumtot,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(xlamden,xlamdentot,1,MPI_REAL8,MPI_SUM
     &                             ,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(xlamdentot,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

         xlam=-xlamnumtot/xlamdentot
c         write(6,"(' iter= ',i5,'  emean,ediag= ',4(1x,e15.7))")
c    &                     iter,emean/conve,ediag/conve
c    &                    ,dabs(emean-ediag),criter
      if(iter.le.nomxiter)then
         if(dabs(emean-ediag).gt.criter)go to 1928
      else
         ifail=ifail+1
      endif
      do i=1,ntot
         eigen(i)=flanz(i,0)
      enddo
      if(idproc.eq.0)then
         write(6,'(9x,"   After ",i5," iterations the energy is")')iter
         write(6,'(9x,"   E= ",e15.7,"  while it should be ",e15.7)')
     &      emean/conve,ediag/conve
      endif                                            
                                            
      return
      end

***************************** valnkry *******************************

      subroutine valnkry(alpha,beta,etrue,nkry,ntrue,nkrymax)
      implicit real*8(a-h,o-z)
      parameter(nnnnmax=100000,ntrmax=20)

*     ************************************
*     *     Evaluation of Eigenvalues    *
*     *      using Lanczos' procedure    *
*     ************************************

      dimension etrue(nkrymax),alpha(nnnnmax),beta(nnnnmax)
*auxiliar matrices

      dimension indori(nnnnmax),indorir(nnnnmax)
      dimension en2(nnnnmax)
      dimension en2r(nnnnmax),alphar(nnnnmax),betar(nnnnmax)
      dimension ereal(nnnnmax),degen(nnnnmax)


      do ikry=1,nkry
         alphar(ikry)=alpha(ikry)
         betar(ikry)=beta(ikry)
      enddo
      betar(1)=0.d0      

*      Diagonalization of the tridiagonal Hamiltonian matrix
 
      call tqli(alpha,beta,nkry,nnnnmax)
      call tqli(alphar,betar,nkry,nnnnmax)
 
**> ordering eigenvalues
 
      do ikry=1,nkry
         en2(ikry)=alpha(ikry)
         en2r(ikry)=alphar(ikry)
         indori(ikry)=ikry
         indorir(ikry)=ikry
      enddo
 
 11   continue
      nchan=0
      do ikry=1,nkry-1
         j=indori(ikry)
         j1=indori(ikry+1)
         if(en2(j).gt.en2(j1))then
            indori(ikry)=j1
            indori(ikry+1)=j
            nchan=nchan+1
         endif
      enddo
      if(nchan.gt.0)go to 11
 
 12   continue
      nchan=0
      do ikry=1,nkry-1
         j=indorir(ikry)
         j1=indorir(ikry+1)
         if(en2r(j).gt.en2r(j1))then
            indorir(ikry)=j1
            indorir(ikry+1)=j
            nchan=nchan+1
         endif
      enddo
      if(nchan.gt.0)go to 12
                                
***>> eliminating ghost states
 
      ireal=1
      ereal(1)=en2(indori(1))
      degen(1)=1
      critener=1.d-8
 
      do ikry=2,nkry
         e=en2(indori(ikry))
         if(dabs(ereal(ireal)-e).lt.critener)then
            degen(ireal)=degen(ireal)+1.d0
         else
            ireal=ireal+1
            degen(ireal)=1.d0
            ereal(ireal)=e
         endif
      enddo
 
      itrue=0
      do i=1,ireal
         if(degen(i).gt.1.d0)then
            itrue=itrue+1
            etrue(itrue)=ereal(i)
         else
            imatch=0
            do ikry=1,nkry
               if(dabs(ereal(i)-en2r(ikry)).lt.critener)then
                 imatch=imatch+1
               endif
            enddo
            if(imatch.eq.0)then
               itrue=itrue+1
               etrue(itrue)=ereal(i)
            endif
         endif
      enddo
      ntrue=itrue
                                          
      return
      end
************************  sincbas  *************************************

      subroutine sincbasnew(rmis,rmax,npun,xmred,hbr
     &                   ,rhod2dvr,rhod1dvr,rhom2dvr)
      implicit real*8(a-h,o-z)

      dimension rhod2dvr(npun,npun),rhod1dvr(npun,npun),rhom2dvr(npun)
      pi=dacos(-1.d0)

*
*     evaluation of Kinetic matrix and hbar/2 mu rho^2
*       for a radial variable using a SINC-DVR representation
*
*           rhod2dvr: -{\hbar^2/ 2 mu} {d^2 / d \rho^2}
*           rhod1dvr:                  {d / d \rho}
*           rhom2dvr:  {\hbar^2/ 2 mu} {1 / \rho^2}
*
*       in a.u.
*
      ah=(rmax-rmis)/dble(npun-1)

**>> -hbar^2/2mu  d^2/d\rho^2

      do 100 iv=1,npun
         do 90 jv=iv,npun

            menos=(jv-iv)
            mas=iv+jv
            sign=dble( (-1)**(menos))
            if(iv.eq.jv)then
               x1=pi*pi/6.d0
               x2=1.d0/( dble(mas)*dble(mas) )
               cint=x1+x2
            else
               x1=1.d0/( dble(menos)*dble(menos) )
               x2=1.d0/( dble(mas)*dble(mas) )
               cint=x1+x2
            endif
            cint=sign*cint/(ah*ah)

            rhod2dvr(iv,jv)=hbr*hbr*cint/xmred
            rhod2dvr(jv,iv)=rhod2dvr(iv,jv)

 90      continue
 100  continue

**>>   d/d\rho

      do iv=1,npun
         do jv=iv,npun

            yj=ah*dble(iv-1)*pi/dble(Rmax-Rmis)
            yk=ah*dble(jv-1)*pi/dble(Rmax-Rmis)
            if(iv.eq.jv)then
               cint=0.d0
            else
               cint=0.d0
               do n=1,npun
               do m=1,npun
                  mas=n+m
                  menos=n-m
                  sign=dble( (-1)**(mas) -1 )
                  x1=dsin(dble(n)*yj)
                  x2=dsin(dble(m)*yk)

                  s1=0.d0
                  s2=0.d0
                  if(mas.ne.0)s1=dble(m)/dble(mas)                  
                  if(menos.ne.0)s2=dble(m)/dble(menos)

                  cint=cint+sign*x1*x2*(s1+s2)
               enddo
               enddo
               cint=-cint*2.d0*ah/((rmax-rmis)*(rmax-rmis))
            endif

            rhod1dvr(iv,jv)=cint
            rhod1dvr(jv,iv)=-rhod1dvr(iv,jv)

         enddo
      enddo

**>>  calculating <1/rho^2>

      do irho=1,npun
         rho=rmis+dble(irho-1)*ah
         rave=hbr*hbr*0.5d0/(rho*rho*xmred)
         rhom2dvr(irho)=rave
      enddo

      return
      end

