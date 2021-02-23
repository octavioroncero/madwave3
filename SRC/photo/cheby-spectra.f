       program chebyspectra
       use ieee_arithmetic
       implicit real*8(a-h,o-y)
       implicit complex*16(z)
       real*8, allocatable :: auto(:)
!
!  absorption spectra from autocorrelation function
!    calculated with a real modified Chebyshev propagator
!       in the unit 5 --> auto
!    output in unit 6  --> spectrum  in atomic units of energy and photodissociation cros section
!
       
********************************************************* Energies in eV
      namelist /inputEspectro/Emin,Emax,ne,xkgam,kmax
*********************************************************Only Jtot is needed
      namelist /inputgridbase/npun1,rmis1,rfin1,npun1min
     &                       ,npun2,rmis2,rfin2
     &                       ,nangu,nangplot
     &     ,Jtot,iparity,inc,nelecmax,iommin,iommax,j0
     &     ,jini,jmax,nvini,nvmax
     &      ,nvref,jref,iomref,ielecref
********************************************************
      open(10,file='input.dat',status='old')
      read(10,nml = inputEspectro)
      close(10)
      allocate(auto(0:kmax))
*energy grid in eV
*       parameter(Emin=-5.d0,Emax=5.d0,ne=100000,xkgam=0.000005)

* kmax (maxinum number of Chebyshev iterations read
*       parameter(kmax=500000)
*     dimension auto(0:kmax)

       facrot=1.d0/3.d0

* conversion constants
       conve1 = 1.197d-4
       CONVL=.52917726D0
       CONVE=4.55633538D-6
       CONVM=.182288853D4
       ev2cm=8065.5d0
       au2eV=27.21113957d0
       zot2au=(1.d0/conve1)*conve
       hbr=convl/dsqrt(convm*conve/conve1)
       pi=dacos(-1.d0)

       cluz_au=137.d0
       epsilon0_au=0.07957747d0 

       CSconstant_au= 1.d0/(cluz_au*epsilon0_au)    ! = 1/(hbar^2 Epsilon_0 c) 
 
* reference energy = Eshift
       open(10,file='../func/eref',status='old')
        read(10,*)eref
       close(10)
       Eshift=eref/(conve1*ev2cm) 
* reading autocorrelation function
       read(5,*)xnorm,ebndini,emindlt,delta2
       auto(0)=1.d0
       do k=1,kmax
          read(5,*,end=1)kk,auto(kk)
          kmaxreal=kk
       enddo
 1     continue

       Estep=(eMAX-eMIN)/dble(ne-1)

* calculating spectrum by modified Chebyshev transform of autocorrelation function 

       enorm=0.d0
       do ie=1,ne
          EE=Emin+dble(ie-1)*Estep
          E=EE*ev2cm*conve1
          Es=(E-emindlt)/delta2
          deno=1.d0-Es*Es
          deno=dsqrt(deno)
          
          espectro=auto(0)*hbr
          factor=2.d0*hbr
          do k=1,kmaxreal
             expo=-dble(k)*dacos(Es)
             rexp=dcos(expo)
             Ckcheby=factor*rexp     *dexp(-dble(k)*xkgam)
             espectro=espectro+Ckcheby*auto(k)
          enddo
          espectro=espectro/(deno*delta2)
c         if(ieee_is_nan(espectro))espectro=0.d0
          if(isnan(espectro))espectro=0.d0
 
          ephoton_au=(ee+eshift-ebndini)/au2eV
          espectro_au=espectro*xnorm/(pi*hbr*zot2au)   ! in units of 1/Energy
          enorm=enorm+espectro_au/xnorm
          Xsection_au=CSconstant_au*ephoton_au*espectro_au
          write(6,'(20(1x,e15.7))')ee,espectro,Xsection_au,ephoton_au
     &                          ,espectro_au,CSconstant_au
       enddo
       write(6,*)'  integral= ',enorm*Estep/au2eV

       stop
       end
