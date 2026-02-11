      program thermal_spectra

      implicit none
      double precision, parameter :: Temperature_K=150.d0
      integer,parameter :: Jtotmax=10,nvmax=1,ne_photon=50001
      double precision,parameter :: Emin_photon_eV=0.5
      double precision,parameter :: Emax_photon_eV=3.d0
      double precision,parameter :: pi=dacos(-1.d0)
      integer :: Jtot,ipar,iv,iiv,Jdelta,Jfinal,ipar_final,ie
      character*60 :: name
      real* 8 :: Energy,total_partition
      real*8,allocatable :: Ebound(:,:,:),partition(:,:,:)
      real*8,allocatable :: sigma(:,:,:,:),Ephoton(:)
      real*8,allocatable :: total_sigma(:,:)
      real*8 :: Ezot,Ecm,Temperature_eV,par,x1,x2,e0,sig_read
      real*8 :: Emin_photon,Emax_photon,deltaE_photon,eee
      integer :: iread_first,iphoton_read

      Temperature_eV=Temperature_K*8.61738569d-5
      Emin_photon=Emin_photon_eV/27.211
      Emax_photon=Emax_photon_eV/27.211
      deltaE_photon=( Emax_photon- Emin_photon)/dble(ne_photon+1)

      allocate( Ebound(0:Jtotmax,nvmax,-1:1)
     &         ,partition(0:Jtotmax,nvmax,-1:1)
     &         ,sigma(ne_photon,0:Jtotmax,nvmax,-1:1)
     &     ,Ephoton(ne_photon)
     &     ,total_sigma(ne_photon,-1:1)
     &     )

      Ebound(:,:,:)=0.d0
      partition(:,:,:)=0.d0
      sigma(:,:,:,:)=0.d0
      Ephoton(:)=0.d0
      total_sigma(:,:)=0.d0
      do ie=1,ne_photon
         Ephoton(ie)= Emin_photon+dble(ie-1)*deltaE_photon
      enddo
      
!     Reading initial bound states energies

      write(6,*)' Reading initial bound states (in eV) '
      e0=1d10
      do Jtot=0,Jtotmax
         ipar=(-1)**Jtot
         if(ipar.gt.0)then
            write(name,"('bndJ',i2.2,'p/niveles.dat')")Jtot
         else
            write(name,"('bndJ',i2.2,'m/niveles.dat')")Jtot
         end if
         write(6,*)' reading file= ',name
         open(10,file=name,status='old')
         do iv=1,nvmax

            read(10,*)iiv,Ezot,Ecm,Ebound(Jtot,iv,ipar)
            e0=dmin1(e0,Ebound(Jtot,iv,ipar))
            
         end do                 ! iv
         close(10)

      end do                    ! Jtot
      
      write(6,*)' lower energy (eV) = ',e0

!     Calculating Parition function

      
      total_partition=0.d0
      do Jtot=0,Jtotmax
         ipar=(-1)**Jtot
         do iv=1,nvmax
            par = dble(2*Jtot+1)
     &           * dexp(-(Ebound(Jtot,iv,ipar)-E0)/Temperature_eV)
            partition(Jtot,iv,ipar)=par
            total_partition=total_partition+par
         end do ! iv

      end do                    ! Jtot
      
!     renormalizing to 1

      write(6,*)' Jtot,ipar,iv, weight '
      do Jtot=0,Jtotmax
         ipar=(-1)**Jtot
         do iv=1,nvmax
            partition(Jtot,iv,ipar) = partition(Jtot,iv,ipar)
     &           / total_partition
            write(6,*)Jtot,ipar,iv,partition(Jtot,iv,ipar)
         end do ! iv

      end do  ! Jtot

!     reading spectra

      iread_first=0
      do Jtot=0,Jtotmax
         ipar=(-1)**Jtot
         do iv=1,nvmax
            do Jdelta=-1,1
               Jfinal=Jtot+Jdelta
               ipar_final=-ipar
               if(Jfinal.ge.0)then
                  if(Jfinal.eq.0.and.ipar_final.lt.0)go to 10

                  if(ipar.gt.0)then
                     write(name,"('J',i3.3,'m-fromJ',i2.2,'p-v',i1
     &                    ,'/espe')")Jfinal,Jtot,iv
                  else if(ipar.lt.0)then
                     write(name,"('J',i3.3,'p-fromJ',i2.2,'m-v',i1
     &                    ,'/espe')")Jfinal,Jtot,iv
                  end if
                  write(6,*)' reading ',name
                  open(10,file=name,status='old')

                  do ie=1,ne_photon
                     read(10,*)x1,sig_read,x2,eee
                     iphoton_read=(eee-Emin_photon)/deltaE_photon+1

                     if(iphoton_read.ge.1.and
     &                        .iphoton_read.le.ne_photon)then
      
                        total_sigma(iphoton_read,Jdelta)=
     &                        total_sigma(iphoton_read,Jdelta) +sig_read
     &                         *dble(Jfinal*2+1)*partition(Jtot,iv,ipar)
                      end if
                  end do

                  close(10)

 10            continue
               end if
           end do
         end do
      end do

      open(10,file='thermal_spectra.dat',status='unknown')
      write(10,*)'# photon Energy (hartree)  p q r spectra '
      do ie=1,ne_photon
         
         write(10,*)Ephoton(ie),(total_sigma(ie,Jdelta),Jdelta=-1,1)

      end do
      
      stop
      end
