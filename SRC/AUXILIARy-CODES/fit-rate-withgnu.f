      implicit real*8(a-h,o-z)
      character*40 name
      parameter(ifit=0)   ! to fit parameters with gnuplot or not)
      parameter(inc=2,nvref=2,jref=0)
      parameter(nTemp=200)
      parameter(nv0=0,nv1=1,jini=0,jmax=30,nelecmax=1)
      dimension ratej(nTemp,jini:jmax),temp(nTemp)
      dimension ratejfit(nTemp,jini:jmax)
      dimension aparam(jini:jmax),bparam(jini:jmax),cparam(jini:jmax)
     &         ,dparam(jini:jmax)

* electronic partition function as a function of temperature (x)

      Qelec(x)=1.d0   !  (2.d0/(2.d0+4.d0*exp(-64.d0/(x*0.629))))

* fit function
*
*      K(T)= c T^b  exp(d/T-a T)
*
* para ver los ficheros 
*   for i in plo*ps ; do echo $i ; gv $i > /dev/null 2>/dev/null ; done


      do ielec=1,nelecmax
      do iv=nv0,nv1
         
         write(name,"('rate.e'i2.2,'.v',i2.2)")ielec,iv
         open(10,file=name,status='old')
         do iT=1,nTemp
            read(10,*)temp(it),(ratej(iT,j),j=jini,jmax)
         enddo
         close(10)

* converting to 10^{-11} cm^3 /s : 6.126
* symmetry factor (inc=1 or 2 depending on the symmetry of reactant BC molecules
* initial rotational factor (2 jref +1)
* electronic partition factor Qe

         convau2cm3s=6.126

         do j=jini,jmax
         do iT=1,nTemp
            t=temp(it)
            Qe=Qelec(T)
            fact=convau2cm3s*Qe*dble(inc)/dble(2*jref+1)
            ratej(iT,j)=ratej(iT,j)*fact
            if(ratej(iT,j).lt.1.d-40)ratej(iT,j)=0.d0
         enddo
         enddo

         write(name,"('rate.cm3s.v'i2.2,'.j',i2.2,'.vp',i2.2)")
     &                     nvref,jref,iv
         open(10,file=name,status='new')
         do iT=1,nTemp
            write(10,'(50(1x,e15.7))')temp(it),(ratej(iT,j),j=jini,jmax)
         enddo
         close(10)
         
         write(name,"('newparameters.v'i2.2,'.j',i2.2,'.jp',i2.2)")
     &                    nvref,jref,iv
         open(10,file=name,status='unknown')

         if(ifit.eq.1)then
            do j=jini,jmax
               open(7,file='pp',status='unknown')
               do iT=1,nTemp
               write(7,*)temp(it),ratej(iT,j)
               enddo
               close(7)
            
               call system('gnuplot gnu.pru')
               open(11,file='qq',status='unknown')
               read(11,*)c,b,a,d
               read(11,*,err=1234)error1,error2
               close(11)

               write(name,'("mv pp.ps plot.v",i2.2,".j"i2.2,".ps")')iv,j
               open(12,file='copia.sh',status='unknown')
               write(12,*)name
               call system('. copia.sh')
               close(12)

               if(error2.lt.1.d-2)then
                  write(10,'(1x,i3,1x,4("&"1x,e15.7,1x),"//")')j,c,b,a,d
                  aparam(j)=a
                  bparam(j)=b
                  cparam(j)=c
                  dparam(j)=d

                  t=temp(it)
                  ratejfit(iT,j)= c*(t**b)+dexp(d/t-a*t)
               endif

 1234          continue
            
            enddo
            close(10)


            write(name,"('ratefit.cm3s.v',i2.2,'.j',i2.2,'vp',i2.2)")
     &               nvref,jref,iv
            open(10,file=name,status='new')
            do iT=1,nTemp
         write(10,'(50(1x,e15.7))')temp(it),(ratejfit(iT,j),j=jini,jmax)
            enddo
            close(10)
        endif

      enddo
      enddo

      stop
      end
   
      
