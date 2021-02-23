      program envia_transition
      implicit real *8(a-h,o-z)
      parameter(nest=21)
      character*50 stateName(nest),dirname,subname
*********************************************************
      namelist /inputgridbase/npun1,rmis1,rfin1,npun1min
     &                       ,npun2,rmis2,rfin2
     &                       ,nangu,nangplot
     &     ,Jtot,iparity,inc,nelecmax,iommin,iommax,j0
     &     ,jini,jmax,nvini,nvmax
     &                       ,nvref,jref,iomref,ielecref
*********************************************************
      namelist /inputbnd/Jtotini,iparini,nvbound,nprocbnd,maxbnddim
*********************************************************
      open(10,file='input1old.dat',status='old')
      read(10,nml = inputgridbase)
      read(10,nml = inputbnd)
      close(10)

* declaring new values

      write(*,*)' Jtotini,iparini,nvbound ?'
      read(*,*)Jtotini,iparini,nvbound
      write(*,*)'Jtot,iparity ?'
      read(*,*)Jtot,iparity
      write(*,*)' end screen input! '
      

      iommax=Jtot
      if( ((-1)**Jtot)*iparity.gt.0)then
             iommin=0
      else
             iommin=1
      endif

      open(11,file='input1.dat',status='new')      
      write(11,nml = inputgridbase)
      write(11,nml = inputbnd)
      close(11)
 
      if(iparini*iparity.ne.-1)then
          write(6,*)' initial of final partity wrong'
          write(6,*)' iparini,ipar= ',iparini,ipar
          call flush(6)
          stop
      elseif(iparini.eq.+1)then
         write(dirname,'("v",i3.3,"J",i2.2,"p-to-Jf",i2.2)')
     &            nvbound,Jtotini,Jtot
      elseif(iparini.eq.-1)then
         write(dirname,'("v",i3.3,"J",i2.2,"m-to-Jf",i2.2)')
     &            nvbound,Jtotini,Jtot
      endif

     
      
* name of electronicstates
      stateName(1)='state2Ap'
      stateName(2)='state3Ap'
      stateName(3)='state4Ap'
      stateName(4)='state5Ap'
      stateName(5)='state6Ap'
      stateName(6)='state7Ap'
      stateName(7)='state8Ap'
      stateName(8)='state9Ap'
      stateName(9)='state10Ap'
      stateName(10)='state11Ap'
      stateName(11)='state12Ap'
      stateName(12)='state13Ap'

      stateName(13)='state1As'
      stateName(14)='state2As'
      stateName(15)='state3As'
      stateName(16)='state4As'
      stateName(17)='state5As'
      stateName(18)='state6As'
      stateName(19)='state7As'
      stateName(20)='state8As'
      stateName(21)='state9As'

* preparing rotational absorption on each final state

      do ie=1,nest
         if(iparini.eq.+1)then
            write(subname,
     &        '("e",i2.2,"v",i3.3,"J",i2.2,"p-to-Jf",i2.2,".com")')
     &           ie, nvbound,Jtotini,Jtot
         elseif(iparini.eq.-1)then
            write(subname,
     &        '("e",i2.2,"v",i3.3,"J",i2.2,"m-to-Jf",i2.2,".com")')
     &            ie,nvbound,Jtotini,Jtot
         endif

         open(10,file='envia.com',status='unknown')

         write(10,*)'cd ',statename(ie)
         write(10,*)'cp ../input0.dat .'
         write(10,*)'cp ../input1.dat .'
         write(10,*)'rm input.dat -f '
         write(10,*)'cat input1.dat input0.dat >> input.dat'
         write(10,*)'mkdir ',dirname
         write(10,*)'cd ',dirname
         write(10,*)'cp ../mad3.out .'
         write(10,*)'cp ../pes* .'
         write(10,*)'cp ../input.dat .'
         write(10,*)'cp ../../submad3.com ',subname
         write(10,*)' qsub ',subname
         write(10,*)'cd ../..'
         close(10)

         call system('. ./envia.com')
      enddo ! ie=1,nest


      stop
      end
