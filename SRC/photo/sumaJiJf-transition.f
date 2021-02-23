      program suma_transitionJiJf
      implicit real *8(a-h,o-z)
      parameter(nest=20)
      character*50 stateName(nest),dirname,filename
      character*50 outname
      parameter(nener=100000,einieV=7.d0,efineV=13.d0)
      parameter(neread=100000)
      dimension Eshift(nest)
      dimension x(neread),f(neread,2)
      dimension espetot(nener),espepar(nener,nest)

*********************************************************
      namelist /inputgridbase/npun1,rmis1,rfin1,npun1min
     &                       ,npun2,rmis2,rfin2
     &                       ,nangu,nangplot
     &     ,Jtot,iparity,inc,nelecmax,iommin,iommax,j0
     &     ,jini,jmax,nvini,nvmax
     &                       ,nvref,jref,iomref,ielecref
*********************************************************
      namelist /inputbnd/igridbnd
     &                  ,Jtotini,iparini,nvbound,nprocbnd,maxbnddim
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
         write(outname,'("espe.v",i3.3,"J",i2.2,"p-to-Jf",i2.2)')
     &            nvbound,Jtotini,Jtot
      elseif(iparini.eq.-1)then
         write(dirname,'("v",i3.3,"J",i2.2,"m-to-Jf",i2.2)')
     &            nvbound,Jtotini,Jtot
         write(outname,'("espe.v",i3.3,"J",i2.2,"m-to-Jf",i2.2)')
     &            nvbound,Jtotini,Jtot
      endif
     
      Cau2eV=27.211d0

* initializing total espectrum

      espetot(:)=0.d0
      espepar(:,:)=0.d0
      estep=(efineV-einieV)/dble(nener-1)
            
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

      stateName(12)='state1As'
      stateName(13)='state2As'
      stateName(14)='state3As'
      stateName(15)='state4As'
      stateName(16)='state5As'
      stateName(17)='state6As'
      stateName(18)='state7As'
      stateName(19)='state8As'
      stateName(20)='state9As'

* preparing rotational absorption on each final state

      do iest=1,nest
         if(statename(iest)(9:9).eq.' ') then
            filename=statename(iest)(1:8)//"/"//dirname(1:16)//"/espe"
            write(6,*)filename
         else
            filename=statename(iest)(1:9)//"/"//dirname(1:16)//"/espe"
            write(6,*)filename
         endif
         call flush(6)
         open(10,file=filename,status='old')
         ireal=0
         x(:)=0.d0
         f(:,:)=0.d0
         ireal=0
         do ieread=1,neread
            ireal=ireal+1
            read(10,*,end=1)a,c,xcs,e
            if(xcs.lt.0.d0)xcs=0.d0
            x(ieread)=e*Cau2eV
            f(ieread,1)=xcs   ! in a.u.
         enddo
 1       continue
         close(10)
         
         call SPLSET(F,X,ireal,neread)

         iold=2
         do ie=1,nener
            e=einieV+dble(ie-1)*estep
            espe=0.d0
            if(e.ge.x(1).and.e.le.x(ireal))then
                 iold=2
                 espe=SPLINQ(F,X,IOLD,ireal,e,neread)
            endif

            if(espe.lt.1.d-90)espe=0.d0
            espetot(ie)=espetot(ie)+espe
            espepar(ie,iest)=espe
         enddo  ! ie=1,nener

         
      enddo ! iest=1,nest

* printig

      open(6,file=outname,status='unknown')
      do ie=1,nener
         e=einieV+dble(ie-1)*estep
         write(6,'(50(1x,e15.7))')e,espetot(ie)
     &           ,(espepar(ie,iest),iest=1,nest)
      enddo

      stop
      end
