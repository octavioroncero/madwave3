      subroutine setxbcpotele(iomdiat,iomatom,sigdiat,sigatom
     &                       ,nelec,nelecmax)
      implicit real*8(a-h,o-z)
      dimension iomdiat(nelecmax),iomatom(nelecmax)
      dimension sigdiat(nelecmax),sigatom(nelecmax)

      call setxbcpot
      return
      end

      subroutine potelebond(rhAhB,rhBAu,costet,potmat,nelec,nelecmax)
      implicit real*8(a-h,o-z)
      dimension potmat(nelecmax,nelecmax)
      dimension der(3)

      return
      end
      subroutine setxbcpot
      implicit real*8(a-h,o-z)

      write(6,"(/,40('-'),//
     & ,10x,'general potential routines '
     & ,//,15x,'  , to compilae madwave3 for all systems only once'
     & ,//,40('-'),//)")

      return
      end

 
