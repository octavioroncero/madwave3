! Electronic coupling routine: 
!   you have to adapt it only to generate pot.out (with colpot.sh)
!   when you want to get Electronic Predissociation halwidths
!               in the time-dependent Golden Rule approach
!  input:  
!        coordinates (in bohr) for A+BC: rBA,rBC,costetABC
!        nelecnelec,,nelecmax: the number of electronic states
!  output: 
!      coup(1:nelecmax) a vector between the initial electronic
!                      state (only one) and the nelecmax electronic states
!                      in which dynamic is  studied

      subroutine ecoupling(rBA,rBC,costetABC,coup,nelec,nelecmax)
      implicit real*8(a-h,o-z)


       
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setcoupini
      implicit real*8(a-h,o-z)
! set the magnitudes used and title of the couplings used in ecoupling routine

      write(6,"(/,40('-'),//
     & ,10x,'Coupling from ............... '
     & ,//,15x,'  , version October, 2022',//,40('-'),//)")

      return
      end


