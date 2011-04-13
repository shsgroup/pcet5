module timers

!===================================================================
!  Timer routines
!-------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2011-04-13 23:49:48 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!  Revision 5.1  2011/02/22 05:13:38  souda
!  Adding new module timers
!
!
!===================================================================

   implicit none
   private

   public :: seconde, second

!===================================================================
contains

   !-- calls system routine ETIME
   !   (works well with PGI)
   function seconde() result(sec)
      implicit none
      real(4), dimension(2) :: time
      real(4) :: etime
      real(8) :: sec
      sec = dble(etime(time))
   end function seconde

   !-- calls intrinsic routine SYSTEM_CLOCK
   !   (works well with both Intel and PGI)
   function second() result(sec)
      implicit none
      real(8) :: sec
      integer :: count, count_rate, count_max
      call system_clock(count, count_rate, count_max)
      sec = real(count,8)/real(count_rate,8)
   end function second

end module timers

