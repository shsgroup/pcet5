module timers

!===================================================================
!  Timer routines
!-------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2011-02-22 05:13:38 $
!  $Revision: 5.1 $
!  $Log: not supported by cvs2svn $
!
!===================================================================

   implicit none
   private

   public :: second

!===================================================================
contains

   function second() result(tic)
      implicit none
      real(4), dimension(2) :: time
      real(4) :: etime
      real(8) :: tic
      tic = dble(etime(time))
      !tic = 0.d0
   end function second

end module timers

