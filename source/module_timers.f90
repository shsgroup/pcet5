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
!===================================================================

   implicit none
   private

   public :: seconde, seconds, second

!===================================================================
contains

   !-- calls system routine ETIME
   !   (works well with PGI)
   function seconde() result(sec)
      implicit none
      real(kind=4), dimension(2) :: time
      real(kind=4) :: etime
      real(kind=8) :: sec
      sec = real(etime(time),8)
   end function seconde

   !-- calls intrinsic routine SYSTEM_CLOCK
   !   (works well with both Intel and PGI)
   function seconds() result(sec)
      implicit none
      real(kind=8) :: sec
      integer :: count, count_rate, count_max
      call system_clock(count, count_rate, count_max)
      sec = real(count,8)/real(count_rate,8)
   end function seconds

   !-- calls intrinsic routine CPU_TIME
   !   (works well with both Intel and PGI)
   function second() result(sec)
      implicit none
      real(kind=8) :: sec
      real(kind=4) :: start
      call cpu_time(start)
      sec = real(start,8)
   end function second


   subroutine timestamp(ichannel)
   !*****************************************************************************80
   !  TIMESTAMP prints the current YMDHMS date as a time stamp.
   !
   !  Example:
   !    31 May 2001   9:45:54.872 AM
   !
   !  Licensing:
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !    18 May 2013
   !
   !  Author:
   !    John Burkardt
   !
   !  Parameters:
   !    ichannel - output channel (standard output if omitted)
   !
   !*****************************************************************************80
      implicit none

      integer(kind=4), optional, intent(in) :: ichannel

      character(len=8) :: ampm
      character (len=9), parameter, dimension(12) :: month = (/ &
       'January  ', 'February ', 'March    ', 'April    ', &
       'May      ', 'June     ', 'July     ', 'August   ', &
       'September', 'October  ', 'November ', 'December ' /)

      integer(kind=4) :: d
      integer(kind=4) :: h
      integer(kind=4) :: m
      integer(kind=4) :: mm
      integer(kind=4) :: n
      integer(kind=4) :: s
      integer(kind=4) :: values(8)
      integer(kind=4) :: y

      call date_and_time (values=values)

      y = values(1)
      m = values(2)
      d = values(3)
      h = values(5)
      n = values(6)
      s = values(7)
      mm = values(8)

      if ( h < 12 ) then
         ampm = 'AM'
      elseif ( h == 12 ) then
         if ( n == 0 .and. s == 0 ) then
            ampm = 'Noon'
         else
            ampm = 'PM'
         endif
      else
         h = h - 12
         if ( h < 12 ) then
            ampm = 'PM'
         elseif ( h == 12 ) then
            if ( n == 0 .and. s == 0 ) then
               ampm = 'Midnight'
            else
               ampm = 'AM'
            endif
         endif
      endif

      if (present(ichannel)) then
         write (ichannel, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
         & d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
      else
         write (*, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
         & d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
      endif

   end subroutine timestamp


end module timers

