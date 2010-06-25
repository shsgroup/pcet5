function second() result(tic)
   implicit none
   real(4), dimension(2) :: time
   real(4) :: etime
   real(8) :: tic
   tic = dble(etime(time))
   !tic = 0.d0
   return
end function second
