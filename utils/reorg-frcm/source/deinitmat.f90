subroutine deinitmat

!===================================================================
!  Deinitializes the matrices on the grid
!===================================================================

   use pardim
   use geosol

   implicit none

   !-- deallocate arrays for matrices on the grid
   !   by calling corresponding module procedures

   call dealloc_geosol

   return

end subroutine deinitmat
