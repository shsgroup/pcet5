subroutine deinitmat

!===================================================================
!  Deinitializes the matrices on the grid
!===================================================================

   use control
   use pardim
   use geosol

   implicit none

   !-- deallocate arrays for matrices on the grid
   !   by calling corresponding module procedures

   call dealloc_geosol
   call dealloc_charge

end subroutine deinitmat
