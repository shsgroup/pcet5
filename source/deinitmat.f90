subroutine deinitmat
!===================================================================
!  Deinitializes the matrices on the grid
!===================================================================
!
!  souda
!  2010/06/25 20:02:35
!  4.1
!  Exp
!  deinitmat.f90,v 4.1 2010/06/25 20:02:35 souda Exp
!  deinitmat.f90,v
!  Revision 4.1  2010/06/25 20:02:35  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2004/01/13 19:32:10  souda
!  Initial PCET-4.0 Release
!
!
!-------------------------------------------------------------------
   use pardim
   use gasmat
   use solmat
   use geogas
   use geosol
   use quantum

   implicit none

   ! deallocate arrays for matrices on the grid
   ! by calling corresponding module procedures

   call dealloc_gasmat
   call dealloc_solmat
   call dealloc_geogas
   call dealloc_geosol
   call dealloc_pquant
   call dealloc_gquant

   return

end subroutine deinitmat
