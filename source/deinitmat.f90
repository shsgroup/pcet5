subroutine deinitmat

!===================================================================
!  Deinitializes the matrices on the grid
!===================================================================
!
!  $Author: souda $
!  $Date: 2011-02-20 00:58:11 $
!  $Revision: 5.3 $
!  $Log: not supported by cvs2svn $
!  Revision 5.2  2010/10/28 21:29:35  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!-------------------------------------------------------------------

   use pardim
   use gasmat
   use solmat
   use geogas
   use geosol
   use quantum
   use feszz_3d, only: reset_feszz3_counter

   implicit none

   ! deallocate arrays for matrices on the grid
   ! by calling corresponding module procedures

   call dealloc_gasmat
   call dealloc_solmat
   call dealloc_geogas
   call dealloc_geosol
   call dealloc_pquant
   call dealloc_gquant
   call deallocate_proton_wavefunctions
   call deallocate_wp_wavefunction
   call reset_feszz3_counter

   return

end subroutine deinitmat
