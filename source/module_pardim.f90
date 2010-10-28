module pardim

!======================================================================
!  Dimensions of arrays
!----------------------------------------------------------------------
!  MAXPNT   - maximum number of grid points along the coordinate
!             of quantum particle (the same for gating coordinate)
!  NELST    - number of basis (EVB) electronic states
!  NPRSTMAX - maximum number of basis vibrational states
!             per electronic state for the proton
!  NGASTMAX - maximum number of gating vibrational states
!             per electron-proton vibronic state
!  MAXSTA   - total number of basis vibronic/gating states
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!======================================================================

   implicit none
   public
   save

   integer :: maxatm
   integer :: maxpnt
   integer :: nelst
   integer :: nprstmax
   integer :: ngastmax
   integer :: maxsta

   public :: init_pardim

contains
  
   subroutine init_pardim
      maxatm   = 250
      maxpnt   = 256
      nelst    =   4
      nprstmax =  80
      ngastmax =  80
      maxsta   = nelst*nprstmax*ngastmax
   end subroutine init_pardim

end module pardim
