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
!  $Date: 2011-02-24 00:49:24 $
!  $Revision: 5.3 $
!  $Log: not supported by cvs2svn $
!  Revision 5.2  2010/10/28 21:29:36  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
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

   subroutine init_pardim(nelst_)
      integer, intent(in) :: nelst_
      maxatm   = 250
      maxpnt   = 256
      if (nelst_.gt.1) then
         nelst = nelst_
      else
         nelst    =   4
      endif
      nprstmax =  100
      ngastmax =  100
      maxsta   = nelst*nprstmax*ngastmax
   end subroutine init_pardim

end module pardim
