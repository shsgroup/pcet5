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
!  souda
!  2010/06/25 20:02:36
!  4.1
!  Exp
!  module_pardim.f90,v 4.1 2010/06/25 20:02:36 souda Exp
!  module_pardim.f90,v
!  Revision 4.1  2010/06/25 20:02:36  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2004/01/13 20:03:51  souda
!  Initial PCET-4.0 Release
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
  
  subroutine init_pardim
     maxatm   = 250
     maxpnt   = 128
     nelst    =   4
     nprstmax =  40
     ngastmax =  40
     maxsta   = nelst*nprstmax*ngastmax
  end subroutine init_pardim

end module pardim
