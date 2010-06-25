module geosol
!======================================================================
!  Geometry of the reaction complex for solvation calculations
!----------------------------------------------------------------------
!  NATSOL - number of atoms
!  IPTSOL - the numbers of atoms in PT interface (Dp-H-Ap)
!  IETSOL - the numbers of atoms in ET interface (De-Ae) (meaningless)
!  LABSOL - atomic numbers of atoms in the complex
!  XYZSOL - Cartesian coordinates of atoms in the complex
!-----------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:36
!  4.1
!  Exp
!  module_geosol.f90,v 4.1 2010/06/25 20:02:36 souda Exp
!  module_geosol.f90,v
!  Revision 4.1  2010/06/25 20:02:36  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2004/01/13 20:02:51  souda
!  Initial PCET-4.0 Release
!
!
!======================================================================
   use pardim

   implicit none
   public
   save

   integer :: natsol
   integer :: iptsol(3), ietsol(2)
   integer, allocatable, dimension(:)   :: labsol
   real(8), allocatable, dimension(:,:) :: xyzsol, chrsol

   !======================================================================
   contains

   subroutine alloc_geosol(ndim)

      implicit none
      integer, intent(in) :: ndim

      real(8) :: zero = 0.d0
      integer :: izero = 0

      if (ndim.le.0) then
         write(*,*) "ERROR in alloc_geosol: Wrong dimension passed..."
         stop
      endif

      if ( .not.allocated(labsol) ) then
         allocate (labsol(ndim))
         labsol = izero
      else
         write(*,*) "ERROR in alloc_geosol: labsol is already allocated..."
         stop
      endif

      if ( .not.allocated(xyzsol) ) then
         allocate (xyzsol(3,ndim))
         xyzsol = zero
      else
         write(*,*) "ERROR in alloc_geosol: xyzsol is already allocated..."
         stop
      endif

      if ( .not.allocated(chrsol) ) then
         allocate (chrsol(nelst,ndim))
         chrsol = zero
      else
         write(*,*) "ERROR in alloc_geosol: chrsol is already allocated..."
         stop
      endif

   end subroutine alloc_geosol


   subroutine dealloc_geosol
      implicit none
      if (allocated(labsol)) deallocate (labsol)
      if (allocated(xyzsol)) deallocate (xyzsol)
      if (allocated(chrsol)) deallocate (chrsol)
   end subroutine dealloc_geosol

end module geosol
